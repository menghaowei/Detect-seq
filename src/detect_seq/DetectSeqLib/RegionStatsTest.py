# _*_ coding: UTF-8 _*_

from __future__ import division, print_function
import math

import multiprocessing
import pysam
import logging
import sys

from scipy import stats
from statsmodels.stats import proportion
import statsmodels.stats.multitest as multi
import numpy as np

from DetectSeqLib.CheckAndLoadFiles import load_reference_fasta_as_dict
from DetectSeqLib.FilterRegions import mpmatLine, find_block_info_and_highest_signal

# Version information START ----------------------------------------------------
VERSION_INFO = \
    """
Author: MENG Howard @ 2020-10-30

Version-01:
  2020-10-30 function part for region statistical test

E-Mail: meng_howard@126.com
"""
# Version information END ------------------------------------------------------

# Function List START ----------------------------------------------------------
COPY_RIGHT = \
    """
Belong to MENG Haowei @ YiLab in Peking University

"""


#################################################################################
# FUN
#################################################################################
def check_mpmat_test_state_by_block_info(mpmat_line_obj,
                                         block_site_min_num_cutoff=1,
                                         block_num_ratio_cutoff=0.8,
                                         block_num_ratio_check_min_num=5,
                                         block_num_max_num_cutoff=15):
    """
    INPUT:
        <mpmat_line_obj>
            obj, mpmatLine obj contain BlockInfo annotation

        <block_site_min_num_cutoff>
            int, non-block site num less than this cutoff will return False

        <block_num_ratio_check_min_num>
            int, non-block site num less than this will not check ratio

        <block_num_ratio_cutoff>
            float, when non-block site num >= <block_num_ratio_check_min_num> check ratio.
                If block site num >= int(<block_num_ratio_check_min_num> * <block_num_ratio_cutoff>), return False.

        <block_num_max_num_cutoff>
            int, block site num >= this value, return False

    RETURN:
        <test_state>
            bool, TRUE means region need to be processed Poisson test
                  False means region will be omitted

        <state_reason>
            str, 'NoSignalSite' 'BlockSiteRatio' 'BlockSiteNum'

    """
    # load block number
    block_site_num = mpmat_line_obj.block_state_list.count(True)
    non_block_site_num = mpmat_line_obj.block_state_list.count(False)

    if non_block_site_num < block_site_min_num_cutoff:
        return False, "NoSignalSite"

    if block_site_num >= block_num_ratio_check_min_num:
        block_site_ratio = block_site_num / 1.0 / len(mpmat_line_obj.block_state_list)
        if block_site_ratio >= block_num_ratio_cutoff:
            return False, "BlockSiteRatio"

    if block_site_num >= block_num_max_num_cutoff:
        return False, "BlockSiteNum"

    return True, None


#################################################################################
# FUN
#################################################################################
def back_indel_shift(info_index_list, cur_index):
    """
    INPUT:
        <info_index_list>
            generated from align.cigar tuples

        <cur_index>
            index related to MD tag in BAM file

    RETURN
        <acc_shift>
    """

    # parse soft clip and insertion
    if len(info_index_list) == 0:
        return 0

    acc_shift = 0
    for info_start_index, info_len in info_index_list:
        if info_start_index >= cur_index:
            return acc_shift

        else:
            acc_shift += info_len

    return acc_shift


#################################################################################
# FUN
#################################################################################
def get_align_mismatch_pairs(align):
    """
    INPUT
        <align>
            pysam AlignedSegment object

    RETURN
        <mismatch_pair_list>
            [ref_index, align_index, ref_base, align_base]

            ref_index is the same coordinate with UCSC genome browser

            When NM == 0, return None
    """
    # Hisat-3n aligner NM will be set as 0, but with several conversions.
    # So I have to comment out this part
    # 2022-07-29

    # No mismatch
    # try:
    #     if align.get_tag("NM") == 0:
    #         return None
    # except:
    #     return None

    # parse soft clip, insertion and deletion
    info_index_list = []
    accu_index = 0

    for cigar_type, cigar_len in align.cigartuples:
        if cigar_type == 1 or cigar_type == 4:
            info_index_list.append((accu_index + 1, cigar_len))

        elif cigar_type == 2:
            info_index_list.append((accu_index + 1, -cigar_len))

        accu_index += cigar_len

    # parse MD tag
    mismatch_pair_list = []
    cur_base = ""
    cur_index = 0
    bases = align.get_tag("MD")

    i = 0
    while i < len(bases):
        base = bases[i]

        if base.isdigit():
            cur_base += base
            i += 1

        else:
            if cur_base != "":
                cur_index += int(cur_base)

            cur_base = ""

            if base == "^":
                i += 1
                del_str = ""

                while (bases[i].isalpha()) and (i < len(bases)):
                    del_str += bases[i]
                    i += 1

                cur_index += len(del_str)

            elif base.isalpha():
                cur_index += 1
                ref_base = base
                i += 1

                # add into list
                fix_index = cur_index + back_indel_shift(info_index_list, cur_index)

                if fix_index < len(align.query_sequence):
                    mismatch_pair_list.append([cur_index + align.reference_start, cur_index - 1, ref_base,
                                               align.query_sequence[fix_index - 1]])
                else:
                    return None

    if len(mismatch_pair_list) == 0:
        return None
    else:
        return mismatch_pair_list


#################################################################################
# FUN
#################################################################################
def get_No_MD_align_mismatch_pairs(align, ref_genome_dict):
    """
    INPUT
        <align>
            pysam AlignedSegment object

        <genome_dict>
            key is like chr1, chr2, ...
            value is chromosome sequence

    RETURN
        <mismatch_pair_list>
            [ref_index, align_index, ref_base, align_base]

            ref_index is the same coordinate with UCSC genome browser

            When NM == 0, return None

    """
    # No mismatch

    # Hisat-3n aligner NM will be set as 0, but with several conversions.
    # So I have to comment out this part
    # 2022-07-29

    # try:
    #     if align.get_tag("NM") == 0:
    #         return None
    # except:
    #     return None

    mismatch_pair_list = []
    for align_idx, ref_idx in align.get_aligned_pairs():
        if (align_idx is not None) and (ref_idx is not None):
            align_base = align.query_sequence[align_idx]
            ref_base = ref_genome_dict[align.reference_name][ref_idx]

            if align_base != ref_base:
                mismatch_pair_list.append([
                    ref_idx + 1,
                    align_idx,
                    ref_base,
                    align_base
                ])

    if len(mismatch_pair_list) == 0:
        return None
    else:
        return mismatch_pair_list


#################################################################################
# FUN
#################################################################################
def analyse_align_mut_state(mpmat_info, align_mismatch_pairs, query_mut_type,
                            site_index_dict, snp_index_dict, block_index_dict):
    """
    INPUT:
        <mpmat_info>
            obj, mpmat info

        <align_mismatch_pairs>
            list, format like
                [[2474644, 81, 'G', 'A'], [2474656, 93, 'G', 'A'], [2474679, 116, 'C', 'T']]

        <query_mut_type>
            str, like "GA"

        <site_index_dict>
            dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list

        <snp_index_dict>
            dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list

        <block_index_dict>
            dict, key is pos as str like '2474644', value is site order in mpmat_info.site_index_list

    RETURN:
        <align_mut_state> str
            "M-M-M" means tandem mutation,
            "M-B-M" means mut, block, mut
            "M-S-M" mean mut, SNV, mut
            "N-N-N" means no mut

        <query_mut_count>
            int, query type mutation count on whole alignment read

        <other_mut_count>
            int, other position or other type mutation in query position

        <total_mut_count>
            int, total number of mutation on whole alignment read

        <region_query_mut_count>
            int, total number of reads with query mutation in mpmat region

    VERSION:
        Final edition date: 2022-07-31
    """

    # var init
    total_mut_count = 0
    other_mut_count = 0
    query_mut_count = 0
    region_query_mut_count = 0

    align_mut_state_list = mpmat_info.mut_key_list[:]

    for site_mis_info in align_mismatch_pairs:
        total_mut_count += 1

        site_index, site_align_pos, from_base, to_base = site_mis_info

        # check block
        block_site_order = block_index_dict.get(str(site_index))

        # check SNP
        snp_site_order = snp_index_dict.get(str(site_index))

        if block_site_order is not None:
            continue

        if snp_site_order is not None:
            continue

        # not in block info list
        site_order = site_index_dict.get(str(site_index))

        # count mutation num
        if (from_base == query_mut_type[0]) and (to_base == query_mut_type[1]):
            query_mut_count += 1

            # check mutation site only in mpmat region
            if site_order is not None:
                region_query_mut_count += 1
                align_mut_state_list[site_order] = "M"

        else:
            other_mut_count += 1

    # align mut state
    align_mut_state = "-".join(align_mut_state_list)

    return align_mut_state, query_mut_count, other_mut_count, total_mut_count, region_query_mut_count


#################################################################################
# FUN
#################################################################################
# make a function
def get_mpmat_region_count(
        in_bam_obj,
        mpmat_info,
        ref_genome_dict,
        query_mut_type,
        query_mut_min_cutoff=1,
        query_mut_max_cutoff=16,
        total_mut_max_cutoff=20,
        other_mut_max_cutoff=16
):
    """
    INPUT:
        <in_bam_obj>
            obj, pysam.AlignmentFile

        <mpmat_info>
            obj, mpmatLine obj

        <ref_genome_dict>
            dict, key is chr_name, value is reference sequence

        <query_mut_min_cutoff>
            int, if mutation number >= query_mut_min_cutoff in the mpmat region,
                will be counted as 'region_mut_count'.

        <query_mut_max_cutoff>
            int, larger than this will be marked.

        <total_mut_max_cutoff>
            int, larger than this will be marked.

        <other_mut_max_cutoff>
            int, larger than this will be marked.

    RETURN:
        <align_count_dict>
            dict, key and value like:
                {
                    "all_align_count": 0,
                    "region_mut_count": 0,
                    "region_non_mut_count": 0,
                    "total_high_mismatch_count": 0,
                    "other_high_mismatch_count": 0,
                    "query_high_mismatch_count": 0,
                    "all_filter_count": 0
                }

        <align_mut_tandem_dict>
            dict, key is  <align_mut_state> (refer to FUN analyse_align_mut_state)
                "M-M-M" means tandem mutation,
                "M-B-M" means mut, block, mut
                "M-S-M" means mut, SNP, mut
                "N-N-N" means no mut

        <align_mut_count_dict>
            dict, key is mutation number like 0,1,2,3,4.... value is align reads count

        <mpmat_info>
            obj, mpmatLine obj, change <mut_key_list> and <mut_key> according to mpmat block info

    VERSION:
        Final edition date: 2022-07-31
    """
    # ---------------------------------------------------------->>>>>>>>>>
    # init var
    # ---------------------------------------------------------->>>>>>>>>>
    # site index dict
    site_index_dict = {}

    # site SNP dict
    snp_site_index_dict = {}

    # site block dict
    block_site_index_dict = {}

    # fix mpmat_info
    for index, site_index in enumerate(mpmat_info.site_index_list):
        site_index_dict[site_index.split("_")[1]] = index

        if mpmat_info.block_info_list[index]:
            block_site_index_dict[site_index.split("_")[1]] = index

        if mpmat_info.SNP_ann_list[index]:
            snp_site_index_dict[site_index.split("_")[1]] = index

    # make non mut key
    non_mut_key = mpmat_info.mut_key
    align_mut_tandem_dict = {non_mut_key: 0}

    # define count dict
    align_count_dict = {
        "all_align_count": 0,
        "region_mut_count": 0,
        "region_non_mut_count": 0,
        "total_high_mismatch_count": 0,
        "other_high_mismatch_count": 0,
        "query_high_mismatch_count": 0,
        "all_filter_count": 0
    }

    # mutation count dict
    align_mut_count_dict = {}
    for mut_num in range(mpmat_info.site_num + 1):
        align_mut_count_dict[mut_num] = 0

    # ---------------------------------------------------------->>>>>>>>>>
    # iter align info
    # ---------------------------------------------------------->>>>>>>>>>
    for align in in_bam_obj.fetch(reference=mpmat_info.chr_name,
                                  start=mpmat_info.chr_start - 1,
                                  end=mpmat_info.chr_end + 1):

        # count total
        align_count_dict["all_align_count"] += 1

        # make sure MD state
        MD_state = True
        try:
            MD_str_tag = align.get_tag("MD")
        except:
            MD_state = False

        # get mismatch pairs
        if MD_state:
            align_mismatch_pairs = get_align_mismatch_pairs(align)
        else:
            align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)

        # analysis of mismatch pairs
        if align_mismatch_pairs is None:
            align_count_dict["region_non_mut_count"] += 1
            align_mut_count_dict[0] += 1
            align_mut_tandem_dict[non_mut_key] += 1

        else:
            align_mut_analyse_res = analyse_align_mut_state(mpmat_info=mpmat_info,
                                                            align_mismatch_pairs=align_mismatch_pairs,
                                                            query_mut_type=query_mut_type,
                                                            site_index_dict=site_index_dict,
                                                            snp_index_dict=snp_site_index_dict,
                                                            block_index_dict=block_site_index_dict)

            if align_mut_tandem_dict.get(align_mut_analyse_res[0]) is None:
                align_mut_tandem_dict[align_mut_analyse_res[0]] = 1
            else:
                align_mut_tandem_dict[align_mut_analyse_res[0]] += 1

            # load mut num
            query_mut_count, other_mut_count, total_mut_count, region_query_mut_count = align_mut_analyse_res[1:]
            align_mut_count_dict[region_query_mut_count] += 1

            # filter high mismatch reads
            if total_mut_count >= total_mut_max_cutoff:
                align_count_dict["total_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            elif other_mut_count >= other_mut_max_cutoff:
                align_count_dict["other_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            elif query_mut_count >= query_mut_max_cutoff:
                align_count_dict["query_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            else:
                if region_query_mut_count == 0:
                    align_count_dict["region_non_mut_count"] += 1

                elif region_query_mut_count >= query_mut_min_cutoff:
                    align_count_dict["region_mut_count"] += 1

    return align_count_dict, align_mut_tandem_dict, align_mut_count_dict


#################################################################################
# Statistics Function part
#################################################################################
def _zstat_generic_parse(value, std_diff, alternative):
    """
    INPUT:
        <value>
            Stats

        <std_diff>
            Std

        <alternative>
            Alternative string

    RETURN:
        z-score, pvalue

    HELP:
        Copied and fixed from statsmodels.stats.weight stats
    """

    zstat = value / std_diff

    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = stats.norm.sf(np.abs(zstat)) * 2

    elif alternative in ['larger', 'l']:
        pvalue = stats.norm.sf(zstat)

    elif alternative in ['smaller', 's']:
        pvalue = stats.norm.cdf(zstat)

    else:
        raise ValueError('invalid alternative')

    return zstat, pvalue


def poisson_test(lambda_1, lambda_2, method='sqrt', alternative='two-sided'):
    """
    INPUT:

        <lambda_1>, <lambda_2>
            Poisson parameter in case1 and case2

        <method>
            Support method:
                'wald': method W1A, wald test, variance based on separate estimates
                'score': method W2A, score test, variance based on estimate under Null
                'sqrt': W5A, based on variance stabilizing square root transformation

                'exact-cond': exact conditional test based on binomial distribution
                'cond-midp': midpoint-pvalue of exact conditional test

        <alternative>
            'two-sided'
                means H0: lambda_1 = lambda_2

            'larger'
                means H0: lambda_1 > lambda_2

            'smaller'
                means H0: lambda_1 < lambda_2

    HELP:

        Reference paper:
            Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
            Biometrical Journal 50 (2008) 2, 2008

        Raw code is copied from html and then fixed realted to DETECT-Seq project requirement
            https://stackoverflow.com/questions/33944914/implementation-of-e-test-for-poisson-in-python
    """

    # calculate stat
    dist = ""
    stat = 0

    if method in ['score']:
        stat = (lambda_1 - lambda_2) / np.sqrt((lambda_1 + lambda_2))
        dist = 'normal'

    elif method in ['wald']:
        stat = (lambda_1 - lambda_2) / np.sqrt((lambda_1 + lambda_2))
        dist = 'normal'

    elif method in ['sqrt']:
        stat = 2 * (np.sqrt(lambda_1 + 3 / 1.0 / 8) - np.sqrt((lambda_2 + 3 / 1.0 / 8)))
        stat /= np.sqrt(2)
        dist = 'normal'

    elif method in ['exact-cond', 'cond-midp']:
        lambda_total = lambda_1 + lambda_2
        stat = None
        pvalue = proportion.binom_test(lambda_1, lambda_total, prop=0.5, alternative=alternative)
        return stat, pvalue

    # return part
    if dist == 'normal':
        return _zstat_generic_parse(stat, 1, alternative)


def run_mpmat_poisson_test(
        mpmat_filename,
        out_poisson_filename,
        ctrl_bam_filename,
        treat_bam_filename,
        ref_genome_fa_filename,
        select_chr_name_list,
        scale_factor_dict=None,
        normalize_scale_factor_dict=None,
        genome_bg_dict=None,
        lambda_bg_method="ctrl_max",
        poisson_method="mutation",
        region_block_mut_num_cutoff=2,
        reads_query_mut_min_cutoff=1,
        reads_query_mut_max_cutoff=16,
        reads_total_mut_max_cutoff=20,
        reads_other_mut_max_cutoff=16,
        log_verbose=3,
        mpmat_filter_col_idx=-1,
        mpmat_block_col_idx=-1
):
    """

    RETURN:
            0, means everything is okay.
    """
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    process_cmd_str = """run_mpmat_poisson_test step with params:
        mpmat_filename={mpmat_filename}
        out_poisson_filename={out_poisson_filename}
        ctrl_bam_filename={ctrl_bam_filename}
        treat_bam_filename={treat_bam_filename}
        ref_genome_fa_filename={ref_genome_fa_filename}
        select_chr_name_list={select_chr_name_list}
        scale_factor_dict={scale_factor_dict}
        normalize_scale_factor_dict={normalize_scale_factor_dict}
        genome_bg_dict={genome_bg_dict}
        lambda_bg_method={lambda_bg_method}
        poisson_method={poisson_method}
        region_block_mut_num_cutoff={region_block_mut_num_cutoff}
        reads_query_mut_min_cutoff={reads_query_mut_min_cutoff}
        reads_query_mut_max_cutoff={reads_query_mut_max_cutoff}
        reads_total_mut_max_cutoff={reads_total_mut_max_cutoff}
        reads_other_mut_max_cutoff={reads_other_mut_max_cutoff}
        log_verbose={log_verbose}""".format(
        mpmat_filename=mpmat_filename,
        out_poisson_filename=out_poisson_filename,
        ctrl_bam_filename=ctrl_bam_filename,
        treat_bam_filename=treat_bam_filename,
        ref_genome_fa_filename=ref_genome_fa_filename,
        select_chr_name_list=select_chr_name_list,
        scale_factor_dict=scale_factor_dict,
        normalize_scale_factor_dict=normalize_scale_factor_dict,
        genome_bg_dict=genome_bg_dict,
        lambda_bg_method=lambda_bg_method,
        poisson_method=poisson_method,
        region_block_mut_num_cutoff=region_block_mut_num_cutoff,
        reads_query_mut_min_cutoff=reads_query_mut_min_cutoff,
        reads_query_mut_max_cutoff=reads_query_mut_max_cutoff,
        reads_total_mut_max_cutoff=reads_total_mut_max_cutoff,
        reads_other_mut_max_cutoff=reads_other_mut_max_cutoff,
        log_verbose=log_verbose)

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # open files
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    try:
        chr_mpmat_file = open(mpmat_filename, "r")
        ctrl_bam = pysam.AlignmentFile(ctrl_bam_filename, "rb")
        treat_bam = pysam.AlignmentFile(treat_bam_filename, "rb")

        if out_poisson_filename == "stdout":
            out_file = sys.stdout
        else:
            out_file = open(out_poisson_filename, "w")

    except:
        raise IOError("Open files error! Please check INPUT and OUTPUT filenames!")

    logging.info("Starting to run Poisson test on \n\t%s" % mpmat_filename)

    logging.debug(process_cmd_str)

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # load FASTA genome
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    ref_genome_dict = load_reference_fasta_as_dict(ref_fasta_path=ref_genome_fa_filename,
                                                   ref_name_list=select_chr_name_list)

    # get pysam Fasta obj
    ref_genome_fa_obj = pysam.FastaFile(ref_genome_fa_filename)

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # iter to run test
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    if mpmat_filter_col_idx == -1:
        parse_mpmat_filter_col_idx = None
    else:
        parse_mpmat_filter_col_idx = mpmat_filter_col_idx

    if mpmat_block_col_idx == -1:
        parse_mpmat_block_col_idx = None
    else:
        parse_mpmat_block_col_idx = mpmat_block_col_idx

    for line_index, line in enumerate(chr_mpmat_file):
        line_list = line.strip().split("\t")
        mpmat_info = mpmatLine(line_list,
                               block_info_index=parse_mpmat_block_col_idx,
                               filter_info_index=parse_mpmat_filter_col_idx)

        # make dict
        mpmat_info_dict = {
            "ctrl_all_count": None,
            "treat_all_count": None,
            "ctrl_mut_count": None,
            "treat_mut_count": None,
            "ctrl_all_count.norm": None,
            "treat_all_count.norm": None,
            "ctrl_mut_count.norm": None,
            "treat_mut_count.norm": None,
            "log2FC_all_count": None,
            "log2FC_mut_count": None,
            "region_count_info": None,
            "test_state": None,
            "pvalue": None,
            "region_site_index": None,
            "region_site_num": None,
            "region_block_site_num": None,
            "region_mut_site_num": None,
            "region_block_state": None,
            "region_highest_site_index": None,
            "region_highest_site_mut_num": None,
            "region_highest_site_cover_num": None,
            "region_highest_site_mut_ratio": None,
        }

        # chr name
        mpmat_chr_name = mpmat_info.chr_name

        # log
        if log_verbose == 0:
            run_report_num = 10000
        elif log_verbose == 1:
            run_report_num = 1000
        else:
            run_report_num = 100

        if line_index % run_report_num == 0:
            logging.info("Running Poisson test on %s, processed line number %s" % (mpmat_chr_name, line_index + 1))

        # check mpmat region state [old-block-methods]
        # if parse_mpmat_block_col_idx is not None:
        #     mpmat_check_res = check_mpmat_test_state_by_block_info(
        #         mpmat_info,
        #         block_site_min_num_cutoff=region_block_site_min_num_cutoff,
        #         block_num_ratio_cutoff=region_block_num_ratio_cutoff,
        #         block_num_ratio_check_min_num=region_block_num_ratio_check_min_num,
        #         block_num_max_num_cutoff=region_block_num_max_num_cutoff
        #     )
        # else:
        #     mpmat_check_res = [True]

        # ---------------------------------------------------------->>>>>>>>>>
        # new-add-step.a hard block and find highest signal
        # ---------------------------------------------------------->>>>>>>>>>
        mpmat_info_block = find_block_info_and_highest_signal(mpmat_info,
                                                              ctrl_bam_obj=ctrl_bam,
                                                              treat_bam_obj=treat_bam,
                                                              ref_genome_obj=ref_genome_fa_obj,
                                                              block_mut_num_cutoff=2)

        # ---------------------------------------------------------->>>>>>>>>>
        # step1. get align count info
        # ---------------------------------------------------------->>>>>>>>>>
        # pipeline
        ctrl_count_res = get_mpmat_region_count(in_bam_obj=ctrl_bam,
                                                mpmat_info=mpmat_info_block,
                                                ref_genome_dict=ref_genome_dict,
                                                query_mut_type=mpmat_info.mut_type,
                                                query_mut_min_cutoff=reads_query_mut_min_cutoff,
                                                query_mut_max_cutoff=reads_query_mut_max_cutoff,
                                                total_mut_max_cutoff=reads_total_mut_max_cutoff,
                                                other_mut_max_cutoff=reads_other_mut_max_cutoff)

        treat_count_res = get_mpmat_region_count(in_bam_obj=treat_bam,
                                                 mpmat_info=mpmat_info_block,
                                                 ref_genome_dict=ref_genome_dict,
                                                 query_mut_type=mpmat_info.mut_type,
                                                 query_mut_min_cutoff=reads_query_mut_min_cutoff,
                                                 query_mut_max_cutoff=reads_query_mut_max_cutoff,
                                                 total_mut_max_cutoff=reads_total_mut_max_cutoff,
                                                 other_mut_max_cutoff=reads_other_mut_max_cutoff)

        ctrl_mpmat_count_dict = ctrl_count_res[0]
        treat_mpmat_count_dict = treat_count_res[0]

        # ---------------------------------------------------------->>>>>>>>>>
        # step2. align count normalization
        # ---------------------------------------------------------->>>>>>>>>>
        # region mut lambda
        # try:
        ctrl_mpmat_lambda = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 \
                            / scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
        treat_mpmat_lambda = treat_mpmat_count_dict["region_mut_count"] / 1.0 \
                             / scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]

        # all reads lambda
        ctrl_mpmat_lambda_all = ctrl_mpmat_count_dict["all_align_count"] / 1.0 \
                                / scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
        treat_mpmat_lambda_all = treat_mpmat_count_dict["all_align_count"] / 1.0 \
                                 / scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]

        # ---------------------------------------------------------->>>>>>>>>>
        # step3. Poisson test
        # ---------------------------------------------------------->>>>>>>>>>
        # global vars
        treat_lambda = 0
        bg_lambda = 0

        # Poisson test
        if lambda_bg_method == "raw":
            if poisson_method == "mutation":
                bg_lambda = ctrl_mpmat_count_dict["region_mut_count"]
                treat_lambda = treat_mpmat_count_dict["region_mut_count"]

            elif poisson_method == "all":
                bg_lambda = ctrl_mpmat_count_dict["all_align_count"]
                treat_lambda = treat_mpmat_count_dict["all_align_count"]

            else:
                logging.error("Set wrong Poisson method!")
                raise ValueError("Set wrong Poisson method!")

        else:
            if poisson_method == "mutation":
                # make mutation lambda list
                ctrl_bg_lambda_list = [
                    genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"],
                    genome_bg_dict["ctrl"]["scale_mut_bg"][mpmat_chr_name],
                    ctrl_mpmat_lambda
                ]

                treat_bg_lambda_list = [
                    genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"],
                    genome_bg_dict["treat"]["scale_mut_bg"][mpmat_chr_name]
                ]

                bg_lambda_list = [
                    genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"],
                    genome_bg_dict["ctrl"]["scale_mut_bg"][mpmat_chr_name],
                    ctrl_mpmat_lambda,
                    genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"],
                    genome_bg_dict["treat"]["scale_mut_bg"][mpmat_chr_name]
                ]

                # set treat lambda value
                treat_lambda = treat_mpmat_lambda

            elif poisson_method == "all":
                # make mutation lambda list
                ctrl_bg_lambda_list = [
                    genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"],
                    genome_bg_dict["ctrl"]["scale_all_bg"][mpmat_chr_name],
                    ctrl_mpmat_lambda_all
                ]

                treat_bg_lambda_list = [
                    genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"],
                    genome_bg_dict["treat"]["scale_all_bg"][mpmat_chr_name]
                ]

                bg_lambda_list = [
                    genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"],
                    genome_bg_dict["ctrl"]["scale_all_bg"][mpmat_chr_name],
                    ctrl_mpmat_lambda_all,
                    genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"],
                    genome_bg_dict["treat"]["scale_all_bg"][mpmat_chr_name]
                ]

                # set treat lambda value
                treat_lambda = treat_mpmat_lambda_all

            else:
                logging.error("Set wrong Poisson method!")
                raise ValueError("Set wrong Poisson method!")

            # select bg lambda
            if lambda_bg_method == "ctrl_max":
                bg_lambda = max(ctrl_bg_lambda_list)

            elif lambda_bg_method == "treat_max":
                bg_lambda = max(treat_bg_lambda_list)

            elif lambda_bg_method == "max":
                bg_lambda = max(bg_lambda_list)

            elif lambda_bg_method == "raw":
                pass

            else:
                logging.error("Set wrong lambda method!")

        # Poisson test
        mut_pvalue = poisson_test(treat_lambda, bg_lambda, alternative="larger", method="sqrt")[1]

        # record state
        state_test = "TestOK"

        # ---------------------------------------------------------->>>>>>>>>>
        # step4. record signal and pvalue
        # ---------------------------------------------------------->>>>>>>>>>
        # make count string
        mut_num_list = []
        ctrl_count_list = []
        treat_count_list = []

        for mut_num in range(mpmat_info.site_num + 1):
            mut_num_list.append(mut_num)
            ctrl_count_list.append(ctrl_count_res[2][mut_num])
            treat_count_list.append(treat_count_res[2][mut_num])

        count_str_list = [mut_num_list, ctrl_count_list, treat_count_list]
        count_str = " ".join([",".join(map(str, x)) for x in count_str_list])
        mpmat_info_dict["region_count_info"] = count_str

        # calculate normalized signal
        mpmat_ctrl_mut_count_norm = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 / normalize_scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
        mpmat_treat_mut_count_norm = treat_mpmat_count_dict["region_mut_count"] / 1.0 / normalize_scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]

        mpmat_ctrl_all_count_norm = ctrl_mpmat_count_dict["all_align_count"] / 1.0 / normalize_scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
        mpmat_treat_all_count_norm = treat_mpmat_count_dict["all_align_count"] / 1.0 / normalize_scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]

        # All -> calculate log2 fold change fix

        if mpmat_ctrl_all_count_norm != 0:
            all_count_FC = mpmat_treat_all_count_norm / 1.0 / mpmat_ctrl_all_count_norm
        else:
            all_count_FC = mpmat_treat_all_count_norm / 1.0 / genome_bg_dict["ctrl"]["norm_scale_all_bg"][mpmat_chr_name]

        if all_count_FC > 0:
            log2FC_all_count = math.log(all_count_FC, 2)
        else:
            log2FC_all_count = "NA"

        # Mut -> calculate log2 fold change fix
        if mpmat_ctrl_mut_count_norm != 0:
            mut_count_FC = mpmat_treat_mut_count_norm / 1.0 / mpmat_ctrl_mut_count_norm
        else:
            mut_count_FC = mpmat_treat_mut_count_norm / 1.0 / genome_bg_dict["ctrl"]["norm_scale_mut_bg"][mpmat_chr_name]

        if mut_count_FC > 0:
            log2FC_mut_count = math.log(mut_count_FC, 2)
        else:
            log2FC_mut_count = "NA"

        # add all info into dict
        mpmat_info_dict["ctrl_all_count"] = ctrl_mpmat_count_dict["all_align_count"]
        mpmat_info_dict["treat_all_count"] = treat_mpmat_count_dict["all_align_count"]

        mpmat_info_dict["ctrl_mut_count"] = ctrl_mpmat_count_dict["region_mut_count"]
        mpmat_info_dict["treat_mut_count"] = treat_mpmat_count_dict["region_mut_count"]

        mpmat_info_dict["ctrl_all_count.norm"] = mpmat_ctrl_all_count_norm
        mpmat_info_dict["treat_all_count.norm"] = mpmat_treat_all_count_norm

        mpmat_info_dict["ctrl_mut_count.norm"] = mpmat_ctrl_mut_count_norm
        mpmat_info_dict["treat_mut_count.norm"] = mpmat_treat_mut_count_norm

        mpmat_info_dict["log2FC_all_count"] = log2FC_all_count
        mpmat_info_dict["log2FC_mut_count"] = log2FC_mut_count

        # add highest info
        mpmat_info_dict["region_site_index"] = ",".join(mpmat_info_block.site_index_list)
        mpmat_info_dict["region_site_num"] = mpmat_info_block.site_num
        mpmat_info_dict["region_block_site_num"] = mpmat_info_block.block_site_num
        mpmat_info_dict["region_mut_site_num"] = mpmat_info_block.site_num - mpmat_info_block.block_site_num
        mpmat_info_dict["region_block_state"] = mpmat_info_block.mut_key
        mpmat_info_dict["region_highest_site_index"] = mpmat_info_block.highest_site_dict["full_site_index"]
        mpmat_info_dict["region_highest_site_mut_num"] = mpmat_info_block.highest_site_dict["mut_num"]
        mpmat_info_dict["region_highest_site_cover_num"] = mpmat_info_block.highest_site_dict["total"]
        mpmat_info_dict["region_highest_site_mut_ratio"] = mpmat_info_block.highest_site_dict["mut_ratio"]

        # record pvalue
        mpmat_info_dict["test_state"] = state_test
        mpmat_info_dict["pvalue"] = mut_pvalue

        # ---------------------------------------------------------->>>>>>>>>>
        # step5. output part
        # ---------------------------------------------------------->>>>>>>>>>
        info_list = line_list[:3]

        # region index
        info_list.append("%s_%s_%s" % (line_list[0], line_list[1], line_list[2]))

        info_list += [
            mpmat_info_dict["region_site_num"],
            mpmat_info_dict["region_block_site_num"],
            mpmat_info_dict["region_mut_site_num"],
            mpmat_info_dict["region_site_index"],
            mpmat_info_dict["region_block_state"],
            mpmat_info_dict["region_highest_site_index"],
            mpmat_info_dict["region_highest_site_mut_num"],
            mpmat_info_dict["region_highest_site_cover_num"],
            mpmat_info_dict["region_highest_site_mut_ratio"],
            mpmat_info_dict["ctrl_all_count"],
            mpmat_info_dict["treat_all_count"],
            mpmat_info_dict["ctrl_mut_count"],
            mpmat_info_dict["treat_mut_count"],
            mpmat_info_dict["ctrl_all_count.norm"],
            mpmat_info_dict["treat_all_count.norm"],
            mpmat_info_dict["ctrl_mut_count.norm"],
            mpmat_info_dict["treat_mut_count.norm"],
            mpmat_info_dict["region_count_info"],
            mpmat_info_dict["log2FC_all_count"],
            mpmat_info_dict["log2FC_mut_count"],
            mpmat_info_dict["test_state"],
            mpmat_info_dict["pvalue"]
        ]

        # write into file
        info_list_str = "\t".join(map(str, info_list))
        out_file.write(info_list_str + "\n")

    # close files
    chr_mpmat_file.close()
    ctrl_bam.close()
    treat_bam.close()

    if out_poisson_filename != "stdout":
        out_file.close()

    # logging
    logging.info("Run Poisson test successfully on .mpmat file \n\t%s" % mpmat_filename)

    return 0


def multi_run_mpmat_poisson_test(
        mpmat_block_split_dict,
        ctrl_bam_filename,
        treat_bam_filename,
        ref_genome_fa_filename,
        scale_factor_dict=None,
        normalize_scale_factor_dict=None,
        genome_bg_dict=None,
        lambda_bg_method="ctrl_max",
        poisson_method="mutation",
        log_verbose=3,
        thread=1,
        **region_filter_args
):
    """
    INPUT:

        <**region_filter_args>
            args list:
                <region_block_mut_num_cutoff>
                <reads_query_mut_min_cutoff>
                <reads_query_mut_max_cutoff>
                <reads_total_mut_max_cutoff>
                <reads_other_mut_max_cutoff>
                <mpmat_filter_info_col_index>
                <mpmat_block_info_col_index>

    RETURN
        <out_poisson_filename_dict>
            dict, contains output mpmat files with block information, and format like

            {
                'chr1': 'chr1.eqdA2zJkHCm0YW1o.AddBlockInfo',
                'chr19': 'chr19.ZIrU6xC2DFcBayE1.AddBlockInfo',
                'chr20': 'chr20.nK2goEB6xh9MzXpD.AddBlockInfo',
                'chr_name_order': ['chr1', 'chr19', 'chr20']
            }

    """
    # ------------------------------------------------------------>>>>>>>>>>
    # log setting
    # ------------------------------------------------------------>>>>>>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ------------------------------------------------------------>>>>>>>>>>
    # params setting
    # ------------------------------------------------------------>>>>>>>>>>
    run_params_dict = {
        "region_block_mut_num_cutoff": 2,
        "reads_query_mut_min_cutoff": 1,
        "reads_query_mut_max_cutoff": 16,
        "reads_total_mut_max_cutoff": 20,
        "reads_other_mut_max_cutoff": 16,
        "mpmat_filter_info_col_index": -1,
        "mpmat_block_info_col_index": -1
    }

    for args_key in region_filter_args:
        if region_filter_args.get(args_key) is not None:
            run_params_dict[args_key] = region_filter_args[args_key]

    # ------------------------------------------------------------>>>>>>>>>>
    # check params
    # ------------------------------------------------------------>>>>>>>>>>
    if scale_factor_dict is None:
        raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'scale_factor_dict' !")

    if normalize_scale_factor_dict is None:
        raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'normalize_scale_factor_dict' !")

    if genome_bg_dict is None:
        raise IOError("FUN <multi_run_mpmat_poisson_test> need set 'genome_bg_dict' !")

    # ------------------------------------------------------------>>>>>>>>>>
    # get chr name order
    # ------------------------------------------------------------>>>>>>>>>>
    chr_name_order_list = mpmat_block_split_dict["chr_name_order"]
    out_poisson_filename_dict = {"chr_name_order": chr_name_order_list}

    # ------------------------------------------------------------>>>>>>>>>>
    # run part
    # ------------------------------------------------------------>>>>>>>>>>
    logging.info("-" * 80)
    logging.info("Starting to run Poisson test...")

    pool = multiprocessing.Pool(processes=thread)

    run_return_info_list = []

    for chr_name in chr_name_order_list:
        # mpmat input and output
        chr_mpmat_old_filename = mpmat_block_split_dict[chr_name]
        chr_poisson_new_filename = chr_mpmat_old_filename + "." + "PoissonResult"
        out_poisson_filename_dict[chr_name] = chr_poisson_new_filename

        run_return_info_list.append(
            pool.apply_async(
                func=run_mpmat_poisson_test,
                args=(
                    chr_mpmat_old_filename,
                    chr_poisson_new_filename,
                    ctrl_bam_filename,
                    treat_bam_filename,
                    ref_genome_fa_filename,
                    [chr_name],
                    scale_factor_dict,
                    normalize_scale_factor_dict,
                    genome_bg_dict,
                    lambda_bg_method,
                    poisson_method,
                    run_params_dict["region_block_mut_num_cutoff"],
                    run_params_dict["reads_query_mut_min_cutoff"],
                    run_params_dict["reads_query_mut_max_cutoff"],
                    run_params_dict["reads_total_mut_max_cutoff"],
                    run_params_dict["reads_other_mut_max_cutoff"],
                    log_verbose,
                    run_params_dict["mpmat_filter_info_col_index"],
                    run_params_dict["mpmat_block_info_col_index"],
                )
            )
        )

    pool.close()
    pool.join()

    # check run state
    final_run_state = 0
    for index, res in enumerate(run_return_info_list):
        run_state = res.get()
        if run_state != 0:
            logging.error("Poisson test error occur with %s!" % chr_name_order_list[index])
            if final_run_state == 0:
                final_run_state = 1

    if final_run_state == 0:
        logging.info("Calculation of Poisson test result. Done!")
        return out_poisson_filename_dict

    else:
        logging.error("Something wrong with Poisson test step!")
        raise RuntimeError()


####################################################################################################
# FDR with BH method
####################################################################################################
def make_qvalue_with_BH_method(pval_list):
    """
    INPUT:
        <pval_list>
            list, may contain 'NA' value

    RETURN
        <fdr_list>
            list, return FDR list with BH method.

    """

    # init vars
    raw_pval_index_dict = {}
    rm_NA_pval_list = []
    run_index = 0

    for index, pval in enumerate(pval_list):
        if pval != "NA":
            rm_NA_pval_list.append(pval)
            raw_pval_index_dict[run_index] = index
            run_index += 1

    FDR_qvalue = multi.multipletests(np.array(rm_NA_pval_list), alpha=0.05, method="fdr_bh", is_sorted=False)
    FDR_qvalue_vec = FDR_qvalue[1]

    return_fdr_list = ["NA"] * len(pval_list)

    for index, fdr in enumerate(FDR_qvalue_vec):
        raw_index = raw_pval_index_dict[index]
        return_fdr_list[raw_index] = fdr

    return return_fdr_list
