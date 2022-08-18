# _*_ coding: UTF-8 _*_
import gzip
import logging
import sys
import multiprocessing

from scipy import stats
from math import log

# Version information START ----------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard @ 2020-10-29

Version-01:
  2020-10-29 function part for check mpmat line

E-Mail: meng_howard@126.com

"""
# Version information END ------------------------------------------------------

# Function List START ----------------------------------------------------------
COPY_RIGHT = \
"""
Belong to MENG Haowei @ YiLab in Peking University

"""
# Function List END ------------------------------------------------------------

# Other Info START -------------------------------------------------------------
# Final edit time: 2020-10-30
# Other Info END ---------------------------------------------------------------


#################################################################################
# class
#################################################################################
class mpmatLine(object):
    """
    INPUT:
        <mpmat_line_list>
            list, info like
            [
                'chr1',
                '629627',
                '629632',
                '4',
                '4',
                '0',
                'chr1_629627_CT,chr1_629628_CT,chr1_629631_CT,chr1_629632_CT',
                '3,5,3,4',
                '37,38,37,38',
                '0.08108,0.13158,0.08108,0.10526',
                'False,False,False,False',
                '0,0,0,0',
                'Pass,Pass,Pass,Pass'
            ]

        <filter_info_index>
            int, col index describe filter info, default=None, start at 1

        <block_info_index>
            int, col index describe block info, default=None, start at 1

    RETURN:
        <mpmatLine obj>
    """

    def __init__(self, mpmat_line_list, filter_info_index=None, block_info_index=None):
        # standard .mpmat file
        self.chr_name = mpmat_line_list[0]
        self.chr_start = int(mpmat_line_list[1])
        self.chr_end = int(mpmat_line_list[2])

        self.site_num = int(mpmat_line_list[3])
        self.mut_site_num = int(mpmat_line_list[4])
        self.SNP_site_num = int(mpmat_line_list[5])

        self.site_index_list = mpmat_line_list[6].split(",")
        self.mut_count_list = list(map(int, mpmat_line_list[7].split(",")))
        self.cover_count_list = list(map(int, mpmat_line_list[8].split(",")))
        self.mut_ratio_list = list(map(float, mpmat_line_list[9].split(",")))
        self.SNP_ann_list = list(map(eval, mpmat_line_list[10].split(",")))
        self.tandem_info_list = mpmat_line_list[11].split(",")

        # region str
        self.region = "%s:%s-%s" % (mpmat_line_list[0], mpmat_line_list[1], mpmat_line_list[2])

        # get mut type
        self.mut_type = mpmat_line_list[6].split(",")[0].split("_")[-1]

        # mutation key
        mut_key_list = ["N"] * self.site_num
        for index, site_index in enumerate(self.site_index_list):
            if self.SNP_ann_list[index]:
                mut_key_list[index] = "S"

        self.mut_key = "-".join(mut_key_list)
        self.mut_key_list = mut_key_list

        # load filter
        if filter_info_index is not None:
            try:
                self.filter_state_list = mpmat_line_list[filter_info_index - 1].split(",")
            except:
                raise IOError("Parsing error occur at <filter_info_index>")

        if block_info_index is not None:
            try:
                self.block_info_list = list(map(eval, mpmat_line_list[block_info_index - 1].split(",")))
                self.block_site_num = self.block_info_list.count(True)
            except:
                raise IOError("Parsing error occur at <block_info_index>")

        # [old-version-block]
        # if block_info_index is not None:
        #     try:
        #         block_info_list =
        #         self.block_name_list = block_info_list[0].split(",")
        #         self.ctrl_hf_state_list = list(map(eval, block_info_list[1].split(",")))
        #         self.ctrl_hf_reason_list = block_info_list[2].split(",")
        #         self.ctrl_binomial_pval_list = list(map(int, block_info_list[3].split(",")))
        #         self.treat_hf_state_list = list(map(eval, block_info_list[4].split(",")))
        #         self.treat_hf_reason_list = block_info_list[5].split(",")
        #         self.block_state_list = list(map(eval, block_info_list[6].split(",")))
        #         self.block_reason_list = block_info_list[7].split(",")
        #     except:
        #         raise IOError("Parsing error occur at <block_info_index>")


#################################################################################
# class
#################################################################################
class siteIndex(object):
    """
    INPUT
        <site_index>
            str, like chr1_10000_CT
        
    RETURN
        siteIndex obj
    """

    def __init__(self, _site_index):
        site_index_split = _site_index.split("_")

        self.chr_name = site_index_split[0]
        self.site_pos = int(site_index_split[1])
        self.site_index = site_index_split[0] + "_" + site_index_split[1]
        self.site_index_raw = _site_index

        if len(site_index_split) > 2:
            self.mut_type = site_index_split[2]


#################################################################################
# class
################################################################################# 
class bmatLine(object):
    """
    INPUT:
        <bmat_line_list>
            list, info like
            [
               'chr1', 
               '10013', 
               'T', 
               '0', 
               '0', 
               '0', 
               '6', 
               '0', 
               '0', 
               '0', 
               '.', 
               '.', 
               '.', 
               '0'
            ]

    RETURN:
        A bmatLine obj
    """

    def __init__(self, bmat_line_list):
        self.chr_name = bmat_line_list[0]
        self.chr_index = int(bmat_line_list[1])
        self.ref_base = bmat_line_list[2]

        self.count_del = int(bmat_line_list[7])
        self.count_insert = int(bmat_line_list[8])
        self.count_ambis = int(bmat_line_list[9])

        self.deletion = bmat_line_list[10]
        self.ambiguous = bmat_line_list[11]
        self.insertion = bmat_line_list[12]
        self.mut_num = int(bmat_line_list[13])

        # make dict
        self.count_dict = {
            "A": int(bmat_line_list[3]),
            "G": int(bmat_line_list[4]),
            "C": int(bmat_line_list[5]),
            "T": int(bmat_line_list[6])
        }


#################################################################################
# FUN
#################################################################################
def cmp_site_index(site_index_a, site_index_b, ref_order_dict):
    """
    INPUT:
        <site_index_a>, <site_index_b>
            str, like chr1_629627_CT
        
        <ref_order>
            dict, format like:
            
            ref_order_dict = {
                'chr1': 0, 
                'chr19': 1, 
                'chr20': 2
            }
        
    RETURN:
        0, site_index_a == site_index_b at position level, ignore mutation info;
        -1, site_a at upstream of site_b;
        1, site_a at downstream of site_b
    """

    site_A = siteIndex(site_index_a)
    site_B = siteIndex(site_index_b)

    if site_A.chr_name == site_B.chr_name:
        if site_A.site_pos == site_B.site_pos:
            return 0

        elif site_A.site_pos < site_B.site_pos:
            return -1

        elif site_A.site_pos > site_B.site_pos:
            return 1

    else:
        site_A_chr_order = ref_order_dict.get(site_A.chr_name)
        site_B_chr_order = ref_order_dict.get(site_B.chr_name)

        if (site_A_chr_order is not None) and (site_B_chr_order is not None):
            if site_A_chr_order < site_B_chr_order:
                return -1

            elif site_A_chr_order > site_B_chr_order:
                return 1

        else:
            raise TypeError("Site index not in your reference!")


#################################################################################
# FUN
################################################################################# 
def query_region_bmat_info(bmat_file, site_index_list, genome_order_dict):
    """
    INPUT:
        <bmat_file>
            file.obj, bmat file handle

        <site_index_list>
            list, like [chr1_20452_CT, chr1_20467_C., chr1_20474_CT]
            
        <genome_order_dict>
            dict, for FUN <cmp_site_index>
    
    RETURN
        <site_dict>
            dict, key is site_index, value bmat line list
        
    """

    # init 
    site_index = siteIndex(site_index_list[0])
    bmat_line = bmat_file.readline()

    # define dict
    query_total_num = len(site_index_list)
    query_site_num = 0
    query_res_dict = {
        "site_index_list": []
    }

    # make init val
    for raw_site_index in site_index_list:
        raw_site_index_obj = siteIndex(raw_site_index)
        query_res_dict["site_index_list"].append(raw_site_index_obj.site_index)
        query_res_dict[raw_site_index_obj.site_index] = None

    # query run
    while bmat_line != "":
        bmat_line_list = bmat_line.strip().split("\t")
        bmat_site_index = "_".join(bmat_line_list[0:3])

        cmp_res = cmp_site_index(site_index.site_index_raw, bmat_site_index, genome_order_dict)

        if cmp_res == 0:
            query_site_num += 1

            # add info into dict 
            # query_res_dict["site_index_list"].append(site_index.site_index)
            query_res_dict[site_index.site_index] = bmatLine(bmat_line_list)

            if query_site_num >= query_total_num:
                break

            else:
                # read new
                site_index = siteIndex(site_index_list[query_site_num])
                bmat_line = bmat_file.readline()

        elif cmp_res == 1:
            bmat_line = bmat_file.readline()

        elif cmp_res == -1:
            query_site_num += 1

            # add info into dict 
            # query_res_dict["site_index_list"].append(site_index.site_index)
            query_res_dict[site_index.site_index] = None

            if query_site_num >= query_total_num:
                break

            else:
                # read new
                site_index = siteIndex(site_index_list[query_site_num])
                bmat_line = bmat_file.readline()

    return query_res_dict


#################################################################################
# FUN
#################################################################################
def site_binomial_test(bmat_line_obj, mut_type, background_pval, return_int=True):
    """
    INPUT:
        <bmat_line_obj>
            obj, bmatLine

        <mut_type>
            str, "CT" mean ref_base is "C" and mut_base is "T"

        <background_pval>
            float, binomial test background pvalue

        <return_int>
            bool, set True will return a close int after -10 * log(pvalue), OR return raw pvalue

    RETURN
        <pval>
            int|float, binomial test pvalue
    """

    # count base info
    ref_base_count = bmat_line_obj.count_dict.get(mut_type[0])
    mut_base_count = bmat_line_obj.count_dict.get(mut_type[1])
    total_base_count = ref_base_count + mut_base_count

    try:
        binom_pval = stats.binom_test(x=mut_base_count,
                                      n=total_base_count,
                                      p=background_pval,
                                      alternative='greater')
    except:
        return None

    if return_int:
        if binom_pval > 1e-127:
            binom_pval_int = int(round(-10 * log(binom_pval, 10), 0))
        else:
            binom_pval_int = 127

        return binom_pval_int

    else:
        return binom_pval


#################################################################################
# FUN
#################################################################################
def site_hard_filter(site_bmat_obj,
                     mut_type,
                     mut_count_cutoff=5,
                     other_mut_count_cutoff=5,
                     other_mut_ratio_cutoff=0.25,
                     cover_min_ratio_check_cutoff=6,
                     cover_up_limit_cutoff=500):
    """
    INPUT:
        <site_bmat_obj>
            obj, bmatLine

        <mut_type>
            str, "CT" mean ref_base is "C" and mut_base is "T"

        <mut_count_cutoff>
            int, site larger than this will be blocked

        <other_mut_count_cutoff>
            int, site other type mutation larger than this will be blocked

        <other_mut_ratio_cutoff>
            float, site other type mutation ratio larger than this will be blocked

        <cover_min_ratio_check_cutoff>
            int, mutation count have to larger than this can process with <other_mut_ratio_cutoff>

        <cover_up_limit_cutoff>
            int, site cover larger than this cutoff will be blocked.

    RETURN:
        <block_state, block_reason>
            tuple, contain two element

            <block_state>
                True, site should be blocked

                False, site should not be blocked

            <block_reason>
                MC, too much query mutation count;
                MR, too high query mutation ratio;
                OC, too much other mutation count;
                OR, too high other mutation ratio;
                TC, too much cover count
                None, when <block_state> is False, <block_reason> is None
    """

    # get number
    ref_base_count = site_bmat_obj.count_dict[mut_type[0]]
    mut_base_count = site_bmat_obj.count_dict[mut_type[1]]
    other_mut_base_count = 0
    total_cover_count = 0

    # check ref base count
    if ref_base_count < cover_min_ratio_check_cutoff:
        return False, None

    # check ref upper cutoff
    if ref_base_count >= cover_up_limit_cutoff:
        return True, "TC"

    # check mut count
    if mut_base_count >= mut_count_cutoff:
        return True, "MC"

    # count ratio
    for base in site_bmat_obj.count_dict:
        total_cover_count += site_bmat_obj.count_dict[base]
        if base not in mut_type:
            other_mut_base_count += site_bmat_obj.count_dict[base]

    other_mut_ratio = other_mut_base_count / 1.0 / total_cover_count

    # check other mut type
    if other_mut_base_count >= other_mut_count_cutoff:
        return True, "OC"

    if other_mut_ratio >= other_mut_ratio_cutoff:
        return True, "OR"

    # final return
    return False, None


#################################################################################
# FUN
#################################################################################
def add_block_info_to_mpmat(
        chr_mpmat_filename,
        chr_bmat_ctrl_filename,
        chr_bmat_treat_filename,
        query_mutation_type,
        chr_ctrl_mut_bg_pval,
        ref_order_dict,
        out_chr_mpmat_filename="stdout",
        hf_ctrl_mut_count_cutoff=3,
        hf_ctrl_other_mut_count_cutoff=5,
        hf_ctrl_other_mut_ratio_cutoff=0.25,
        hf_ctrl_cover_min_ratio_check_cutoff=6,
        hf_ctrl_cover_up_limit_cutoff=500,
        ctrl_binomial_cutoff_int=30,
        hf_treat_mut_count_cutoff=5000,
        hf_treat_other_mut_count_cutoff=50,
        hf_treat_other_mut_ratio_cutoff=0.6,
        hf_treat_cover_min_ratio_check_cutoff=10,
        hf_treat_cover_up_limit_cutoff=5000,
        treat_site_mut_min_cutoff=1,
        log_verbose=3
):
    """
    HELP:

    """

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.basicConfig(level=10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # part I open file
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    try:
        chr_mpmat_file = open(chr_mpmat_filename, "r")

        # open ctrl bmat file
        if (chr_bmat_ctrl_filename[-3:] == ".gz") or (chr_bmat_ctrl_filename[-5:] == ".gzip"):
            chr_bmat_ctrl_file = gzip.open(chr_bmat_ctrl_filename, "r")
        else:
            chr_bmat_ctrl_file = open(chr_bmat_ctrl_filename, "r")

        # open treat bmat file
        if (chr_bmat_treat_filename[-3:] == ".gz") or (chr_bmat_treat_filename[-5:] == ".gzip"):
            chr_bmat_treat_file = gzip.open(chr_bmat_treat_filename, "r")
        else:
            chr_bmat_treat_file = open(chr_bmat_treat_filename, "r")

        if out_chr_mpmat_filename == "stdout":
            out_chr_mpmat_file = sys.stdout
        else:
            out_chr_mpmat_file = open(out_chr_mpmat_filename, "w")
    except:
        raise IOError("Load or Open files error!")

    logging.info("Start to analysis sites block info. \n\t%s" % chr_mpmat_filename)
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # part II iteration and add block info
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    for line_index, line in enumerate(chr_mpmat_file):

        line_list = line.strip().split("\t")
        mpmat_info = mpmatLine(line_list)

        # log
        if log_verbose == 0:
            run_report_num = 1000
        elif log_verbose == 1:
            run_report_num = 100
        elif log_verbose == 2:
            run_report_num = 50
        else:
            run_report_num = 20

        if line_index % run_report_num == 0:
            logging.info("Running mpmat block step on %s \n\tProcessed line number %s" % (chr_mpmat_filename, line_index + 1))

        # check mut type
        if mpmat_info.mut_type not in query_mutation_type:
            raise ValueError("mpmatLine.mut_type doesn't match query mutation type!")

        # ---------------------------------------------------------->>>>>>>>>>
        # query site bmat info
        # ---------------------------------------------------------->>>>>>>>>>
        query_site_idx_list = mpmat_info.site_index_list

        ctrl_bmat_res_dict = query_region_bmat_info(bmat_file=chr_bmat_ctrl_file,
                                                    site_index_list=query_site_idx_list,
                                                    genome_order_dict=ref_order_dict)

        treat_bmat_res_dict = query_region_bmat_info(bmat_file=chr_bmat_treat_file,
                                                     site_index_list=query_site_idx_list,
                                                     genome_order_dict=ref_order_dict)

        # ---------------------------------------------------------->>>>>>>>>>
        # ctrl Hard filter
        # ---------------------------------------------------------->>>>>>>>>>
        ctrl_site_hard_filter_state = []
        ctrl_site_hard_filter_reason = []

        for site_idx in ctrl_bmat_res_dict["site_index_list"]:
            if ctrl_bmat_res_dict[site_idx] is None:
                ctrl_site_hard_filter_state.append(False)

                # NS means 'Non sequencing data'
                ctrl_site_hard_filter_reason.append("NS")

            else:
                hard_filter_res = site_hard_filter(site_bmat_obj=ctrl_bmat_res_dict[site_idx],
                                                   mut_type=mpmat_info.mut_type,
                                                   mut_count_cutoff=hf_ctrl_mut_count_cutoff,
                                                   other_mut_count_cutoff=hf_ctrl_other_mut_count_cutoff,
                                                   other_mut_ratio_cutoff=hf_ctrl_other_mut_ratio_cutoff,
                                                   cover_min_ratio_check_cutoff=hf_ctrl_cover_min_ratio_check_cutoff,
                                                   cover_up_limit_cutoff=hf_ctrl_cover_up_limit_cutoff)

                ctrl_site_hard_filter_state.append(hard_filter_res[0])
                ctrl_site_hard_filter_reason.append(hard_filter_res[1])

        # ---------------------------------------------------------->>>>>>>>>>
        # ctrl Binomial test
        # ---------------------------------------------------------->>>>>>>>>>
        binom_test_pval = []
        binom_test_res = 0

        for index, site_idx in enumerate(ctrl_bmat_res_dict["site_index_list"]):
            if not ctrl_site_hard_filter_state[index]:
                if ctrl_site_hard_filter_reason[index] != "NS":
                    binom_test_res = site_binomial_test(bmat_line_obj=ctrl_bmat_res_dict[site_idx],
                                                        mut_type=mpmat_info.mut_type,
                                                        background_pval=chr_ctrl_mut_bg_pval,
                                                        return_int=True)

            binom_test_pval.append(binom_test_res)

        # ---------------------------------------------------------->>>>>>>>>>
        # treat Hard filter
        # ---------------------------------------------------------->>>>>>>>>>
        treat_site_hard_filter_state = []
        treat_site_hard_filter_reason = []

        for site_idx in treat_bmat_res_dict["site_index_list"]:
            if treat_bmat_res_dict[site_idx] is None:
                treat_site_hard_filter_state.append(True)

                # NS means 'Non sequencing data'
                treat_site_hard_filter_reason.append("NS")

            else:
                hard_filter_res = site_hard_filter(site_bmat_obj=treat_bmat_res_dict[site_idx],
                                                   mut_type=mpmat_info.mut_type,
                                                   mut_count_cutoff=hf_treat_mut_count_cutoff,
                                                   other_mut_count_cutoff=hf_treat_other_mut_count_cutoff,
                                                   other_mut_ratio_cutoff=hf_treat_other_mut_ratio_cutoff,
                                                   cover_min_ratio_check_cutoff=hf_treat_cover_min_ratio_check_cutoff,
                                                   cover_up_limit_cutoff=hf_treat_cover_up_limit_cutoff)

                treat_site_hard_filter_state.append(hard_filter_res[0])
                treat_site_hard_filter_reason.append(hard_filter_res[1])

        # ---------------------------------------------------------->>>>>>>>>>
        # make block state
        # ---------------------------------------------------------->>>>>>>>>>
        block_state_list = [False] * mpmat_info.site_num
        block_reason_list = ["None"] * mpmat_info.site_num

        for index in range(mpmat_info.site_num):
            try:
                if ctrl_site_hard_filter_state[index]:
                    block_state_list[index] = True
                    # block reason: ctrl hard filter
                    block_reason_list[index] = "CHF"
                    continue
            except:
                print(ctrl_site_hard_filter_state)
                print(ctrl_bmat_res_dict)

            if binom_test_pval[index] > ctrl_binomial_cutoff_int:
                block_state_list[index] = True
                # block reason: ctrl binomial test
                block_reason_list[index] = "CBT"
                continue

            # block by treat info
            site_idx = treat_bmat_res_dict["site_index_list"][index]

            if treat_bmat_res_dict[site_idx] is None:
                block_state_list[index] = True
                # block reason: treat non cover
                block_reason_list[index] = "TNC"
                continue

            else:
                site_mut_count = treat_bmat_res_dict[site_idx].count_dict[mpmat_info.mut_type[-1]]
                if site_mut_count < treat_site_mut_min_cutoff:
                    block_state_list[index] = True
                    # block reason: treat non mut
                    block_reason_list[index] = "TNM"
                    continue

            if treat_site_hard_filter_state[index]:
                block_state_list[index] = True
                # block reason: treat hard filter
                block_reason_list[index] = "THF"
                continue

        # ---------------------------------------------------------->>>>>>>>>>
        # out block info
        # ---------------------------------------------------------->>>>>>>>>>
        out_block_info_list = [
            "CHF,CHR,CBT,THF,THR,FBS,FBR"
        ]

        info_list = [
            ctrl_site_hard_filter_state,
            ctrl_site_hard_filter_reason,
            binom_test_pval,
            treat_site_hard_filter_state,
            treat_site_hard_filter_reason,
            block_state_list,
            block_reason_list
        ]

        for info in info_list:
            out_block_info_list.append(",".join(map(str, info)))

        # out to file
        line_list.append(" ".join(out_block_info_list))
        out_chr_mpmat_file.write("\t".join(line_list) + "\n")

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # close file
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    chr_mpmat_file.close()
    chr_bmat_ctrl_file.close()
    chr_bmat_treat_file.close()

    if out_chr_mpmat_filename != "stdout":
        out_chr_mpmat_file.close()

    logging.info("Done! \n\t%s" % chr_mpmat_filename)

    return 1


#################################################################################
# FUN
#################################################################################
def multi_mpmat_block_site(
        mpmat_split_dict,
        ctrl_bmat_split_dict,
        treat_bmat_split_dict,
        query_mutation_type,
        ctrl_mut_bg_pval_dict,
        ref_order_dict,
        thread=1,
        log_verbose=3,
        **filter_args):
    """
    INPUT:

        <**filter_args>
            args list:
                <hf_ctrl_mut_count_cutoff>
                <hf_ctrl_other_mut_count_cutoff>
                <hf_ctrl_other_mut_ratio_cutoff>
                <hf_ctrl_cover_min_ratio_check_cutoff>
                <hf_ctrl_cover_up_limit_cutoff>
                <ctrl_binomial_cutoff_int>
                <hf_treat_mut_count_cutoff>
                <hf_treat_other_mut_count_cutoff>
                <hf_treat_other_mut_ratio_cutoff>
                <hf_treat_cover_min_ratio_check_cutoff>
                <hf_treat_cover_up_limit_cutoff>
                <treat_site_mut_min_cutoff>

    RETURN
        <out_mpmat_filename_dict>
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
        "hf_ctrl_mut_count_cutoff": 3,
        "hf_ctrl_other_mut_count_cutoff": 5,
        "hf_ctrl_other_mut_ratio_cutoff": 0.25,
        "hf_ctrl_cover_min_ratio_check_cutoff": 6,
        "hf_ctrl_cover_up_limit_cutoff": 500,
        "ctrl_binomial_cutoff_int": 30,
        "hf_treat_mut_count_cutoff": 5000,
        "hf_treat_other_mut_count_cutoff": 50,
        "hf_treat_other_mut_ratio_cutoff": 0.6,
        "hf_treat_cover_min_ratio_check_cutoff": 10,
        "hf_treat_cover_up_limit_cutoff": 5000,
        "treat_site_mut_min_cutoff": 1
    }

    for args_key in filter_args:
        if filter_args.get(args_key) is not None:
            run_params_dict[args_key] = filter_args[args_key]

    # ------------------------------------------------------------>>>>>>>>>>
    # var initial
    # ------------------------------------------------------------>>>>>>>>>>
    # check mpmat chr_name and bmat chr_name
    chr_name_order_list = []
    for mpmat_chr_name in mpmat_split_dict["chr_name_order"]:
        if ref_order_dict.get(mpmat_chr_name) is not None:
            chr_name_order_list.append(mpmat_chr_name)

    out_mpmat_filename_dict = {"chr_name_order": chr_name_order_list}
    # ------------------------------------------------------------>>>>>>>>>>
    # run part
    # ------------------------------------------------------------>>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    logging.info("Start to calculate genome block info...")

    pool = multiprocessing.Pool(processes=thread)

    run_return_info_list = []

    for chr_name in chr_name_order_list:
        # mpmat input and output
        chr_mpmat_old_filename = mpmat_split_dict[chr_name]
        chr_mpmat_new_filename = chr_mpmat_old_filename + "." + "AddBlockInfo"
        out_mpmat_filename_dict[chr_name] = chr_mpmat_new_filename

        # bmat
        ctrl_bmat_filename = ctrl_bmat_split_dict[chr_name]
        treat_bmat_filename = treat_bmat_split_dict[chr_name]

        run_return_info_list.append(
            pool.apply_async(
                func=add_block_info_to_mpmat,
                args=(
                    chr_mpmat_old_filename,
                    ctrl_bmat_filename,
                    treat_bmat_filename,
                    query_mutation_type,
                    ctrl_mut_bg_pval_dict[chr_name]["query_mut_bg_pval"],
                    ref_order_dict,
                    chr_mpmat_new_filename,
                    run_params_dict["hf_ctrl_mut_count_cutoff"],
                    run_params_dict["hf_ctrl_other_mut_count_cutoff"],
                    run_params_dict["hf_ctrl_other_mut_ratio_cutoff"],
                    run_params_dict["hf_ctrl_cover_min_ratio_check_cutoff"],
                    run_params_dict["hf_ctrl_cover_up_limit_cutoff"],
                    run_params_dict["ctrl_binomial_cutoff_int"],
                    run_params_dict["hf_treat_mut_count_cutoff"],
                    run_params_dict["hf_treat_other_mut_count_cutoff"],
                    run_params_dict["hf_treat_other_mut_ratio_cutoff"],
                    run_params_dict["hf_treat_cover_min_ratio_check_cutoff"],
                    run_params_dict["hf_treat_cover_up_limit_cutoff"],
                    run_params_dict["treat_site_mut_min_cutoff"],
                )
            )
        )

    pool.close()
    pool.join()

    # check run state
    final_run_state = 0
    for index, res in enumerate(run_return_info_list):
        run_state = res.get()
        if run_state != 1:
            logging.error("Blocking error occur with %s!" % chr_name_order_list[index])
            if final_run_state == 0:
                final_run_state = 1

    if final_run_state == 0:
        logging.info("Calculation of genome block info. Done!")
        return out_mpmat_filename_dict

    elif final_run_state == 1:
        logging.error("Something wrong with block step!")
        raise RuntimeError("")


#################################################################################
# FUN
#################################################################################
def query_site_mpileup_info(site_idx_list,
                            bam_obj,
                            genome_obj,
                            ignore_overlaps=True,
                            min_base_quality=20,
                            min_mapping_quality=20,
                            dist_cutoff=10000):
    """
    INPUT
        <site_idx_list>
            list, format like:

                ["chr1_125178576_CT", "chr1_125178578_CT", "chr1_125178580_CT", "chr1_125178588_CA"]

            The site coordinate is related to the UCSC genome browser index.

            Site index have to be sorted.

        <bam_obj>
            pysam.AlignmentFile() obj

        <genome_obj>
            genome obj, which is created by pysam.FastaFile()

    RETURN
        dict
            key:
                site_index only key chr_name and chr_pos like "chr1_125178576" rather than "chr1_125178576_CT"

            value:
                A, T, C, G, N, other count
    """

    # -------------------------------------------------------->>>>>>>>
    # check site index list
    # -------------------------------------------------------->>>>>>>>
    if len(site_idx_list) == 0:
        return None

    chr_name = site_idx_list[0].split("_")[0]
    region_start = region_end = int(site_idx_list[0].split("_")[1])
    site_dist = 0

    for site_idx in site_idx_list:
        if site_idx.split("_")[0] != chr_name:
            raise IOError("<site_idx_list> all site index have to on the same chromosome.")

    if len(site_idx_list) >= 2:
        region_start = int(site_idx_list[0].split("_")[1])
        region_end = int(site_idx_list[-1].split("_")[1])
        site_dist = region_end - region_start

    if site_dist >= dist_cutoff:
        raise IOError("<site_idx_list> mpileup region is too large, which could take a lot of time!")

    # -------------------------------------------------------->>>>>>>>
    # make raw dict
    # -------------------------------------------------------->>>>>>>>
    site_dict = {}
    for site_idx in site_idx_list:
        key = "_".join(site_idx.split("_")[:2])
        site_dict[key] = {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0, "total": 0}

    # -------------------------------------------------------->>>>>>>>
    # run mpileup
    # -------------------------------------------------------->>>>>>>>
    mpileup_extend_length = 10

    mpileup_iter = bam_obj.pileup(chr_name,
                                  region_start - mpileup_extend_length,
                                  region_end + mpileup_extend_length,
                                  fastafile=genome_obj,
                                  ignore_overlaps=ignore_overlaps,
                                  min_base_quality=min_base_quality,
                                  min_mapping_quality=min_mapping_quality,
                                  stepper="samtools")

    for pileup in mpileup_iter:
        run_index = "%s_%s" % (chr_name, pileup.reference_pos + 1)

        if (pileup.reference_pos + 1) > region_end:
            break

        if site_dict.get(run_index) is not None:
            base_list = list(map(str.upper, pileup.get_query_sequences()))
            site_dict[run_index]["A"] = base_list.count("A")
            site_dict[run_index]["T"] = base_list.count("T")
            site_dict[run_index]["C"] = base_list.count("C")
            site_dict[run_index]["G"] = base_list.count("G")
            site_dict[run_index]["N"] = base_list.count("N")
            site_dict[run_index]["total"] = len(base_list)

    return site_dict


#################################################################################
# FUN
#################################################################################
def find_block_info_and_highest_signal(
        mpmat_info,
        ctrl_bam_obj,
        treat_bam_obj,
        ref_genome_obj,
        block_mut_num_cutoff=2
):
    """
    INPUT:
        <ctrl_bam_obj> and <treat_bam_obj>
            pysam.AlignmentFile obj

        <ref_genome_obj>
            pysam.Fasta obj

        <mpmat_info>
            mpmatLine obj

        <block_mut_num_cutoff>
            int, site in ctrl sample contain mutation number >= this cutoff will be blocked in the following steps

    RETURN
        <back_mpmat_line>
            mpmatLine obj with additional features
                1. mpmatLine.block_info_list
                    list, contain True or False, True means need to block in the following steps.

                2. mpmatLine.highest_site_dict
                    dict, format like
                        {'A': 0,
                         'T': 6,
                         'C': 96,
                         'G': 0,
                         'N': 0,
                         'total': 102,
                         'site_index': 'chr1_121885716',
                         'full_site_index': 'chr1_121885716_CT',
                         'mut_num': 6,
                         'mut_ratio': 0.058823529411764705}

                3. block_site_num
                    int, count block site number

    """

    # init params
    if not hasattr(mpmat_info, "block_info_list"):
        mpmat_info.block_info_list = [False] * mpmat_info.site_num
        mpmat_info.block_site_num = 0

    mut_type_base = mpmat_info.mut_type[1]

    highest_treat_mut_site_dict = {
        "A": 0,
        "T": 0,
        "G": 0,
        "C": 0,
        "N": 0,
        "total": 0,
        "site_index": "",
        "full_site_index": "",
        "mut_num": 0,
        "mut_ratio": 0,
        "run_state": False
    }

    # query ctrl and treat pileup info
    collect_ctrl_site_cover_dict = query_site_mpileup_info(site_idx_list=mpmat_info.site_index_list,
                                                           bam_obj=ctrl_bam_obj,
                                                           genome_obj=ref_genome_obj,
                                                           ignore_overlaps=True,
                                                           min_base_quality=20,
                                                           min_mapping_quality=20,
                                                           dist_cutoff=10000)

    collect_treat_site_cover_dict = query_site_mpileup_info(site_idx_list=mpmat_info.site_index_list,
                                                            bam_obj=treat_bam_obj,
                                                            genome_obj=ref_genome_obj,
                                                            ignore_overlaps=True,
                                                            min_base_quality=20,
                                                            min_mapping_quality=20,
                                                            dist_cutoff=10000)

    # add block info and find highest signal
    for run_idx, full_site_index in enumerate(mpmat_info.site_index_list):
        if mpmat_info.SNP_ann_list[run_idx]:
            mpmat_info.block_info_list[run_idx] = True
            continue

        # fix site index
        site_index = "_".join(full_site_index.split("_")[:2])

        # query site cover in ctrl sample
        ctrl_site_cover_dict = collect_ctrl_site_cover_dict.get(site_index)

        if ctrl_site_cover_dict is not None:
            ctrl_site_mut_num = ctrl_site_cover_dict.get(mut_type_base)
        else:
            continue

        if ctrl_site_mut_num is None:
            continue

        if ctrl_site_mut_num >= block_mut_num_cutoff:
            mpmat_info.block_info_list[run_idx] = True

        if not mpmat_info.block_info_list[run_idx]:
            treat_site_cover_dict = collect_treat_site_cover_dict.get(site_index)
            if not highest_treat_mut_site_dict["run_state"]:
                highest_treat_mut_site_dict = treat_site_cover_dict.copy()
                highest_treat_mut_site_dict["site_index"] = site_index
                highest_treat_mut_site_dict["full_site_index"] = full_site_index
                highest_treat_mut_site_dict["run_state"] = True

            if treat_site_cover_dict[mut_type_base] > highest_treat_mut_site_dict[mut_type_base]:
                highest_treat_mut_site_dict = treat_site_cover_dict.copy()
                highest_treat_mut_site_dict["site_index"] = site_index
                highest_treat_mut_site_dict["full_site_index"] = full_site_index
                highest_treat_mut_site_dict["run_state"] = True

    # calculate highest info and signal
    if highest_treat_mut_site_dict["total"] > 0:
        highest_treat_mut_site_dict["mut_num"] = highest_treat_mut_site_dict[mut_type_base]
        highest_treat_mut_site_dict["mut_ratio"] = highest_treat_mut_site_dict["mut_num"] / 1.0 / highest_treat_mut_site_dict["total"]
    else:
        highest_treat_mut_site_dict["mut_num"] = 0
        highest_treat_mut_site_dict["mut_ratio"] = 0.0
        highest_treat_mut_site_dict["run_state"] = False

    # add dict into mpmat info
    mpmat_info.highest_site_dict = highest_treat_mut_site_dict.copy()

    # fix mut_key_list
    for index, site_index in enumerate(mpmat_info.site_index_list):
        if mpmat_info.block_info_list[index]:
            mpmat_info.mut_key_list[index] = "B"

        if mpmat_info.SNP_ann_list[index]:
            mpmat_info.mut_key_list[index] = "S"

    # fix mut_key
    mpmat_info.mut_key = "-".join(mpmat_info.mut_key_list)

    # add block_site_num
    mpmat_info.block_site_num = mpmat_info.block_info_list.count(True)

    return mpmat_info
