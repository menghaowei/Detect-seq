#! /Users/meng/menghw_HD/miniconda3/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import sys
import os
import subprocess
import logging

import string
import random

import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Align.substitution_matrices import Array

from typing import Dict, Union, Any

# Version information START --------------------------------------------------
VERSION_INFO = """
Author: MENG Howard

Version-06-TALE:
    2022-02-12
        TALE sequence alignment based on V6 code

Version-06:
    2022-01-05
        Test gap open penalty, especially in PAM distal side

    2021-12-07
        Add signal dist penalty;
        Support .bed format table
        Support Python3.x
        Reconstruction all code!
        
Version-05:
    2021-11-26
        Add signal dist penalty;
        Support .bed format table 

Version-04:
    2020-09-02
        Fix output region_id 'nan' error.

Version-03:
    2020-05-22
        Add select alignment mode function

Version-02:
    2020-03-06
        Fix d=0 issue and add alignment strategy option

Version-01:
    2019-12-03 
        mpmat file convert to alignment result table (.art) file 

E-Mail: meng_howard@126.com

"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = """

There will be a lot of parameters to set,

but default parameters can handle output figure in most cases.

"""

PIPELINE = """

INPUT
    <mpmat> or <BED> format 
    
    <fasta>
    
    <alignment settings>
    
    <id_prefix>
    
OUTPUT
    <alignment_result_table>
        format like:

            
DESIGN
    0. use bedtools getfasta to get fasta temp files
    
    1. read all mpmat and fa file 
    
    2. seed alignment (+/-)
    
    3. extend alignment (+/-)
    
    4. select the best alignment
    
    5. calculate output info (.art) file
"""


# Learning Part END-----------------------------------------------------------


###############################################################################
# function part 
###############################################################################
def _log_cmd_str(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <full_cmd_str>
            str, record full command str info
    """
    full_cmd_str = """python mpmat-to-art.py
    --input {input}
    --query_seq {query_seq}
    --output {output}                    
    --reference {reference}
    --extend_method {extend_method}
    --extend_length {extend_length}
    --input_filetype {input_filetype}
    --input_header {input_header}
    --mpmat_fwd_mut_type {mpmat_fwd_mut_type}
    --mpmat_rev_mut_type {mpmat_rev_mut_type}    
    --bedtools_path {bedtools_path}
    --align_settings {align_settings}
    --sub_mat_alphabet {sub_mat_alphabet}
    --sub_mat_specific {sub_mat_specific}
    --N0_T_award {N0_T_award}
    --temp_dir {temp_dir}
    --verbose {verbose}""".format(
        input=args.input,
        query_seq=args.query_seq,
        output=args.output,
        reference=args.reference,
        extend_method=args.extend_method,
        extend_length=args.extend_length,
        input_filetype=args.input_filetype,
        input_header=args.input_header,
        mpmat_fwd_mut_type=args.mpmat_fwd_mut_type,
        mpmat_rev_mut_type=args.mpmat_rev_mut_type,
        bedtools_path=args.bedtools_path,
        align_settings=args.align_settings,
        sub_mat_alphabet=args.sub_mat_alphabet,
        sub_mat_specific=args.sub_mat_specific,
        N0_T_award=args.N0_T_award,
        temp_dir=args.temp_dir,
        verbose=args.verbose
    )

    return full_cmd_str


def find_highest_mut_site(mpmat_line_list, find_method="score"):
    """
    INPUT:
        <mpmat_row>
            mpmat_row, pandas Series obj, with default column names

        <find_method>
            1. ratio
                return 5' site with highest mutation ratio;

            2. mut_count
                return 5' site with highest mutation reads count;

            3. score
                    score = mut_count * ratio ^ 2
                return 5' site with highest score;
    """
    site_num = int(mpmat_line_list[3])

    # fetch site index
    site_index_list = mpmat_line_list[6].split(",")

    # fetch SNP state
    snp_ann_list = mpmat_line_list[10].split(",")

    # fetch count and ratio
    mut_base_num_list = list(map(int, mpmat_line_list[7].split(",")))
    cover_base_num_list = list(map(int, mpmat_line_list[8].split(",")))
    mut_ratio_list = list(map(float, mpmat_line_list[9].split(",")))

    # select site index
    highest_site_index = site_index_list[0]
    highest_index = 0

    if find_method == "ratio":
        highest_index = mut_ratio_list.index(max(mut_ratio_list))

    elif find_method == "mut_count":
        highest_index = mut_base_num_list.index(max(mut_base_num_list))

    elif find_method == "score":
        mut_score_list = []
        for i in range(site_num):
            if snp_ann_list[i] != "True":
                mut_score = mut_ratio_list[i] * mut_ratio_list[i] * mut_base_num_list[i]
                mut_score_list.append(mut_score)
            else:
                mut_score_list.append(0)
        highest_index = mut_score_list.index(max(mut_score_list))

    # get site_index
    highest_site_index = site_index_list[highest_index]

    # return part
    chr_name, chr_pos = highest_site_index.split("_")[:2]
    mut_count = mut_base_num_list[highest_index]
    mut_ratio = mut_ratio_list[highest_index]
    cover_count = cover_base_num_list[highest_index]

    return chr_name, int(chr_pos), highest_site_index, mut_count, cover_count, mut_ratio


def get_mut_type_count(site_index_list_str):
    """
    INPUT
        <site_index_list_str>
            A str like 'chr1_143265617_CT,chr1_143265637_CT,chr1_143265643_CT,chr1_143265655_CT'

    RETURN
        count mut_type and return the max count mut_type
    """

    mut_type_dict = {}

    site_index_list = site_index_list_str.split(",")

    for site_index in site_index_list:
        mut_type = site_index.split("_")[-1]
        if mut_type_dict.get(mut_type) is None:
            mut_type_dict[mut_type] = 1
        else:
            mut_type_dict[mut_type] += 1

    return sorted(mut_type_dict.items(), key=lambda i: i[1], reverse=True)[0][0]


def make_temp_fasta(input_filename,
                    input_filetype,
                    genome_fasta,
                    temp_dir,
                    mut_type_rev="GA",
                    bedtools_path="bedtools",
                    header_state=False,
                    extend_region_size=20,
                    extend_region_method="highest_site",
                    log_verbose=3):
    """
    INPUT:
        <input_filename>
            str, .bed file or .mpmat file

        <input_filetype>
            str, 'bed' or 'mpmat'

        <genome_fasta>
            str, reference genome fasta file

        <temp_dir>
            str, temp dir to store temp files

        <header_state>
            bool, if input file contains a header line

        <extend_region_size>
            int,

        <select_region_method>
            1. region
            2. upstream_site
            3. downstream_site
            4. highest_site

    RETURN:
        If successful return
            (0, temp_fasta_file, temp_bed_file)

        If unsuccessful return
            (1, None, None)

    """
    # logging
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # make temp BED file
    temp_bed_basename = "temp" + "_" + "".join(
        random.sample(string.ascii_letters + string.digits, 16)) + ".bed"
    in_fun_temp_bed_filename = os.path.abspath(os.path.join(temp_dir, temp_bed_basename))
    temp_bed_file = open(in_fun_temp_bed_filename, "w")

    with open(input_filename, "r") as in_file:
        if header_state:
            line = in_file.readline()

        for line in in_file:
            line_list = line.split("\t")

            region_chr = line_list[0]
            region_start = int(line_list[1])
            region_end = int(line_list[2])
            region_index = "%s_%s_%s" % (region_chr, region_start, region_end)

            # get region strand
            strand_info = "+"
            if input_filetype == "bed":
                if len(line_list) >= 6:
                    if line_list[5] == "-":
                        strand_info = "-"

            elif input_filetype == "mpmat":
                region_mut_type = get_mut_type_count(line_list[6])
                if region_mut_type == mut_type_rev:
                    strand_info = "-"

            # extend region
            if extend_region_method == "region":
                region_extend_start = region_start - extend_region_size
                region_extend_end = region_end + extend_region_size

            elif extend_region_method == "upstream_site":
                if strand_info == "+":
                    region_extend_start = region_start - extend_region_size
                    region_extend_end = region_start + extend_region_size
                else:
                    region_extend_start = region_end - extend_region_size
                    region_extend_end = region_end + extend_region_size

            elif extend_region_method == "downstream_site":
                if strand_info == "+":
                    region_extend_start = region_end - extend_region_size
                    region_extend_end = region_end + extend_region_size
                else:
                    region_extend_start = region_start - extend_region_size
                    region_extend_end = region_start + extend_region_size

            elif extend_region_method == "highest_site":
                if input_filetype == "mpmat":
                    site_chr_name, site_chr_pos = find_highest_mut_site(mpmat_line_list=line_list, find_method="score")
                    region_extend_start = site_chr_pos - extend_region_size
                    region_extend_end = site_chr_pos + extend_region_size

                else:
                    region_extend_start = region_start - extend_region_size
                    region_extend_end = region_end + extend_region_size

            else:
                raise KeyError("Wrong extend method!")

            # output info into bed
            bed_out_list = [region_chr, region_extend_start, region_extend_end, region_index, 1, strand_info]
            temp_bed_file.write("\t".join(map(str, bed_out_list)) + "\n")

    temp_bed_file.close()

    # log
    logging.info("Output temp bed with filename: %s\n" % in_fun_temp_bed_filename)

    # use bedtools to make fasta file
    temp_fasta_filename = in_fun_temp_bed_filename + ".fa"

    # only report sequence same to the genome fwd
    tool_bedtools_cmd_fmt = "{bedtools} getfasta -fi {ref_fasta}  -fo {out_fasta} -bed {input_bed} -name+ "
    tool_bedtools_cmd = tool_bedtools_cmd_fmt.format(
        bedtools=bedtools_path,
        ref_fasta=genome_fasta,
        out_fasta=temp_fasta_filename,
        input_bed=in_fun_temp_bed_filename
    )

    # bedtools getfasta is 0-based index system, while UCSC is 1-based index system
    # the real index
    # chr_start = bedtools.start + 1
    # chr_end = bedtools.end

    # log
    logging.info("bedtools command:\n%s\n" % tool_bedtools_cmd)

    # run bedtools
    get_fasta_runcode = subprocess.call(tool_bedtools_cmd, shell=True)

    # return
    if get_fasta_runcode != 0:
        os.remove(temp_fasta_filename)
        os.remove(in_fun_temp_bed_filename)
        logging.error("Error! Unsuccessfully get region FASTA seq!")
        return 1, None, None

    else:
        logging.info("Successfully get region FASTA seq!\n")
        return 0, temp_fasta_filename, in_fun_temp_bed_filename


# ---------------------------------------------------------------------->>>>>>
# alignment related functions
# ---------------------------------------------------------------------->>>>>>
def _make_np_substitution_mat(mat_match_score=5, mat_mismatch_score=-4, alphabet="AGCTN", add_sub_info=""):
    """
    INPUT:
        <mat_match_score>
            int, default = 5

        <mat_mismatch_score>
            int, default = -4

        <alphabet>
            str, default = "AGCTN", use to build up score matrix

        <add_info>
            str, format like ref:query=score, A:G=5,G:A=3..., comma as sperator

    RETURN:
        <substitution_matrix>
            dict, alignment matrix for aligner
    """
    # raw matrix
    score_matrix = Array(alphabet=alphabet, dims=2, dtype="float")

    for i1, letter1 in enumerate(alphabet):
        for i2, letter2 in enumerate(alphabet):
            if letter1 == letter2:
                score_matrix[i1, i2] = mat_match_score

            elif letter1 != letter2:
                score_matrix[i1, i2] = mat_mismatch_score

    # fix align score
    if add_sub_info != "":
        try:
            add_sub_info_list = add_sub_info.split(",")
            for sub_info in add_sub_info_list:
                mut_info, align_score = sub_info.split("=")
                from_base, to_base = mut_info.split(":")
                score_matrix[from_base.upper(), to_base.upper()] = float(align_score)
        except:
            Warning("Something wrong with add_info, " + "\n")
            raise IOError("Please set add info as format ref:query=score, A:G=5,G:A=3..." + "\n")

    return score_matrix


def degen_align_obj_analysis(alignment, N0_T_award_score=5, degeneracy_base="AG"):
    """
    INPUT:
        <alignment>
            obj, alignment result

    RETURN:
        <info_dict>
            int, gap count
            int, match count
            int, mismatch count
            list, gap position [0 based index]
            list, mismatch position [0 based index]
            float, alignment score, higher means better alignment
            obj, alignment obj

    INFO:
        The gap count should be the num of gap contain in qeury alignment region.

        e.g.

        AGTGGTAGACATT
        ------||X|--|
        ------AGCC--T


        return = {
            "gap_count" : 2,
            "match_count" : 4,
            "mismatch_count" : 1,
            "gap_pos" : [10, 11],
            "mismatch_pos" : [8],
            "alignment_score" : [N0_T_award_score] + [raw alignment score]
            "alignment" : alignment obj,
            "meta": {
                "N0_T_award_score" : N0_T_award_score,
                "degeneracy_base" : degeneracy_base,
            }
        }
    """

    # init dict
    info_dict = {
        "gap_count": 0,
        "match_count": 0,
        "mismatch_count": 0,
        "degen_match_count": 0,
        "degen_count": 0,
        "degen_mismatch_count": 0,
        "gap_pos": [],
        "mismatch_pos": [],
        "degen_match_pos": [],
        "degen_mismatch_pos": [],
        "alignment_score": alignment.score,
        "alignment": alignment,
        "meta": {
            "N0_T_award_score": N0_T_award_score,
            "degeneracy_base": degeneracy_base
        }
    }

    # fetch alignment info
    alignment_list = str(alignment).split("\n")
    query_length = len(alignment.query)

    # load alignment
    ref_list = alignment_list[0]
    info_list = alignment_list[1]
    query_list = alignment_list[2]

    fix_info_list = list(info_list)

    # count init value
    count_state = False
    count_query_base = 0

    query_index = 0

    # counting
    for run_index, info_base in enumerate(info_list):
        if (not count_state) and (query_list[run_index] != "-"):
            count_state = True

            # record N0 info
            if info_dict.get("N0_base") is None:
                info_dict["N0_base"] = alignment_list[0][run_index]

                # make N0 T award
                if info_dict["N0_base"] == "T":
                    info_dict["alignment_score"] += N0_T_award_score

        if count_state:
            if info_base == "|":
                info_dict["match_count"] += 1
                info_dict["degen_match_count"] += 1
                query_index += 1

            elif info_base == "X" or info_base == ".":
                info_dict["mismatch_count"] += 1
                info_dict["mismatch_pos"].append(query_index)
                query_index += 1

                if (ref_list[run_index] == degeneracy_base[0]) and (query_list[run_index] == degeneracy_base[1]):
                    info_dict["degen_count"] += 1
                    info_dict["degen_match_count"] += 1
                    info_dict["degen_match_pos"].append(query_index)
                    fix_info_list[run_index] = "."

                else:
                    info_dict["degen_mismatch_count"] += 1
                    info_dict["degen_mismatch_pos"].append(query_index)
                    fix_info_list[run_index] = "X"

            elif info_base == "-":
                info_dict["gap_count"] += 1
                info_dict["gap_pos"].append(query_index)

                if query_list[run_index] != "-":
                    query_index += 1

        if query_index >= query_length:
            break

    # add dict
    info_dict["alignment_degen_path"] = "".join(fix_info_list)

    if query_index < query_length:
        info_dict["alignment_score"] -= 20

    return info_dict


def run_step_align(align_ref_seq, align_query_seq, aligner,
                   align_window_size=25,
                   align_step_size=13,
                   N0_T_award_score=5,
                   degeneracy_base="AG"):
    """
    INPUT:
        <align_ref_seq>
            str, reference sequence

        <align_query_seq>
            str, align sequence

        <aligner>
            obj, Bio PairwiseAligner with setting score matrix

    RETURN
        <align_res_list>
            list, each item is also a list like:
                [
                    str.full_ref_seq,
                    str.full_path,
                    str.full_query_seq,
                    int.rel_align_start_idx,
                    int.rel_align_end_idx,
                    int.align_res.score,
                    int.start_idx,
                    int.end_idx,
                    obj.align_res
                ]

    """
    # convert Seq obj
    align_ref_seq = Bio.Seq.Seq(str(align_ref_seq))
    align_query_seq = Bio.Seq.Seq(str(align_query_seq))

    # return list
    align_res_list = []

    # make align steps
    region_idx_list = []
    start_idx = 0
    while start_idx < len(align_ref_seq):
        end_idx = start_idx + align_window_size

        if end_idx <= len(align_ref_seq):
            region_idx_list.append((start_idx, end_idx))
        else:
            region_idx_list.append((start_idx, len(align_ref_seq)))
            break

        start_idx += align_step_size

    # step alignment
    for start_idx, end_idx in region_idx_list:
        region_ref_seq = align_ref_seq[start_idx: end_idx]

        # select the best hit in one alignment
        align_res_dict_list = []
        for step_align_res in aligner.align(region_ref_seq, align_query_seq):
            # this params is designed for SpCas9 seed region alignment
            align_res_dict = degen_align_obj_analysis(step_align_res,
                                                      N0_T_award_score=N0_T_award_score,
                                                      degeneracy_base=degeneracy_base)

            align_res_dict_list.append(align_res_dict)

        # sort
        align_res_dict_list = sorted(align_res_dict_list,
                                     key=lambda i: (i["alignment_score"],
                                                    i["degen_match_count"],
                                                    i["match_count"]),
                                     reverse=True)

        step_best_align_dict = align_res_dict_list[0]

        step_align_res = align_res_dict_list[0]["alignment"]
        align_info_list = str(step_align_res).split("\n")
        align_count = 0

        # fix align path
        fix_align_info_path = align_res_dict_list[0]["alignment_degen_path"]

        path_start_idx = 0
        path_end_idx = len(align_info_list[0])

        path_start_state = False
        path_end_state = False

        for path_idx, path_info in enumerate(align_info_list[1]):
            # query_info = align_info_list[2][path_idx]
            query_info = fix_align_info_path[path_idx]

            if (not path_start_state) and (query_info != "-"):
                path_start_idx = path_idx
                path_start_state = True

            if path_start_state:
                if align_info_list[2][path_idx] != "-":
                    align_count += 1

            if align_count >= len(align_query_seq):
                if path_start_state and (not path_end_state) and query_info == "-":
                    path_end_idx = path_idx
                    path_end_state = True

            if path_start_state and path_end_state:
                break

        # make full length alignment result
        full_ref_seq = align_ref_seq[:start_idx] + "".join(align_info_list[0]) + align_ref_seq[end_idx - 1:]
        # full_path = "-" * start_idx + "".join(align_info_list[1]) + "-" * (len(align_ref_seq) - end_idx)
        full_path = "-" * start_idx + "".join(fix_align_info_path) + "-" * (len(align_ref_seq) - end_idx)
        full_query_seq = "-" * start_idx + "".join(align_info_list[2]) + "-" * (len(align_ref_seq) - end_idx)

        # make relative index
        rel_align_start_idx = path_start_idx + start_idx
        rel_align_end_idx = path_end_idx + start_idx

        # make dict
        step_best_align_dict["align_info_ref"] = str(full_ref_seq)
        step_best_align_dict["align_info_path"] = str(full_path)
        step_best_align_dict["align_info_query"] = str(full_query_seq)
        step_best_align_dict["align_info_start_idx"] = rel_align_start_idx
        step_best_align_dict["align_info_end_idx"] = rel_align_end_idx
        step_best_align_dict["step_start_idx"] = start_idx
        step_best_align_dict["step_end_idx"] = end_idx

        # back list
        align_res_list.append(step_best_align_dict)

    return align_res_list


def get_dist_to_signal(
        signal_chr_start,
        signal_chr_end,
        signal_strand,
        in_align_chr_start,
        in_align_chr_end,
        align_reverse=False):
    """
    INPUT
        <signal_chr_start>
            int,
        <signal_chr_end>
            int,
        <align_chr_start>
            int,
        <align_chr_end>
            int,
        <align_reverse>
            bool, True means aligned to the - strand

    RETURN:
        the distance between aligned TALE 3' site to signal 5' site

    INFO:

        Example1:

            REF - - - - - - - - - - - - 5'C C C C C C C 3' -  -  -  -
            TALE  5'-CGATGGCTATTT-3'
                                  1 2 3   4
            RETURN dist=4

        Example2:

            REF - - - - - - - - - - - - 3'G G G G G G G 5' -  -  -  -
            TALE  5'- T T G T T T T-3'
                               -1 0 1 2 3 4 5 6 7 8 9 10
            RETURN dist=10
    """

    if signal_strand == "+":
        if not align_reverse:
            return signal_chr_start - in_align_chr_end

        else:
            return in_align_chr_start - signal_chr_start

    else:
        if not align_reverse:
            return signal_chr_end - in_align_chr_end

        else:
            return in_align_chr_start - signal_chr_end


def get_align_chrom_coordinate(align_res_dict,
                               ref_chr_start,
                               ref_chr_end,
                               align_ref_start,
                               align_ref_end,
                               reverse_state=False):
    """
    INPUT:
        <align_res_dict>
            dict, format like

            {
                'N0_base': 'T',
                 'align_info_end_idx': 83,
                 'align_info_path': '-----------------------------------------------------------------------.||||.||||||-------------------------------------------------------',
                 'align_info_query': '-----------------------------------------------------------------------CGATGGCTATTT-------------------------------------------------------',
                 'align_info_ref': 'AGTTCAAACGAGGATTCAGGTTTTAGGCCAAATAAGAGATAACATTCTTGAATTCATAATTGCTAGAATCATGATGACTATTTGAAAGGATATATAATTTATCCCTTAATAACTGAAAGGTGTCATCCAAAACTGTTGT',
                 'align_info_start_idx': 71,
                 'alignment': <Bio.Align.PairwiseAlignment object at 0x7fb451d74e20>,
                 'alignment_score': 54.0,
                 'degen_count': 1,
                 'degen_match_count': 11,
                 'degen_match_pos': [5],
                 'degen_mismatch_count': 1,
                 'degen_mismatch_pos': [0],
                 'gap_count': 0,
                 'gap_pos': [],
                 'match_count': 10,
                 'meta': {'N0_T_award_score': 5, 'degeneracy_base': 'AG'},
                 'mismatch_count': 2,
                 'mismatch_pos': [0, 5],
                 'step_end_idx': 85,
                 'step_start_idx': 60
            }

        <ref_chr_start>
            int, UCSC chromosome index of ref sequence start

        <ref_chr_end>
            int, UCSC chromosome index of ref sequence end

        <align_ref_start>
            int, align index on ref seq 0-based

        <align_ref_end>
            int, align index on ref seq 0-based

        <reverse_state>
            bool, if align to the reverse strand

    """

    total_gap_count = align_res_dict["align_info_ref"][align_ref_start:align_ref_end].count("-")

    if not reverse_state:
        back_align_chr_start = ref_chr_start + align_ref_start
        back_align_chr_end = ref_chr_start + align_ref_end - total_gap_count - 1

    else:
        back_align_chr_end = ref_chr_end - align_ref_start
        back_align_chr_start = ref_chr_end - align_ref_end + 1 + total_gap_count

    return back_align_chr_start, back_align_chr_end


def report_TALE_align_res(TALE_align_dict):
    """
    INPUT:
        <TALE_align_dict>
            dict,

            {
                'N0_base': 'T',
                 'align_info_end_idx': 83,
                 'align_info_path': '-----------------------------------------------------------------------.||||.||||||-------------------------------------------------------',
                 'align_info_query': '-----------------------------------------------------------------------CGATGGCTATTT-------------------------------------------------------',
                 'align_info_ref': 'AGTTCAAACGAGGATTCAGGTTTTAGGCCAAATAAGAGATAACATTCTTGAATTCATAATTGCTAGAATCATGATGACTATTTGAAAGGATATATAATTTATCCCTTAATAACTGAAAGGTGTCATCCAAAACTGTTGT',
                 'align_info_start_idx': 71,
                 'alignment': <Bio.Align.PairwiseAlignment object at 0x7fb451d74e20>,
                 'alignment_score': 54.0,
                 'degen_count': 1,
                 'degen_match_count': 11,
                 'degen_match_pos': [5],
                 'degen_mismatch_count': 1,
                 'degen_mismatch_pos': [0],
                 'gap_count': 0,
                 'gap_pos': [],
                 'match_count': 10,
                 'meta': {'N0_T_award_score': 5, 'degeneracy_base': 'AG'},
                 'mismatch_count': 2,
                 'mismatch_pos': [0, 5],
                 'step_end_idx': 85,
                 'step_start_idx': 60
            }

    RETURN:
        <align_out_list>
            list, .art format list
            "chrom",
            "start",
            "end",
            "region_index",
            "align_chr_name",
            "align_chr_start",
            "align_chr_end",
            "align_strand",
            "align_dist_to_signal",
            "align_N0_base",
            "align_total_match",
            "align_total_mismatch",
            "align_degen_total_match",
            "align_degen_total_mismatch",
            "align_degen_num",
            "align_total_gap",
            "align_score",
            "align_target_seq",
            "align_info_state",
            "align_query_seq"
    """

    # load index
    start_idx = TALE_align_dict["align_info_start_idx"]
    end_idx = TALE_align_dict["align_info_end_idx"]

    # make list
    return_out_list = [
        TALE_align_dict["region_info"][0],
        TALE_align_dict["region_info"][1],
        TALE_align_dict["region_info"][2],
        TALE_align_dict["region_info"][3],
        TALE_align_dict["align_chr_name"],
        TALE_align_dict["align_chr_start"],
        TALE_align_dict["align_chr_end"],
        TALE_align_dict["align_chr_strand"],
        TALE_align_dict["align_dist_to_signal"],
        TALE_align_dict["N0_base"],
        TALE_align_dict["match_count"],
        TALE_align_dict["mismatch_count"],
        TALE_align_dict["degen_match_count"],
        TALE_align_dict["degen_mismatch_count"],
        TALE_align_dict["degen_count"],
        TALE_align_dict["gap_count"],
        TALE_align_dict["alignment_score"],
        TALE_align_dict["align_info_ref"][start_idx:end_idx],
        TALE_align_dict["align_info_path"][start_idx:end_idx],
        TALE_align_dict["align_info_query"][start_idx:end_idx]
    ]

    return return_out_list


###############################################################################
# main part
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="From .mpmat to alignment result table (.art) file")

    parser.add_argument("-i", "--input",
                        help=".mpmat table or .bed table, generated from <pmat-merge.py> or <mpmat-select.py>",
                        required=True)

    parser.add_argument("-q", "--query_seq",
                        help="TALE_seq with N0 base info", required=True)

    parser.add_argument("-o", "--output",
                        help="Output alignment result table (.art) filename, default=stdout", default="stdout")

    parser.add_argument("-r", "--reference",
                        help="Reference genome fasta file", required=True)

    parser.add_argument("-m", "--extend_method",
                        help="Select and extend region of mpmat file to get FASTA file, "
                             "<region> <upstream_site> <downstream_site> <highest_site>, "
                             "default=region", default="region")

    parser.add_argument("-e", "--extend_length",
                        help="Region extend length, default=40", default=40, type=int)

    parser.add_argument("--input_filetype",
                        help="'mpmat' or 'bed', default='mpmat'", default="mpmat")

    parser.add_argument("--input_header",
                        help="If .mpmat file contain header, default=False", default="False")

    parser.add_argument("--mpmat_fwd_mut_type",
                        help="mpmat file fwd strand mutation type, default=CT", default="CT")

    parser.add_argument("--mpmat_rev_mut_type",
                        help="mpmat file fwd strand mutation type, default=GA", default="GA")

    parser.add_argument("--bedtools_path",
                        help="Software <bedtools> PATH, default considered <bedtools> is already included in current "
                             "PATH environment.",
                        default="bedtools")

    parser.add_argument("--align_settings",
                        help="Set <align_match_score> <align_mismatch_score> <align_gap_open_score> "
                             "<align_gap_extension_score>, default=5,-4,-24,-8",
                        default="5,-4,-24,-8")

    parser.add_argument("--sub_mat_specific",
                        help="Set specific alignment score, format like ref:query=score, default='A:G=3'",
                        default="A:G=3")

    parser.add_argument("--sub_mat_alphabet",
                        help="substitution matrix alphabet, default=ATCG",
                        default="AGCT")

    parser.add_argument("--N0_T_award",
                        help="TALE prefer 'T' at N0 position, there will be a award score if a 'T' at N0, default=10",
                        default=10, type=int)

    parser.add_argument("--temp_dir",
                        help="default is the same dir to --input", default=None)

    parser.add_argument("--verbose",
                        help="A larger value means more details, default=3", default=3, type=int)

    ARGS = parser.parse_args()

    # ------------------------------------------------------------------------>>>>>>>
    # logging
    # ------------------------------------------------------------------------>>>>>>>
    logging.basicConfig(level=(4 - ARGS.verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ------------------------------------------------------------------------>>>>>>>
    # get temp dir check output file dir
    # ------------------------------------------------------------------------>>>>>>>
    if ARGS.output == "stdout":
        temp_out_dir = os.getcwd()
    else:
        temp_out_dir = os.path.dirname(ARGS.output)

    ARGS.temp_dir = temp_out_dir

    if ARGS.input_header == "True":
        ARGS.input_header = True
    else:
        ARGS.input_header = False

    # ------------------------------------------------------------------------>>>>>>>
    # output full cmd
    # ------------------------------------------------------------------------>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write(_log_cmd_str(ARGS) + "\n")

    # ------------------------------------------------------------------------>>>>>>>
    # make temp fasta
    # ------------------------------------------------------------------------>>>>>>>
    # get fasta temp
    sys.stderr.write("-" * 80 + "\n")
    logging.info("From input file to generate FASTA file...\n")
    temp_fa_state, temp_fa_filename, temp_bed_filename = make_temp_fasta(ARGS.input,
                                                                         ARGS.input_filetype,
                                                                         ARGS.reference,
                                                                         temp_out_dir,
                                                                         mut_type_rev=ARGS.mpmat_rev_mut_type,
                                                                         bedtools_path=ARGS.bedtools_path,
                                                                         header_state=ARGS.input_header,
                                                                         extend_region_size=ARGS.extend_length,
                                                                         extend_region_method=ARGS.extend_method,
                                                                         log_verbose=ARGS.verbose)

    if temp_fa_state != 0:
        logging.error("Make FASTA file unsuccessfully!")
        raise IOError()
    else:
        logging.info("Make FASTA file successfully!")

    # ------------------------------------------------------------------------>>>>>>>
    # make substitution matrix
    # ------------------------------------------------------------------------>>>>>>>
    align_match, align_mismatch, align_gap_open, align_gap_extension = map(int, str(ARGS.align_settings).split(","))

    align_score_matrix = _make_np_substitution_mat(mat_match_score=align_match,
                                                   mat_mismatch_score=align_mismatch,
                                                   alphabet=ARGS.sub_mat_alphabet,
                                                   add_sub_info=ARGS.sub_mat_specific)

    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("Aligner substitution matrix:\n")
    sys.stderr.write(str(align_score_matrix))

    # ------------------------------------------------------------------------>>>>>>>
    # make aligner
    # ------------------------------------------------------------------------>>>>>>>
    # reference Table 6.1: Meta-attributes of the pairwise aligner objects.
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec118
    TALE_aligner = Bio.Align.PairwiseAligner()
    TALE_aligner.substitution_matrix = align_score_matrix
    TALE_aligner.open_gap_score = align_gap_open
    TALE_aligner.extend_gap_score = align_gap_extension
    TALE_aligner.query_left_gap_score = 0
    TALE_aligner.query_right_gap_score = 0
    TALE_aligner.mode = "global"

    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("TALE aligner:" + "\n")
    sys.stderr.write(str(TALE_aligner))
    sys.stderr.write("-" * 80 + "\n")

    # ------------------------------------------------------------------------>>>>>>>
    # output header
    # ------------------------------------------------------------------------>>>>>>>
    out_header_list = [
        "chrom",
        "start",
        "end",
        "region_index",
        "align_chr_name",
        "align_chr_start",
        "align_chr_end",
        "align_strand",
        "align_dist_to_signal",
        "align_N0_base",
        "align_total_match",
        "align_total_mismatch",
        "align_degen_total_match",
        "align_degen_total_mismatch",
        "align_degen_num",
        "align_total_gap",
        "align_score",
        "align_target_seq",
        "align_info_state",
        "align_query_seq"
    ]

    # ------------------------------------------------------------------------>>>>>>>
    # open files
    # ------------------------------------------------------------------------>>>>>>>
    if ARGS.output == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(ARGS.output, "w")

    region_file = open(ARGS.input, "r")
    temp_bed = open(temp_bed_filename, "r")
    temp_fa = Bio.SeqIO.parse(temp_fa_filename, format="fasta")

    # output header line
    out_file.write("\t".join(out_header_list) + "\n")

    # ------------------------------------------------------------------------>>>>>>>
    # make query sequence
    # ------------------------------------------------------------------------>>>>>>>
    query_seq = Bio.Seq.Seq(str(ARGS.query_seq).upper())

    # ------------------------------------------------------------------------>>>>>>>
    # alignment
    # ------------------------------------------------------------------------>>>>>>>
    if ARGS.input_header:
        region_file.readline()

    for fa_idx, fa_info in enumerate(temp_fa):
        bed_list = temp_bed.readline().strip().split("\t")
        region_list = region_file.readline().strip().split("\t")

        # load genome coordinate of bedtools getfasta
        bed_chr_start = int(bed_list[1]) + 1
        bed_chr_end = int(bed_list[2])

        # load signal region coordinate
        region_chr_start = int(region_list[1])
        region_chr_end = int(region_list[2])
        region_signal_strand = bed_list[5]

        # load fasta seq
        fasta_seq_fwd = fa_info.seq.upper()
        fasta_seq_rev = Bio.Seq.reverse_complement(fasta_seq_fwd)

        # step alignment
        step_window_size = 50
        step_move_size = 10

        step_align_res_fwd = run_step_align(align_ref_seq=fasta_seq_fwd,
                                            align_query_seq=query_seq,
                                            aligner=TALE_aligner,
                                            align_window_size=step_window_size,
                                            align_step_size=step_move_size,
                                            N0_T_award_score=ARGS.N0_T_award,
                                            degeneracy_base="AG")

        step_align_res_rev = run_step_align(align_ref_seq=fasta_seq_rev,
                                            align_query_seq=query_seq,
                                            aligner=TALE_aligner,
                                            align_window_size=step_window_size,
                                            align_step_size=step_move_size,
                                            N0_T_award_score=ARGS.N0_T_award,
                                            degeneracy_base="AG")

        # select the best for each strand
        best_align_dict_fwd = sorted(step_align_res_fwd,
                                     key=lambda i: (i["alignment_score"],
                                                    i["degen_match_count"],
                                                    i["match_count"]),
                                     reverse=True)[0]

        best_align_dict_rev = sorted(step_align_res_rev,
                                     key=lambda i: (i["alignment_score"],
                                                    i["degen_match_count"],
                                                    i["match_count"]),
                                     reverse=True)[0]

        # select the best
        best_align_dict = None
        best_align_reverse_state = None

        if best_align_dict_fwd["alignment_score"] > best_align_dict_rev["alignment_score"]:
            best_align_dict = best_align_dict_fwd
            best_align_reverse_state = False

        else:
            best_align_dict = best_align_dict_rev
            best_align_reverse_state = True

        # get genome coordinate
        align_chr_start, align_chr_end = get_align_chrom_coordinate(best_align_dict,
                                                                    ref_chr_start=bed_chr_start,
                                                                    ref_chr_end=bed_chr_end,
                                                                    align_ref_start=best_align_dict[
                                                                        "align_info_start_idx"],
                                                                    align_ref_end=best_align_dict["align_info_end_idx"],
                                                                    reverse_state=best_align_reverse_state)

        # calculate signal to dist penalty
        align_dist_to_signal = get_dist_to_signal(signal_chr_start=region_chr_start,
                                                  signal_chr_end=region_chr_end,
                                                  signal_strand=region_signal_strand,
                                                  in_align_chr_start=align_chr_start,
                                                  in_align_chr_end=align_chr_end,
                                                  align_reverse=best_align_reverse_state)

        # add info into dict
        best_align_dict["align_chr_name"] = bed_list[0]
        best_align_dict["align_chr_start"] = align_chr_start
        best_align_dict["align_chr_end"] = align_chr_end
        best_align_dict["align_dist_to_signal"] = align_dist_to_signal
        best_align_dict["region_info"] = bed_list[:4]

        if not best_align_reverse_state:
            best_align_dict["align_chr_strand"] = "+"
        else:
            best_align_dict["align_chr_strand"] = "-"

        out_list = report_TALE_align_res(best_align_dict)

        # output
        out_file.write("\t".join(list(map(str, out_list))) + "\n")

    # ------------------------------------------------------------------------>>>>>>>
    # close and  remove temp files
    # ------------------------------------------------------------------------>>>>>>>
    if ARGS.output != "stdout":
        out_file.close()

    temp_bed.close()
    region_file.close()

    # remove temp file
    os.remove(temp_bed_filename)
    os.remove(temp_fa_filename)

    # ------------------------------------------------------------------------>>>>>>>
    # final log
    # ------------------------------------------------------------------------>>>>>>>
    logging.info("Done!")

# 2022-02-13 by Meng Haowei
