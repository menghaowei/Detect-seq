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

# Version information START --------------------------------------------------
from typing import Dict, Union, Any

VERSION_INFO = """
Author: MENG Howard

Version-07:
    For Cas9 based sgRNA alignment

Version-06:
    2022-03-31
        fix bug with penalty params

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
    --find_highest_signal_method {find_highest_signal_method}
    --input_filetype {input_filetype}
    --input_header {input_header}
    --mpmat_fwd_mut_type {mpmat_fwd_mut_type}
    --mpmat_rev_mut_type {mpmat_rev_mut_type}    
    --PAM_index {PAM_index}
    --seed_index {seed_index}
    --distal_index {distal_index}
    --bedtools_path {bedtools_path}
    --align_settings {align_settings}
    --align_step_window_size {align_step_window_size}
    --align_step_move_size {align_step_move_size}
    --sub_mat_alphabet {sub_mat_alphabet}
    --sub_mat_specific {sub_mat_specific}
    --dna_bulge_penalty {dna_bulge_penalty}
    --rna_bulge_penalty {rna_bulge_penalty}
    --PAM_type_penalty {PAM_type_penalty}
    --dna_bulge_cmp_weight {dna_bulge_cmp_weight}
    --rna_bulge_cmp_weight {rna_bulge_cmp_weight}
    --mismatch_cmp_weight {mismatch_cmp_weight}
    --dist_to_signal_penalty_breaks {dist_to_signal_penalty_breaks}
    --dist_to_signal_penalty_k {dist_to_signal_penalty_k}
    --dist_to_signal_penalty_offset {dist_to_signal_penalty_offset}
    --temp_dir {temp_dir}
    --verbose {verbose}""".format(
        input=args.input,
        query_seq=args.query_seq,
        output=args.output,
        reference=args.reference,
        extend_method=args.extend_method,
        extend_length=args.extend_length,
        find_highest_signal_method=args.find_highest_signal_method,
        input_filetype=args.input_filetype,
        input_header=args.input_header,
        mpmat_fwd_mut_type=args.mpmat_fwd_mut_type,
        mpmat_rev_mut_type=args.mpmat_rev_mut_type,
        PAM_index=args.PAM_index,
        seed_index=args.seed_index,
        distal_index=args.distal_index,
        bedtools_path=args.bedtools_path,
        align_settings=args.align_settings,
        align_step_window_size=args.align_step_window_size,
        align_step_move_size=args.align_step_move_size,
        sub_mat_alphabet=args.sub_mat_alphabet,
        sub_mat_specific=args.sub_mat_specific,
        dna_bulge_penalty=args.dna_bulge_penalty,
        rna_bulge_penalty=args.rna_bulge_penalty,
        PAM_type_penalty=args.PAM_type_penalty,
        dna_bulge_cmp_weight=args.dna_bulge_cmp_weight,
        rna_bulge_cmp_weight=args.rna_bulge_cmp_weight,
        mismatch_cmp_weight=args.mismatch_cmp_weight,
        dist_to_signal_penalty_breaks=args.dist_to_signal_penalty_breaks,
        dist_to_signal_penalty_k=args.dist_to_signal_penalty_k,
        dist_to_signal_penalty_offset=args.dist_to_signal_penalty_offset,
        temp_dir=args.temp_dir,
        verbose=args.verbose
    )

    return full_cmd_str


def find_highest_mut_site(mpmat_line, find_method="score"):
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
    mpmat_line_list = mpmat_line.strip().split("\t")

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
                    header_state="False",
                    extend_region_size=20,
                    extend_region_method="highest_site",
                    highest_method="score",
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
    temp_bed_filename = os.path.abspath(os.path.join(temp_dir, temp_bed_basename))
    temp_bed_file = open(temp_bed_filename, "w")

    with open(input_filename, "r") as in_file:
        if header_state != "False":
            in_file.readline()

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
                    find_highest_res = find_highest_mut_site(mpmat_line=line, find_method=highest_method)
                    site_chr_pos = find_highest_res[1]
                    region_extend_start = site_chr_pos - extend_region_size
                    region_extend_end = site_chr_pos + extend_region_size

                else:
                    region_extend_start = region_start - extend_region_size
                    region_extend_end = region_end + extend_region_size

            else:
                raise KeyError("Wrong extend method!")

            # output info into bed
            bed_out_list = [region_chr, region_extend_start, region_extend_end, region_index]
            temp_bed_file.write("\t".join(map(str, bed_out_list)) + "\n")

    temp_bed_file.close()

    # log
    logging.info("Output temp bed with filename: %s\n" % temp_bed_filename)

    # use bedtools to make fasta file
    temp_fasta_filename = temp_bed_filename + ".fa"

    tool_bedtools_cmd_fmt = "{bedtools} getfasta -fi {ref_fasta}  -fo {out_fasta} -bed {input_bed} -name+"
    tool_bedtools_cmd = tool_bedtools_cmd_fmt.format(
        bedtools=bedtools_path,
        ref_fasta=genome_fasta,
        out_fasta=temp_fasta_filename,
        input_bed=temp_bed_filename
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
        os.remove(temp_bed_filename)
        logging.error("Error! Unsuccessfully get region FASTA seq!")
        return 1, None, None

    else:
        logging.info("Successfully get region FASTA seq!\n")
        return 0, temp_fasta_filename, temp_bed_filename


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


def weight_align_obj_analysis(alignment,
                              gap_weight_coef=10,
                              gap_start_score=1,
                              mis_weight_coef=1,
                              mis_start_score=5):
    """
    INPUT:
        <alignment>
            obj, alignment result

        <*_weight_coef>
            int, end_score = start score * weight_coef;

        <*_start_score>
            int, start score

    RETURN:
        <info_dict>
            int, gap count
            int, match count
            int, mismatch count
            list, gap position [0 based index]
            list, mismatch position [0 based index]
            float, gap penalty score
            float, mismatch penalty score
            obj, alignment obj

    INFO:
        The gap count should be the num of gap contain in sgRNA alignment region.

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
            "gap_penalty_score" : 1 + 10 * (10 / 13) + 1 + 10 * (11 / 13),
            "mismatch_penalty_score" : 1 + 10 * (8 / 13),
            "alignment" : alignment obj,
            "meta": {
                "gap_weight_coef" : gap_weight_coef,
                "gap_start_score" : gap_start_score,
                "mis_weight_coef" : mis_weight_coef,
                "mis_start_score" : mis_start_score
            }
        }
    """

    # init dict
    info_dict = {
        "gap_count": 0,
        "match_count": 0,
        "mismatch_count": 0,
        "gap_pos": [],
        "mismatch_pos": [],
        "gap_penalty_score": 0,
        "mismatch_penalty_score": 0,
        "alignment": alignment,
        "meta": {
            "gap_weight_coef": gap_weight_coef,
            "gap_start_score": gap_start_score,
            "mis_weight_coef": mis_weight_coef,
            "mis_start_score": mis_start_score
        }
    }

    # fetch alignment info
    alignment_list = str(alignment).split("\n")
    query_length = len(alignment.query)

    # load alignment
    info_list = alignment_list[1]
    query_list = alignment_list[2]

    # count init value
    count_state = False
    count_query_base = 0

    query_index = 0

    # counting
    for run_index, info_base in enumerate(info_list):
        if (not count_state) and (query_list[run_index] != "-"):
            count_state = True

        if count_state:
            if info_base == "|":
                info_dict["match_count"] += 1
                query_index += 1

            elif info_base == "X" or info_base == ".":
                info_dict["mismatch_count"] += 1
                info_dict["mismatch_pos"].append(query_index)
                query_index += 1

            elif info_base == "-":
                info_dict["gap_count"] += 1
                info_dict["gap_pos"].append(query_index)

                if query_list[run_index] != "-":
                    query_index += 1

        if query_index >= query_length:
            break

    # count gap score
    step_len = len(alignment.query) - 1
    gap_end_score = gap_start_score * gap_weight_coef
    gap_step_score = (gap_end_score - gap_start_score) / 1.0 / step_len

    if info_dict["gap_count"] > 0:
        for gap_pos in info_dict["gap_pos"]:
            info_dict["gap_penalty_score"] += gap_start_score + gap_pos * gap_step_score

    # count mismatch score
    mis_end_score = mis_start_score * mis_weight_coef
    mis_step_score = (mis_end_score - mis_start_score) / 1.0 / step_len

    if info_dict["mismatch_count"] > 0:
        for mis_pos in info_dict["mismatch_pos"]:
            info_dict["mismatch_penalty_score"] += mis_start_score + mis_pos * mis_step_score

    return info_dict


def run_step_align(align_ref_seq, align_query_seq, aligner,
                   align_window_size=23,
                   align_step_size=12, **args):
    """
    INPUT:
        <align_ref_seq>
            str, reference sequence

        <align_query_seq>
            str, align sequence, support "N" and "R" degeneracy, N="A/T/C/G" R="A/G"

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
    # params default
    if args.get("gap_weight_coef") is None:
        align_gap_weight_coef = 1
    else:
        align_gap_weight_coef = int(args.get("gap_weight_coef"))

    if args.get("gap_start_score") is None:
        align_gap_start_score = 24
    else:
        align_gap_start_score = int(args.get("gap_start_score"))

    if args.get("mis_weight_coef") is None:
        align_mis_weight_coef = 1
    else:
        align_mis_weight_coef = int(args.get("mis_weight_coef"))

    if args.get("mis_start_score") is None:
        align_mis_start_score = 5
    else:
        align_mis_start_score = int(args.get("mis_start_score"))

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
            align_res_dict = weight_align_obj_analysis(step_align_res,
                                                       gap_weight_coef=align_gap_weight_coef,
                                                       gap_start_score=align_gap_start_score,
                                                       mis_weight_coef=align_mis_weight_coef,
                                                       mis_start_score=align_mis_start_score)

            align_res_dict["total_penalty_score"] = align_res_dict["gap_penalty_score"] + \
                                                    align_res_dict["mismatch_penalty_score"]

            align_res_dict_list.append(align_res_dict)

        # sort
        align_res_dict_list = sorted(align_res_dict_list,
                                     key=lambda i: (i["total_penalty_score"],
                                                    i["gap_penalty_score"],
                                                    i["mismatch_penalty_score"])
                                     )

        step_align_res = align_res_dict_list[0]["alignment"]
        align_info_list = str(step_align_res).split("\n")
        align_count = 0

        path_start_idx = 0
        path_end_idx = len(align_info_list[0])

        path_start_state = False
        path_end_state = False

        for path_idx, path_info in enumerate(align_info_list[1]):
            query_info = align_info_list[2][path_idx]

            if (not path_start_state) and (query_info != "-"):
                path_start_idx = path_idx
                path_start_state = True

            if path_start_state:
                if query_info != "-":
                    align_count += 1

            if align_count >= len(align_query_seq):
                if path_start_state and (not path_end_state) and query_info == "-":
                    path_end_idx = path_idx
                    path_end_state = True

            if path_start_state and path_end_state:
                break

        # make full length alignment result
        full_ref_seq = align_ref_seq[:start_idx] + "".join(align_info_list[0]) + align_ref_seq[end_idx - 1:]
        full_path = "-" * start_idx + "".join(align_info_list[1]) + "-" * (len(align_ref_seq) - end_idx)
        full_query_seq = "-" * start_idx + "".join(align_info_list[2]) + "-" * (len(align_ref_seq) - end_idx)

        # make relative index
        rel_align_start_idx = path_start_idx + start_idx
        rel_align_end_idx = path_end_idx + start_idx

        # back list
        align_res_list.append(
            [
                str(full_ref_seq),
                str(full_path),
                str(full_query_seq),
                rel_align_start_idx,
                rel_align_end_idx,
                step_align_res.score,
                start_idx,
                end_idx,
                step_align_res
            ])

    return align_res_list


def weight_extend_align_obj_analysis(extend_align_res,
                                     align_query_seq,
                                     pam_pos_idx=(21, 23),
                                     dna_bulge_coef=10,
                                     dna_bulge_start_score=1,
                                     rna_bulge_coef=20,
                                     rna_bulge_start_score=2,
                                     mis_coef=10,
                                     mis_start_score=1,
                                     pam_type_penalty=(0, 25, 50)
                                     ):
    """
    INPUT:
        <extend_align_res>
            list, info like:
                [
                    str(full_align_ref_seq),
                    str(full_align_path),
                    str(full_align_query_seq),
                    rel_align_start_idx,
                    rel_align_end_idx
                ]

        <pam_pos_idx>
            tuple, 1-based index of PAM position, for SpCas9 is (21,23)

        <*_coef>
            int, end_score = start score * coef;

        <*_start_score>
            int, start score

        <pam_type_penalty>
            tuple, default=(0,25,50) means the penalty score for NGG, NAG and NotNRG

    RETURN:
        <info_dict>
            int, match count
            int, mismatch count
            int, DNA bulge count
            int, RNA bulge count
            int, total gap count
            list, mismatch position [1 based index]
            list, DNA bulge position [1 based index], 1 means after the 1st base
            list, RNA bulge position [1 based index], 1 means after the 1st base
            str, PAM type (NGG, NAG, NoPAM)
            str, PAM seq
            list, alignment info
            int, target start index [0 based]
            int, target end index [0 based]
            float, dna bulge penalty score
            float, rna bulge penalty score
            float, gap penalty score
            float, mismatch penalty score
            float, PAM type penalty score

    INFO:
        The gap count should be the num of gap contain in sgRNA alignment region.

        e.g.

        DNA bulge = 2, RNA bulge = 0;

        AGTGGTAGACATT
        ------||X|--|
        ------AGCC--T

        DNA bulge = 2, RNA bulge = 1;
        AGTGGTA-ACATT
        ------|-X|--|
        ------AGCC--T

        return = {

        }
    """

    # init dict
    info_dict = {
        "match_count": 0,
        "mismatch_count": 0,
        "dna_bulge_count": 0,
        "rna_bulge_count": 0,
        "gap_count": 0,
        "mismatch_pos": [],
        "dna_bulge_pos": [],
        "rna_bulge_pos": [],
        "PAM_type": None,
        "PAM_seq": "",
        "alignment_info": extend_align_res[:3],
        "target_start_index": extend_align_res[3],
        "target_end_index": extend_align_res[4],
        "dna_bulge_penalty_score": 0,
        "rna_bulge_penalty_score": 0,
        "gap_penalty_score": 0,
        "mismatch_penalty_score": 0,
        "PAM_type_penalty_score": 0
    }

    # fetch alignment info
    target_list = extend_align_res[0]
    info_list = extend_align_res[1]
    query_list = extend_align_res[2]

    # length
    query_length = len(align_query_seq)

    # count init value
    count_state = False

    # index init
    query_index = 0

    # counting
    for run_index, info_base in enumerate(info_list):
        if (not count_state) and (query_list[run_index] != "-"):
            count_state = True

        if count_state:
            # count PAM
            if pam_pos_idx[0] <= (query_index + 1) <= pam_pos_idx[1]:
                info_dict["PAM_seq"] += target_list[run_index]
                query_index += 1

            else:
                if info_base == "|":
                    info_dict["match_count"] += 1
                    query_index += 1

                elif info_base == "X" or info_base == ".":
                    info_dict["mismatch_count"] += 1
                    info_dict["mismatch_pos"].append(query_index + 1)
                    query_index += 1

                elif info_base == "-":
                    info_dict["gap_count"] += 1

                    if target_list[run_index] == "-":
                        info_dict["rna_bulge_count"] += 1
                        info_dict["rna_bulge_pos"].append(query_index)
                        query_index += 1
                    else:
                        info_dict["dna_bulge_count"] += 1
                        info_dict["dna_bulge_pos"].append(query_index)

        if query_index >= query_length:
            break

    # PAM
    if info_dict["PAM_seq"][-2:] == "AG":
        info_dict["PAM_type"] = "NAG"
        info_dict["PAM_type_penalty_score"] = pam_type_penalty[1]

    elif info_dict["PAM_seq"][-2:] == "GG":
        info_dict["PAM_type"] = "NGG"
        info_dict["PAM_type_penalty_score"] = pam_type_penalty[0]

    else:
        info_dict["PAM_type"] = "NotNRG"
        info_dict["PAM_type_penalty_score"] = pam_type_penalty[2]

        # count step score
    step_len = query_length - 1

    # score dict
    score_dict = {
        "dna_bulge_start_score": dna_bulge_start_score,
        "dna_bulge_step_score": (dna_bulge_start_score * dna_bulge_coef - dna_bulge_start_score) / 1.0 / step_len,
        "rna_bulge_start_score": rna_bulge_start_score,
        "rna_bulge_step_score": (rna_bulge_start_score * rna_bulge_coef - rna_bulge_start_score) / 1.0 / step_len,
        "mismatch_start_score": mis_start_score,
        "mismatch_step_score": (mis_start_score * mis_coef - mis_start_score) / 1.0 / step_len,
    }

    # count item
    score_item_list = ["mismatch", "dna_bulge", "rna_bulge"]

    # calculate penalty score
    for score_item in score_item_list:
        key_count = score_item + "_count"
        key_pos = score_item + "_pos"
        key_score = score_item + "_penalty_score"
        key_start = score_item + "_start_score"
        key_step = score_item + "_step_score"

        if info_dict[key_count] > 0:
            for pos_idx in info_dict[key_pos]:
                info_dict[key_score] += score_dict[key_start] + pos_idx * score_dict[key_step]

    # total gap penalty
    info_dict["gap_penalty_score"] = info_dict["dna_bulge_penalty_score"] + info_dict["rna_bulge_penalty_score"]

    return info_dict


def run_extend_align(step_align_res_list,
                     align_ref_seq,
                     seed_query_seq,
                     extend_query_seq,
                     run_extend_aligner,
                     extend_length=18,
                     extend_direction="left",
                     **args):
    """
    INPUT:
        <step_align_res_list>
            list, generated from FUN <step_align>

        <align_ref_seq>
            str, full length reference sequence

        <seed_query_seq>
            str, support N and R degeneracy.

        <extend_query_seq>
            str, support N and R degeneracy.
                <extend_query_seq> + <seed_query_seq> = <full_query_seq>

        <extend_length>
            int, Default=18, suitable for SpCas9;
                 A suggestion extend_length = 5 + len(extend_query_seq)

        <extend_aligner>
            aligner, obj

    FUTURE:
        <extend_direction>

    RETURN:
        <full_query_align_list>
            list, each item is also a list like:
            [
                str(full_align_ref_seq),
                str(full_align_path),
                str(full_align_query_seq),
                rel_align_start_idx,
                rel_align_end_idx
            ]

            Fetch alignment info:
                print(full_align_ref_seq[rel_align_start_idx:rel_align_end_idx])
                print(full_align_path[rel_align_start_idx:rel_align_end_idx])
                print(full_align_query_seq[rel_align_start_idx:rel_align_end_idx])
    """
    # set default params
    # params default
    if args.get("gap_weight_coef") is None:
        align_gap_weight_coef = 1
    else:
        align_gap_weight_coef = int(args.get("gap_weight_coef"))

    if args.get("gap_start_score") is None:
        align_gap_start_score = 24
    else:
        align_gap_start_score = int(args.get("gap_start_score"))

    if args.get("mis_weight_coef") is None:
        align_mis_weight_coef = 1
    else:
        align_mis_weight_coef = int(args.get("mis_weight_coef"))

    if args.get("mis_start_score") is None:
        align_mis_start_score = 5
    else:
        align_mis_start_score = int(args.get("mis_start_score"))

    # init vars
    full_query_seq = Bio.Seq.Seq(str(extend_query_seq + seed_query_seq))
    align_ref_seq = Bio.Seq.Seq(str(align_ref_seq))
    extend_length = max(extend_length, len(seed_query_seq) + 5)

    # make aligner
    run_extend_aligner = run_extend_aligner

    # store align results
    full_query_align_list = []

    for step_align_res in step_align_res_list:
        align_start_idx = step_align_res[3]
        align_end_idx = step_align_res[4]

        if extend_direction == "left":
            if align_start_idx < len(extend_query_seq):
                continue

            elif align_start_idx < extend_length:
                fetch_start_idx = 0

            else:
                fetch_start_idx = align_start_idx - extend_length

            fetch_align_seq = align_ref_seq[fetch_start_idx:align_end_idx]

            # sort align result
            extend_align_dict_list = []
            for extend_align_res in run_extend_aligner.align(fetch_align_seq, full_query_seq):
                extend_align_dict = weight_align_obj_analysis(extend_align_res,
                                                              gap_weight_coef=align_gap_weight_coef,
                                                              gap_start_score=align_gap_start_score,
                                                              mis_weight_coef=align_mis_weight_coef,
                                                              mis_start_score=align_mis_start_score)

                extend_align_dict["total_penalty_score"] = extend_align_dict["gap_penalty_score"] + \
                                                           extend_align_dict["mismatch_penalty_score"]

                extend_align_dict_list.append(extend_align_dict)

            # select the best hit in one alignment
            extend_align_dict_list = sorted(extend_align_dict_list,
                                            key=lambda i: (i["total_penalty_score"],
                                                           i["gap_count"],
                                                           i["gap_penalty_score"],
                                                           i["mismatch_penalty_score"]))

            extend_align_res_select = extend_align_dict_list[0]["alignment"]

            # get reference pos index
            # analysis alignment result
            align_info_list = str(extend_align_res_select).split("\n")
            align_count = 0

            path_start_idx = 0
            path_end_idx = len(align_info_list[0])

            path_start_state = False
            path_end_state = False

            for path_idx, path_info in enumerate(align_info_list[1]):
                query_info = align_info_list[2][path_idx]

                if (not path_start_state) and (query_info != "-"):
                    path_start_idx = path_idx
                    path_start_state = True

                if path_start_state:
                    if query_info != "-":
                        align_count += 1

                if align_count >= len(full_query_seq):
                    if path_start_state and (not path_end_state) and query_info == "-":
                        path_end_idx = path_idx
                        path_end_state = True

                if path_start_state and path_end_state:
                    break

            # make full length alignment result
            full_align_ref_seq = align_ref_seq[:fetch_start_idx] + \
                                 "".join(align_info_list[0]) + \
                                 align_ref_seq[align_end_idx - 1:]

            full_align_path = "-" * fetch_start_idx + \
                              "".join(align_info_list[1]) + \
                              "-" * (len(align_ref_seq) - align_end_idx)

            full_align_query_seq = "-" * fetch_start_idx + \
                                   "".join(align_info_list[2]) + \
                                   "-" * (len(align_ref_seq) - align_end_idx)

            # make relative index
            rel_align_start_idx = path_start_idx + fetch_start_idx
            rel_align_end_idx = path_end_idx + fetch_start_idx

            # back list
            full_query_align_list.append(
                [
                    str(full_align_ref_seq),
                    str(full_align_path),
                    str(full_align_query_seq),
                    rel_align_start_idx,
                    rel_align_end_idx
                ])

    return full_query_align_list


def get_align_chrom_coordinate(align_res,
                               ref_chr_start,
                               ref_chr_end,
                               align_ref_start,
                               align_ref_end,
                               reverse_state=False):
    """
    INPUT:
        <align_res>
            list,
            [
                str(full_align_ref_seq),
                str(full_align_path),
                str(full_align_query_seq),
                rel_align_start_idx,
                rel_align_end_idx
            ]

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
    total_gap_count = align_res[0][align_ref_start:align_ref_end].count("-")

    if not reverse_state:
        return_align_chr_start = ref_chr_start + align_ref_start
        return_align_chr_end = ref_chr_start + align_ref_end - total_gap_count - 1

    else:
        return_align_chr_end = ref_chr_end - align_ref_start
        return_align_chr_start = ref_chr_end - align_ref_end + 1 + total_gap_count

    return return_align_chr_start, return_align_chr_end


def get_dist_to_signal(
        signal_chr_start,
        signal_chr_end,
        align_chr_start,
        align_chr_end,
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
        alignment dist to signal

    INFO:
        - - - GGCACTGCGGNRG  -  -  -  -
        3 2 1 0           0 -1 -2 -3 -4

    """

    if not align_reverse:
        if align_chr_start > signal_chr_end:
            return signal_chr_end - align_chr_start

        elif align_chr_start < signal_chr_start:
            return signal_chr_start - align_chr_start

        else:
            return 0

    else:
        if align_chr_end > signal_chr_end:
            return align_chr_end - signal_chr_end

        elif align_chr_end < signal_chr_start:
            return align_chr_end - signal_chr_start

        else:
            return 0


def signal_dist_penalty(dist_to_signal,
                        region_break_list=(-10, -5, 0, 8, 15),
                        region_k_list=(-4, -2, 0, 0, 2, 4),
                        region_offset_list=(0, 0, 0, 0, 0, 0)):
    """
    INPUT
        <dist_to_signal>
            int, negative num means signal at upstream;
                 positive num means signal at downstream

        <region_break_list>
            list or tuple, n breaks will create n+1 regions from -Inf to +Inf

        <region_k_list>
            list or tuple, each k match relative region and penalty score in such region = k * dist length

        <region_offset_list>
            list or tuple, each offset match relative region

        penalty score = k * dist length + offset

    RETURN
        <penalty_score>

    INFO
        This function can calculate penalty score related to the distance bewteen signal and alignment.

    """
    region_count = len(region_break_list) + 1
    region_len_list = [0] * region_count

    # get region 0 start
    region_zero_idx = -1
    region_dist_idx = -1

    for region_idx, region_break in enumerate(region_break_list):
        if region_zero_idx == -1:
            if 0 < region_break:
                region_zero_idx = region_idx

        if region_dist_idx == -1:
            if dist_to_signal < region_break:
                region_dist_idx = region_idx

        if region_idx > 0:
            region_len_list[region_idx] = region_break - region_break_list[region_idx - 1]

    if region_zero_idx == -1:
        region_zero_idx = region_count - 1

    if region_dist_idx == -1:
        region_dist_idx = region_count - 1

    penalty_score = 0
    accum_dist = 0
    if region_dist_idx >= region_zero_idx:
        for run_idx in range(region_zero_idx, region_dist_idx):
            penalty_score += region_k_list[run_idx] * region_len_list[run_idx] + region_offset_list[run_idx]
            accum_dist += region_len_list[run_idx]

        penalty_score += region_k_list[region_dist_idx] * (dist_to_signal - accum_dist) + region_offset_list[
            region_dist_idx]

    else:
        for run_idx in range(region_dist_idx + 1, region_zero_idx):
            penalty_score += -1 * region_k_list[run_idx] * region_len_list[run_idx] + region_offset_list[run_idx]
            accum_dist += region_len_list[run_idx]

        penalty_score += -1 * region_k_list[region_dist_idx] * (accum_dist - dist_to_signal) + region_offset_list[
            region_dist_idx]

    return penalty_score


def report_align_res(input_align_dict,
                     pam_proximal_idx=(13, 20),
                     pam_distal_idx=(1, 8)):
    """
    INPUT:
        <align_dict>
            dict, generated from sorted FUN <weight_extend_align_obj_analysis>

        <pam_proximal_idx>
            tuple, 1-based index of query_seq, (13,20) for SpCas9

        <pam_distal_idx>
            tuple, 1-based index of query_seq, (1,8) for SpCas9

    RETURN:
        <align_out_list>
            list, .art format list

            # "chrom",
            # "start",
            # "end",
            # "region_index",
            # "align_chr_name",
            # "align_chr_start",
            # "align_chr_end",
            # "align_strand",
            # "align dist to region"
            # "PAM_type",
            # "PAM_seq",
            # "align_total_match",
            # "align_total_mismatch",
            # "align_total_gap",
            # "dna_bulge_num",
            # "rna_bulge_num",
            # "head_mismatch",
            # "head_gap",
            # "seed_mismatch",
            # "seed_gap",
            # "align_penalty_score",
            # "align_target_seq",
            # "align_info_state",
            # "align_query_seq"

    """

    # count
    head_gap_count = 0
    head_mis_count = 0
    seed_gap_count = 0
    seed_mis_count = 0

    if input_align_dict["dna_bulge_count"] != 0:
        for bulge_idx in input_align_dict["dna_bulge_pos"]:
            if pam_distal_idx[0] <= (bulge_idx + 1) <= pam_distal_idx[1]:
                head_gap_count += 1

            if pam_proximal_idx[0] <= (bulge_idx + 1) <= pam_proximal_idx[1]:
                seed_gap_count += 1

    if input_align_dict["rna_bulge_count"] != 0:
        for bulge_idx in input_align_dict["rna_bulge_pos"]:
            if pam_distal_idx[0] <= (bulge_idx + 1) <= pam_distal_idx[1]:
                head_gap_count += 1

            if pam_proximal_idx[0] <= (bulge_idx + 1) <= pam_proximal_idx[1]:
                seed_gap_count += 1

    if input_align_dict["mismatch_count"] != 0:
        for mis_idx in input_align_dict["mismatch_pos"]:
            if pam_distal_idx[0] <= (mis_idx + 1) <= pam_distal_idx[1]:
                head_mis_count += 1

            if pam_proximal_idx[0] <= (mis_idx + 1) <= pam_proximal_idx[1]:
                seed_mis_count += 1

    # make output list
    start_idx = input_align_dict["target_start_index"]
    end_idx = input_align_dict["target_end_index"]

    return_out_list = [
        input_align_dict["region_chr_name"],
        input_align_dict["region_chr_start"],
        input_align_dict["region_chr_end"],
        input_align_dict["region_index"],
        input_align_dict["region_chr_name"],
        input_align_dict["align_chr_start"],
        input_align_dict["align_chr_end"],
        input_align_dict["align_strand"],
        input_align_dict["dist_align_to_signal"],
        input_align_dict["PAM_type"],
        input_align_dict["PAM_seq"],
        input_align_dict["match_count"],
        input_align_dict["mismatch_count"],
        input_align_dict["gap_count"],
        input_align_dict["dna_bulge_count"],
        input_align_dict["rna_bulge_count"],
        head_mis_count,
        head_gap_count,
        seed_mis_count,
        seed_gap_count,
        round(input_align_dict["total_penalty"], 5),
        input_align_dict["alignment_info"][0][start_idx:end_idx],
        input_align_dict["alignment_info"][1][start_idx:end_idx],
        input_align_dict["alignment_info"][2][start_idx:end_idx]
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
                        help="sgRNA sequence with PAM (NGG/NAG) sequence, support NGG, NRG", required=True)

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

    parser.add_argument("--find_highest_signal_method",
                        help="Use 'ratio', 'mut_count' or 'score' to find highest signal"
                             "default=score", default="score")

    parser.add_argument("--input_filetype",
                        help="'mpmat' or 'bed', default='mpmat'", default="mpmat")

    parser.add_argument("--input_header",
                        help="If .mpmat file contain header, default=False", default="False")

    parser.add_argument("--mpmat_fwd_mut_type",
                        help="mpmat file fwd strand mutation type, default=AG", default="AG")

    parser.add_argument("--mpmat_rev_mut_type",
                        help="mpmat file fwd strand mutation type, default=TC", default="TC")

    parser.add_argument("--PAM_index",
                        help="PAM index, for SpCas9 is 21,23", default="21,23")

    parser.add_argument("--seed_index",
                        help="seed index, for SpCas9 is 13,20", default="13,20")

    parser.add_argument("--distal_index",
                        help="distal index, for SpCas9 is 1,8", default="1,8")

    parser.add_argument("--bedtools_path",
                        help="Software <bedtools> PATH, default considered <bedtools> is already included in current "
                             "PATH environment.",
                        default="bedtools")

    parser.add_argument("--align_settings",
                        help="Set <align_match_score> <align_mismatch_score> <align_gap_open_score> "
                             "<align_gap_extension_score>, default=5,-4,-12,-8",
                        default="5,-4,-12,-8")

    parser.add_argument("--align_step_window_size",
                        help="Step alignment window size, default=30", type=int,
                        default=30)

    parser.add_argument("--align_step_move_size",
                        help="Step alignment move size, default=10", type=int,
                        default=10)

    parser.add_argument("--sub_mat_alphabet",
                        help="substitution matrix alphabet, default=ATCGRN",
                        default="AGCTRN")

    parser.add_argument("--sub_mat_specific",
                        help="Set specific alignment score, format like ref:query=score,"
                             "Default=A:R=0,G:R=5,R:R=-4,A:N=0,G:N=0,C:N=0,T:N=0",
                        default="A:R=0,G:R=5,A:N=0,G:N=0,C:N=0,T:N=0")

    parser.add_argument("--dna_bulge_penalty",
                        help="A gap open and extend penalty on query sequence, default=12,8",
                        default="12,8", type=str)

    parser.add_argument("--rna_bulge_penalty",
                        help="A gap open and extend penalty on ref sequence, default=20,8",
                        default="20,8", type=str)

    parser.add_argument("--PAM_type_penalty",
                        help="PAM type penalty score for NGG, NAG and not NRG, default=0,4,12",
                        default="0,4,12")

    parser.add_argument("--dna_bulge_cmp_weight",
                        help="A weighted score for alignment comparison, "
                             "score=index * (start * coef - start) / query.len，"
                             "default=coef,start=1,12",
                        default="1,12")

    parser.add_argument("--rna_bulge_cmp_weight",
                        help="A weighted score for alignment comparison, "
                             "score=index * (start * coef - start) / query.len，"
                             "default=coef,start=1,16",
                        default="1,16")

    parser.add_argument("--mismatch_cmp_weight",
                        help="A weighted score for alignment comparison, "
                             "score=index * (start * coef - start) / query.len，"
                             "default=coef,start=2,5",
                        default="2,5")

    # parser.add_argument("--dist_to_signal_ref_point",
    #                     help="Dist to signal ref point, can be set as "
    #                          "'upstream', 'downstream', 'highest_site' or 'region_mid', default=upstream",
    #                     default="upstream")

    parser.add_argument("--dist_to_signal_penalty_breaks",
                        help="Dist to signal region penalty breaks, default=-10,-5,0,8,15",
                        default="-10,-5,0,8,15")

    parser.add_argument("--dist_to_signal_penalty_k",
                        help="Dist to signal region penalty k, default=-4,-2,0,0,2,4",
                        default="-4,-2,0,0,2,4")

    parser.add_argument("--dist_to_signal_penalty_offset",
                        help="Dist to signal region penalty offset, dist_to_signal_penalty=index * k + offset, "
                             "default=10,0,0,0,0,10",
                        default="10,0,0,0,0,10")

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
                                                                         highest_method=ARGS.find_highest_signal_method,
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
    # make step aligner
    # ------------------------------------------------------------------------>>>>>>>
    dna_bulge_gap_open, dna_bulge_gap_extension = map(int, str(ARGS.dna_bulge_penalty).split(","))
    rna_bulge_gap_open, rna_bulge_gap_extension = map(int, str(ARGS.rna_bulge_penalty).split(","))

    dna_bulge_gap_open = -1 * dna_bulge_gap_open
    dna_bulge_gap_extension = -1 * dna_bulge_gap_extension

    rna_bulge_gap_open = -1 * rna_bulge_gap_open
    rna_bulge_gap_extension = -1 * rna_bulge_gap_extension

    step_aligner = Bio.Align.PairwiseAligner()
    step_aligner.substitution_matrix = align_score_matrix
    step_aligner.open_gap_score = align_gap_open
    step_aligner.extend_gap_score = align_gap_extension
    step_aligner.query_left_gap_score = 0
    step_aligner.query_right_gap_score = 0
    step_aligner.mode = "global"

    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("Step aligner:" + "\n")
    sys.stderr.write(str(step_aligner))

    # reference Table 6.1: Meta-attributes of the pairwise aligner objects.
    # http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec118
    extend_aligner = Bio.Align.PairwiseAligner()
    extend_aligner.substitution_matrix = align_score_matrix
    extend_aligner.open_gap_score = align_gap_open
    extend_aligner.extend_gap_score = align_gap_extension
    extend_aligner.query_left_gap_score = 0
    extend_aligner.query_right_gap_score = 0
    extend_aligner.query_internal_open_gap_score = rna_bulge_gap_open
    extend_aligner.query_internal_extend_gap_score = rna_bulge_gap_extension
    extend_aligner.target_open_gap_score = dna_bulge_gap_open
    extend_aligner.target_extend_gap_score = dna_bulge_gap_extension
    extend_aligner.mode = "global"

    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("Extend aligner:" + "\n")
    sys.stderr.write(str(extend_aligner))
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
        "PAM_type",
        "PAM_seq",
        "align_total_match",
        "align_total_mismatch",
        "align_total_gap",
        "dna_bulge_num",
        "rna_bulge_num",
        "head_mismatch",
        "head_gap",
        "seed_mismatch",
        "seed_gap",
        "align_penalty_score",
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
    # load alignment penalty args
    # ------------------------------------------------------------------------>>>>>>>
    seed_start, seed_end = map(int, str(ARGS.seed_index).split(","))
    PAM_start, PAM_end = map(int, str(ARGS.PAM_index).split(","))
    head_start, head_end = map(int, str(ARGS.distal_index).split(","))
    pam_type_penalty_score = tuple(map(int, str(ARGS.PAM_type_penalty).split(",")))

    # load DNA and RNA bulge weight score
    dna_bulge_weight_coef, dna_bulge_start_score = map(int, str(ARGS.dna_bulge_cmp_weight).split(","))
    rna_bulge_weight_coef, rna_bulge_start_score = map(int, str(ARGS.rna_bulge_cmp_weight).split(","))
    mismatch_weight_coef, mismatch_start_score = map(int, str(ARGS.mismatch_cmp_weight).split(","))
    gap_weight_coef = dna_bulge_weight_coef
    gap_start_score = dna_bulge_start_score

    dist_signal_breaks_list = tuple(map(int, str(ARGS.dist_to_signal_penalty_breaks).split(",")))
    dist_signal_k_list = tuple(map(int, str(ARGS.dist_to_signal_penalty_k).split(",")))
    dist_signal_offset_list = tuple(map(int, str(ARGS.dist_to_signal_penalty_offset).split(",")))

    # ------------------------------------------------------------------------>>>>>>>
    # make query sequence
    # ------------------------------------------------------------------------>>>>>>>
    query_seq = Bio.Seq.Seq(ARGS.query_seq)
    seed_seq = query_seq[seed_start - 1:]
    extend_seq = query_seq[:seed_start - 1]

    # ------------------------------------------------------------------------>>>>>>>
    # alignment
    # ------------------------------------------------------------------------>>>>>>>
    if ARGS.input_header != "False":
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

        # load fasta seq
        fasta_seq_fwd = fa_info.seq.upper()
        fasta_seq_rev = Bio.Seq.reverse_complement(fasta_seq_fwd)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # step alignment
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        step_window_size = ARGS.align_step_window_size
        step_move_size = ARGS.align_step_move_size

        step_align_res_fwd = run_step_align(align_ref_seq=fasta_seq_fwd,
                                            align_query_seq=seed_seq,
                                            aligner=step_aligner,
                                            align_window_size=step_window_size,
                                            align_step_size=step_move_size,
                                            gap_weight_coef=gap_weight_coef,
                                            gap_start_score=gap_start_score,
                                            mis_weight_coef=mismatch_weight_coef,
                                            mis_start_score=mismatch_start_score)

        step_align_res_rev = run_step_align(align_ref_seq=fasta_seq_rev,
                                            align_query_seq=seed_seq,
                                            aligner=step_aligner,
                                            align_window_size=step_window_size,
                                            align_step_size=step_move_size,
                                            gap_weight_coef=gap_weight_coef,
                                            gap_start_score=gap_start_score,
                                            mis_weight_coef=mismatch_weight_coef,
                                            mis_start_score=mismatch_start_score)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # extend alignment
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        extend_align_res_fwd = run_extend_align(step_align_res_list=step_align_res_fwd,
                                                align_ref_seq=fasta_seq_fwd,
                                                seed_query_seq=seed_seq,
                                                extend_query_seq=extend_seq,
                                                run_extend_aligner=extend_aligner,
                                                extend_direction="left",
                                                gap_weight_coef=gap_weight_coef,
                                                gap_start_score=gap_start_score,
                                                mis_weight_coef=mismatch_weight_coef,
                                                mis_start_score=mismatch_start_score)

        extend_align_res_rev = run_extend_align(step_align_res_list=step_align_res_rev,
                                                align_ref_seq=fasta_seq_rev,
                                                seed_query_seq=seed_seq,
                                                extend_query_seq=extend_seq,
                                                run_extend_aligner=extend_aligner,
                                                extend_direction="left",
                                                gap_weight_coef=gap_weight_coef,
                                                gap_start_score=gap_start_score,
                                                mis_weight_coef=mismatch_weight_coef,
                                                mis_start_score=mismatch_start_score)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # analysis alignment results
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        align_analysis_list = []

        # fwd alignment
        for align_res in extend_align_res_fwd:
            # create analysis dict
            align_dict = weight_extend_align_obj_analysis(extend_align_res=align_res,
                                                          align_query_seq=query_seq,
                                                          pam_pos_idx=(PAM_start, PAM_end),
                                                          dna_bulge_coef=dna_bulge_weight_coef,
                                                          dna_bulge_start_score=dna_bulge_start_score,
                                                          rna_bulge_coef=rna_bulge_weight_coef,
                                                          rna_bulge_start_score=rna_bulge_start_score,
                                                          mis_coef=mismatch_weight_coef,
                                                          mis_start_score=mismatch_start_score,
                                                          pam_type_penalty=pam_type_penalty_score)
            # get genome coordinate
            align_chr_start, align_chr_end = get_align_chrom_coordinate(align_res,
                                                                        ref_chr_start=bed_chr_start,
                                                                        ref_chr_end=bed_chr_end,
                                                                        align_ref_start=align_res[3],
                                                                        align_ref_end=align_res[4],
                                                                        reverse_state=False)

            # calculate signal to dist penalty
            align_dist_to_signal = get_dist_to_signal(signal_chr_start=region_chr_start,
                                                      signal_chr_end=region_chr_end,
                                                      align_chr_start=align_chr_start,
                                                      align_chr_end=align_chr_end,
                                                      align_reverse=False)

            dist_penalty_score = signal_dist_penalty(dist_to_signal=align_dist_to_signal,
                                                     region_break_list=dist_signal_breaks_list,
                                                     region_k_list=dist_signal_k_list,
                                                     region_offset_list=dist_signal_offset_list)

            # calculate total penalty
            align_dict["align_strand"] = "+"
            align_dict["align_chr_start"] = align_chr_start
            align_dict["align_chr_end"] = align_chr_end

            align_dict["region_chr_name"] = bed_list[0]
            align_dict["region_chr_start"] = region_chr_start
            align_dict["region_chr_end"] = region_chr_end
            align_dict["region_index"] = bed_list[3]

            align_dict["dist_penalty"] = dist_penalty_score
            align_dict["dist_align_to_signal"] = align_dist_to_signal

            align_dict["total_penalty"] = align_dict["gap_penalty_score"] + \
                                          align_dict["mismatch_penalty_score"] + \
                                          align_dict["dist_penalty"] + \
                                          align_dict["PAM_type_penalty_score"]

            # add into list
            align_analysis_list.append(align_dict)

        # rev alignment
        for align_res in extend_align_res_rev:
            # create analysis dict
            align_dict = weight_extend_align_obj_analysis(extend_align_res=align_res,
                                                          align_query_seq=query_seq,
                                                          pam_pos_idx=(PAM_start, PAM_end),
                                                          dna_bulge_coef=dna_bulge_weight_coef,
                                                          dna_bulge_start_score=dna_bulge_start_score,
                                                          rna_bulge_coef=rna_bulge_weight_coef,
                                                          rna_bulge_start_score=rna_bulge_start_score,
                                                          mis_coef=mismatch_weight_coef,
                                                          mis_start_score=mismatch_start_score,
                                                          pam_type_penalty=pam_type_penalty_score)
            # get genome coordinate
            align_chr_start, align_chr_end = get_align_chrom_coordinate(align_res,
                                                                        ref_chr_start=bed_chr_start,
                                                                        ref_chr_end=bed_chr_end,
                                                                        align_ref_start=align_res[3],
                                                                        align_ref_end=align_res[4],
                                                                        reverse_state=True)

            # calculate signal to dist penalty
            align_dist_to_signal = get_dist_to_signal(signal_chr_start=region_chr_start,
                                                      signal_chr_end=region_chr_end,
                                                      align_chr_start=align_chr_start,
                                                      align_chr_end=align_chr_end,
                                                      align_reverse=True)

            dist_penalty_score = signal_dist_penalty(dist_to_signal=align_dist_to_signal,
                                                     region_break_list=dist_signal_breaks_list,
                                                     region_k_list=dist_signal_k_list,
                                                     region_offset_list=dist_signal_offset_list)

            # calculate total penalty
            align_dict["align_strand"] = "-"
            align_dict["align_chr_start"] = align_chr_start
            align_dict["align_chr_end"] = align_chr_end

            align_dict["region_chr_name"] = bed_list[0]
            align_dict["region_chr_start"] = region_chr_start
            align_dict["region_chr_end"] = region_chr_end
            align_dict["region_index"] = bed_list[3]

            align_dict["dist_penalty"] = dist_penalty_score
            align_dict["dist_align_to_signal"] = align_dist_to_signal

            align_dict["total_penalty"] = align_dict["gap_penalty_score"] + \
                                          align_dict["mismatch_penalty_score"] + \
                                          align_dict["dist_penalty"] + \
                                          align_dict["PAM_type_penalty_score"]

            # add into list
            align_analysis_list.append(align_dict)

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # select the best alignment
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        align_analysis_list_sort = sorted(align_analysis_list,
                                          key=lambda i: (i["total_penalty"],
                                                         i["gap_penalty_score"],
                                                         i["mismatch_penalty_score"],
                                                         i["dist_penalty"]))

        # select the best alignment
        best_align_dict = align_analysis_list_sort[0]

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # output best alignment
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        out_list = report_align_res(input_align_dict=best_align_dict,
                                    pam_proximal_idx=(seed_start, seed_end),
                                    pam_distal_idx=(head_start, head_end))

        out_list_str = "\t".join(map(str, out_list))
        out_file.write(out_list_str + "\n")

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


# 2021-12-07 by Meng Haowei
