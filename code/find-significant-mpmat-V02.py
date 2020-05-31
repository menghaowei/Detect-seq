#! /Users/meng/menghw_HD/anaconda2/bin/python
# #! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
    2019-12-30
        Input ctrl json and treat json respectiv

    2019-12-29
        Fix output part 

    2019-12-28 
        Find significant mpmat signal

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""

Design pipeline:

Input:
    1. ctrl_BAM treat_BAM
    2. efficient genome size 
    3. mpmat region 
    4. reference
    5. mut_type_list
    6. mut_num_cutoff

Output:
    1. scale factor on each chromosome
    2. region raw mutant reads 
    3. region raw non-mutant reads 
    4. region fix mutant reads 
    5. region fix non-mutant reads
    6. log2FoldChange 
    7. pvalue 
    8. qvalue 

"""
# Learning Part END-----------------------------------------------------------


import argparse
import os
import sys

import random
import string
import time
import json

import pysam 
from Bio import  SeqIO

from scipy import stats
import statsmodels.stats.multitest as multi
import numpy as np
import math

from DetectSeqLib.FindMutRegion import *

###############################################################################
# function part 
###############################################################################
def _count_dict_sum_up(count_dict, key_ref_name_list="All", ignore_key_list = ["total"]):
    """
    INPUT 
        <count_dict>
            key = chr1, chr2, chr3 ...
            value = dict and key is 
            
                  non_mut_align_count
                  all_align_count
                  total_high_mismatch_count
                  query_high_mismatch_count
                  mut_align_count
                  other_high_mismatch_coun
                  all_filter_count
                  
        <key_ref_name_list>
            use those key to sum up and count 
            
        <ignal_key_list>
            ignore  
    
    RETURN 
        Dict like <count_dict> and add 'total' key
    """
    
    if key_ref_name_list == "All":
        key_ref_name_list = count_dict.keys()
    
    total_dict = {}
    inner_key_list = count_dict[key_ref_name_list[0]].keys()

    for inner_key in inner_key_list:
        total_dict[inner_key] = 0

    for ref_key in key_ref_name_list:
        if ref_key not in ignore_key_list:
            for ref_inner_key in inner_key_list:
                total_dict[ref_inner_key] += count_dict[ref_key][ref_inner_key]
    
    count_dict["total"] = total_dict.copy()
    
    return(count_dict)


def _calculate_scale_dict(meta_data_dict, to_large_state = False, scale_reads_count=None):
    """
    INPUT:
        <meta_data_dict>
            A dict generated from raw BAM file, and data structure like
                meta_count_dict = {
                    "meta_data" : {
                        "create_date" : ,
                        "ctrl_BAM" : ,
                        "treat_BAM" : ,
                        "reference" : ,
                        "threads" : ,
                        "query_mutation_type" : ,
                        "query_mutation_min_cutoff" : ,
                        "query_mutation_max_cutoff" : ,
                        "other_mismatch_max_cutoff" : ,
                        "total_mismatch_max_cutoff" : ,
                        "filter_log" : 
                    },
                    "ctrl" : ,
                    "treat" : ,
                }
            
        <to_large_state>
            Default is False
                Means scale the larger sample up to the smaller sample. 
                Usually, scaling down will bring down background noise, which is good for enrichment detection.
    
    RETURN:
        <scale_dict>
        
        
    INFO:
        Final changes, 2019-12-29
    """
    
    # set dict
    scale_dict = {"ctrl" : {}, "treat" : {}}
    
    select_key_list = [
        "all_align_count",
        "mut_align_count",
        "non_mut_align_count"
    ]
    
    if  scale_reads_count != None:
        ctrl_total_scale = (meta_data_dict["ctrl"]["total"]["all_align_count"]) / 1.0 /  int(scale_reads_count) 
        treat_total_scale = (meta_data_dict["treat"]["total"]["all_align_count"]) / 1.0 / int(scale_reads_count) 
    
    for ref_name in meta_data_dict["ctrl"].keys():
        # make dict 
        if scale_dict["ctrl"].get(ref_name) == None:
            scale_dict["ctrl"][ref_name] = {}

        if scale_dict["treat"].get(ref_name) == None:
            scale_dict["treat"][ref_name] = {}
            
        for key in meta_data_dict["ctrl"][ref_name].keys():
            ctrl_count = meta_data_dict["ctrl"][ref_name][key]
            treat_count = meta_data_dict["treat"][ref_name][key]
            
            if scale_reads_count == None:
                if to_large_state:
                    scale_value = max(ctrl_count, treat_count)
                else:
                    scale_value = min(ctrl_count, treat_count)

                if scale_value == 0:
                    ctrl_scale_factor = None
                    treat_scale_factor = None
                else:
                    ctrl_scale_factor = ctrl_count / 1.0 / scale_value
                    treat_scale_factor = treat_count / 1.0 / scale_value
            else:
                ctrl_scale_factor = ctrl_total_scale
                treat_scale_factor = treat_total_scale
            
            # make new dict key
            if key in select_key_list:
                scale_dict_key = key + ".scale_factor"
                scale_dict["ctrl"][ref_name][scale_dict_key] = ctrl_scale_factor
                scale_dict["treat"][ref_name][scale_dict_key] = treat_scale_factor

    # return part 
    return scale_dict
        

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="A tool to find significant enrichment mpmat region")

    parser.add_argument("-i", "--mpmat_table",
        help=".mpmat table, which can be generated from <pmat-merge.py> code",required=True)

    parser.add_argument("-c", "--ctrl_BAM",
        help="Control BAM file",required=True)

    parser.add_argument("-t", "--treat_BAM",
        help="Treatment BAM file",required=True)

    parser.add_argument("-r", "--reference",
        help="Genome FASTA file",required=True)

    parser.add_argument("-g", "--genome_effective_json",
        help="Genome effective sequence information in .json format, can be generated from script <count-effective-genome>",required=True)

    parser.add_argument("-m", "--ctrl_BAM_json",
        help="Ctrl BAM json, meta data information in .json format, generated from <calculate-mut-stats.py>",required=True)

    parser.add_argument("-n", "--treat_BAM_json",
        help="Treat BAM json, meta data information in .json format, generated from <calculate-mut-stats.py>",required=True)

    parser.add_argument("-o", "--output",
        help="Output file, default=stdout",default="stdout")

    parser.add_argument("--region_mutation_min_cutoff",
        help="No less than this cutoff consider as a region mutant read, default=1",default=1, type=int)    

    parser.add_argument("--query_mutation_type",
        help="Reads mutation type, which considered contain useful information, default=CT,GA",default="CT,GA")

    parser.add_argument("--query_mutation_min_cutoff",
        help="An alignment contain mutation number lower than this, considered as 'non_mut_align', or as 'mut_align', default=2",default=2, type=int)

    parser.add_argument("--query_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'query_high_mismatch', which often caused by Bismark mapping error, default=16",default=16, type=int)

    parser.add_argument("--other_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'other_high_mismatch', default=12",default=12, type=int)

    parser.add_argument("--total_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'total_high_mismatch', default=24",default=24, type=int)

    parser.add_argument("--to_large",
        help="Means scale the larger sample up to the smaller sample. Usually, \
        scaling down will bring down background noise, \
        which is good for enrichment detection. default=False",default=False, type=bool)

    parser.add_argument("--scale_reads_count",
        help="Scaled final output region signal, default=100e6 Usually, default=100e6",default=int(100e6), type=int)

    parser.add_argument("--mpmat_header",
        help="If set True, means mpmat table contain header line, default=False",default=False, type=bool)

    parser.add_argument("--lambda_method",
        help="Can be set as 'ctrl_max', 'treat_max', 'max', default=ctrl_max",default="ctrl_max")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    mpmat_filename = ARGS.mpmat_table
    ctrl_bam_filename = ARGS.ctrl_BAM
    treat_bam_filename = ARGS.treat_BAM
    
    query_mut_type_list = ARGS.query_mutation_type.split(",")
    region_mut_min_cutoff = ARGS.region_mutation_min_cutoff
    query_mut_min_cutoff = ARGS.query_mutation_min_cutoff
    query_mut_max_cutoff = ARGS.query_mutation_max_cutoff
    other_mut_max_cutoff = ARGS.other_mutation_max_cutoff
    total_mut_max_cutoff = ARGS.total_mutation_max_cutoff
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make full cmd
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    full_cmd_str = "python find-significant-mpmat.py \n \
    --mpmat_table {mpmat_table} \n \
    --ctrl_BAM {ctrl_BAM} \n \
    --treat_BAM {treat_BAM} \n \
    --reference {reference} \n \
    --genome_effective_json {genome_effective_json} \n \
    --ctrl_BAM_json {ctrl_BAM_json} \n \
    --treat_BAM_json {treat_BAM_json} \n \
    --output {output} \n \
    --query_mutation_type {query_mutation_type} \n \
    --region_mutation_min_cutoff {region_mutation_min_cutoff} \n \
    --query_mutation_min_cutoff {query_mutation_min_cutoff} \n \
    --query_mutation_max_cutoff {query_mutation_max_cutoff} \n \
    --other_mutation_max_cutoff {other_mutation_max_cutoff} \n \
    --total_mutation_max_cutoff {total_mutation_max_cutoff} \n \
    --scale_reads_count {scale_reads_count} \n \
    --lambda_method {lambda_method} \n \
    --mpmat_header {mpmat_header} \
    ".format(
        mpmat_table = ARGS.mpmat_table,
        ctrl_BAM = ARGS.ctrl_BAM,
        treat_BAM = ARGS.treat_BAM,
        reference = ARGS.reference,
        genome_effective_json = ARGS.genome_effective_json,
        ctrl_BAM_json = ARGS.ctrl_BAM_json,
        treat_BAM_json = ARGS.treat_BAM_json,
        output = ARGS.output,
        query_mutation_type = ARGS.query_mutation_type,
        region_mutation_min_cutoff = ARGS.region_mutation_min_cutoff,
        query_mutation_min_cutoff = ARGS.query_mutation_min_cutoff,
        query_mutation_max_cutoff = ARGS.query_mutation_max_cutoff,
        other_mutation_max_cutoff = ARGS.other_mutation_max_cutoff,
        total_mutation_max_cutoff = ARGS.total_mutation_max_cutoff,
        scale_reads_count = ARGS.scale_reads_count,
        lambda_method = ARGS.lambda_method,
        mpmat_header = ARGS.mpmat_header
    )

    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(full_cmd_str)
    sys.stderr.write("\n" + "-" * 80 + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check input bam file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    check_res = check_input_bam_file(ctrl_bam_filename)

    if check_res == 1:
        raise IOError("BAM not exisit!")

    elif check_res == 2:
        raise IOError("BAM not sorted by coordinate!")

    elif check_res == 3:
        raise IOError("BAM doesn't contain index file!")


    check_res = check_input_bam_file(treat_bam_filename)

    if check_res == 1:
        raise IOError("BAM not exisit!")

    elif check_res == 2:
        raise IOError("BAM not sorted by coordinate!")

    elif check_res == 3:
        raise IOError("BAM doesn't contain index file!")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check other file exist
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if not os.path.exists(ARGS.mpmat_table):
        raise IOError("--mpmat_table does not exist! ")

    if not os.path.exists(ARGS.reference):
        raise IOError("--reference does not exist! ")

    if not os.path.exists(ARGS.genome_effective_json):
        raise IOError("--genome_effective_json does not exist! ")

    if not os.path.exists(ARGS.ctrl_BAM_json):
        raise IOError("--ctrl_BAM_json does not exist! ")

    if not os.path.exists(ARGS.treat_BAM_json):
        raise IOError("--treat_BAM_json does not exist! ")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load json file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load meta data dict
    meta_data_dict = {}

    # load ctrl 
    ctrl_meta_dict = json.load(open(ARGS.ctrl_BAM_json, "rb"))
    meta_data_dict["ctrl"] = ctrl_meta_dict["count_dict"].copy()

    # load treat
    treat_meta_dict = json.load(open(ARGS.treat_BAM_json, "rb"))
    meta_data_dict["treat"] = treat_meta_dict["count_dict"].copy()    

    # test key 
    if treat_meta_dict["count_dict"].keys() != ctrl_meta_dict["count_dict"].keys():
        raise ValueError("Ctrl json and Treat json file not match!")

    # add total key to the dict 
    meta_data_dict["ctrl"] = _count_dict_sum_up(meta_data_dict["ctrl"])
    meta_data_dict["treat"] = _count_dict_sum_up(meta_data_dict["treat"])

    # load genome count json 
    genome_base_count_dict = None
    with open(ARGS.genome_effective_json,"r") as genome_info_file:
        genome_base_count_dict = json.load(genome_info_file)

    # make scale factor dict
    scale_factor_dict = _calculate_scale_dict(meta_data_dict, to_large_state=False).copy()

    # make normalize scale factor
    nomalize_scale_factor_dict = _calculate_scale_dict(meta_data_dict, to_large_state=False, scale_reads_count=ARGS.scale_reads_count).copy()

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load genome as dict
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Start to load genome ... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    genome_dict = load_reference_fasta_as_dict(ref_fasta_path=ARGS.reference, ref_name_list="All")

    sys.stderr.write("Start to load genome. Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # calculate genome and chromosome background
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    bg_dict = {
        "ctrl" : {}, 
        "treat" : {}
    }

    # genome background
    ctrl_bg_genome_lambda = meta_data_dict["ctrl"]["total"]["mut_align_count"] / 1.0 / scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]  /  (genome_base_count_dict["total"]["A"] + genome_base_count_dict["total"]["T"] + genome_base_count_dict["total"]["C"] + genome_base_count_dict["total"]["G"]) * 300
    treat_bg_genome_lambda = meta_data_dict["treat"]["total"]["mut_align_count"] / 1.0 / scale_factor_dict["treat"]["total"]["all_align_count.scale_factor"]  /  (genome_base_count_dict["total"]["A"] + genome_base_count_dict["total"]["T"] + genome_base_count_dict["total"]["C"] + genome_base_count_dict["total"]["G"]) * 300

    bg_dict["ctrl"]["genome_bg"] = ctrl_bg_genome_lambda
    bg_dict["treat"]["genome_bg"] = treat_bg_genome_lambda

    # chromosome backround
    for ref_name in meta_data_dict["ctrl"].keys():
        ctrl_bg_chr_lambda = meta_data_dict["ctrl"][ref_name]["mut_align_count"] / 1.0 /  scale_factor_dict["ctrl"][ref_name]["all_align_count.scale_factor"]   / (genome_base_count_dict[ref_name]["A"] + genome_base_count_dict[ref_name]["T"] + genome_base_count_dict[ref_name]["C"] + genome_base_count_dict[ref_name]["G"]) * 300
        treat_bg_chr_lambda = meta_data_dict["treat"][ref_name]["mut_align_count"] / 1.0 /  scale_factor_dict["treat"][ref_name]["all_align_count.scale_factor"]   / (genome_base_count_dict[ref_name]["A"] + genome_base_count_dict[ref_name]["T"] + genome_base_count_dict[ref_name]["C"] + genome_base_count_dict[ref_name]["G"]) * 300
        
        bg_dict["ctrl"][ref_name] = ctrl_bg_chr_lambda
        bg_dict["treat"][ref_name] = treat_bg_chr_lambda

    # print background info to stderr
    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("Out background info:  \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
    sys.stderr.write(json.dumps(bg_dict, encoding="utf-8", sort_keys=True, indent=4, separators=(', ', ': ')) + "\n")
    sys.stderr.write("-" * 80 + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load mpmat and calculate 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ctrl_bam_obj = pysam.AlignmentFile(ctrl_bam_filename, "rb")
    treat_bam_obj = pysam.AlignmentFile(treat_bam_filename, "rb")
    
    mpmat_file = open(mpmat_filename, "rb")
    mpmat_header = None
    mpmat_header_list = []

    if ARGS.mpmat_header:
        mpmat_header = mpmat_file.readline()
        mpmat_header_list = mpmat_header.strip().split("\t")

    # store pvalue in each region
    poisson_pvalue_list = []

    mpmat_info_dict = {
        "ctrl_all_count" : [],
        "treat_all_count" : [],
        "ctrl_mut_count" : [],
        "treat_mut_count" : [],
        "ctrl_all_count.norm" : [],
        "treat_all_count.norm" : [],
        "ctrl_mut_count.norm" : [],
        "treat_mut_count.norm" : [],
        "log2FC_all_count" : [],
        "log2FC_mut_count" : []
    }

    for line_idx, line in enumerate(mpmat_file):
        if line_idx % 100 == 0:
            sys.stderr.write("Processing line number: %s \t %s \n" % (line_idx + 1, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        line_list = line.strip().split("\t")

        try:
            mpmat_chr_name = line_list[0]
            mpmat_start = int(line_list[1])
            mpmat_end = int(line_list[2])
        except:
            raise ValueError("Input mpmat file have to contain chr_name, region_start, region_end in the first three columns, and separated by \\t")

        ctrl_mpmat_count_dict = \
            get_mpmat_region_count(
                in_bam_obj = ctrl_bam_obj,
                mpmat_chr_name = mpmat_chr_name,
                mpmat_region_start = mpmat_start,
                mpmat_region_end = mpmat_end,
                ref_genome_dict = genome_dict,
                query_mismatch_type_list = query_mut_type_list, 
                other_mismatch_type_list = "All",
                region_mismatch_type_list = query_mut_type_list,
                region_query_mut_min_cutoff = region_mut_min_cutoff,
                query_mut_max_cutoff = query_mut_max_cutoff,
                total_mut_max_cutoff = total_mut_max_cutoff,
                other_mut_max_cutoff = other_mut_max_cutoff
            )

        treat_mpmat_count_dict = \
            get_mpmat_region_count(
                in_bam_obj = treat_bam_obj,
                mpmat_chr_name = mpmat_chr_name,
                mpmat_region_start = mpmat_start,
                mpmat_region_end = mpmat_end,
                ref_genome_dict = genome_dict,
                query_mismatch_type_list = query_mut_type_list, 
                other_mismatch_type_list = "All",
                region_mismatch_type_list = query_mut_type_list,
                region_query_mut_min_cutoff = region_mut_min_cutoff,
                query_mut_max_cutoff = query_mut_max_cutoff,
                total_mut_max_cutoff = total_mut_max_cutoff,
                other_mut_max_cutoff = other_mut_max_cutoff
            )

        # region lambda
        ctrl_mpmat_lambda = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 /  scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"] 
        treat_mpmat_lambda = treat_mpmat_count_dict["region_mut_count"] / 1.0 /  scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"] 

        # all reads lambda 
        ctrl_mpmat_lambda_all = ctrl_mpmat_count_dict["all_align_count"] / 1.0 /  scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"] 
        treat_mpmat_lambda_all = treat_mpmat_count_dict["all_align_count"] / 1.0 /  scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"] 

        # make lambda list 
        ctrl_bg_lambda_list = [
            bg_dict["ctrl"]["genome_bg"],
            bg_dict["ctrl"][mpmat_chr_name],
            ctrl_mpmat_lambda    
        ]

        treat_bg_lambda_list = [
            bg_dict["treat"]["genome_bg"],
            bg_dict["treat"][mpmat_chr_name]
        ]

        bg_lambda_list = [
            bg_dict["ctrl"]["genome_bg"],
            bg_dict["ctrl"][mpmat_chr_name],
            ctrl_mpmat_lambda,
            bg_dict["treat"]["genome_bg"],
            bg_dict["treat"][mpmat_chr_name]
        ]

        # select bg lambda
        if ARGS.lambda_method == "ctrl_max":
            bg_lambda = max(ctrl_bg_lambda_list)

        elif ARGS.lambda_method == "treat_max":
            bg_lambda = max(treat_bg_lambda_list)

        elif ARGS.lambda_method == "max":
            bg_lambda = max(bg_lambda_list)

        # mut reads get pvalue 
        mut_pvalue = poisson_test(treat_mpmat_lambda, bg_lambda, alternative="larger", method="score")[1]
        poisson_pvalue_list.append(mut_pvalue)

        # all reads get pvalue 
        # all_pvalue = poisson_test(treat_mpmat_lambda_all, ctrl_mpmat_lambda_all, alternative="larger", method="score")[1]

        # calculate normalized signal
        mpmat_ctrl_mut_count_nrom = ctrl_mpmat_count_dict["all_align_count"] / 1.0 / nomalize_scale_factor_dict["ctrl"][mpmat_chr_name]["all_align_count.scale_factor"]
        mpmat_ctrl_all_count_nrom = ctrl_mpmat_count_dict["region_mut_count"] / 1.0 / nomalize_scale_factor_dict["ctrl"][mpmat_chr_name]["mut_align_count.scale_factor"]

        mpmat_treat_mut_count_nrom = treat_mpmat_count_dict["all_align_count"] / 1.0 / nomalize_scale_factor_dict["treat"][mpmat_chr_name]["all_align_count.scale_factor"]
        mpmat_treat_all_count_nrom = treat_mpmat_count_dict["region_mut_count"] / 1.0 / nomalize_scale_factor_dict["treat"][mpmat_chr_name]["mut_align_count.scale_factor"]

        # calculate log2 fold change
        # if mpmat_ctrl_all_count_nrom != 0:
        #     log2FC_all_count = math.log(mpmat_treat_all_count_nrom / 1.0 / mpmat_ctrl_all_count_nrom, 2)
        # else:
        #     log2FC_all_count = math.log(mpmat_treat_all_count_nrom / 1.0, 2)

        # if mpmat_ctrl_mut_count_nrom != 0:
        #     log2FC_mut_count = math.log(mpmat_treat_mut_count_nrom / 1.0 / mpmat_ctrl_mut_count_nrom, 2)
        # else:
        #     log2FC_mut_count = math.log(mpmat_treat_mut_count_nrom / 1.0, 2)

        # All -> calculate log2 fold change fix 
        if mpmat_ctrl_all_count_nrom != 0:
            all_count_FC = mpmat_treat_all_count_nrom / 1.0 / mpmat_ctrl_all_count_nrom
        else:
            all_count_FC = mpmat_treat_all_count_nrom / 1.0

        if all_count_FC > 0:
            log2FC_all_count = math.log(all_count_FC, 2)
        else:
            log2FC_all_count = "NA"

        # Mut -> calculate log2 fold change fix 
        if mpmat_ctrl_mut_count_nrom != 0:
            mut_count_FC = mpmat_treat_mut_count_nrom / 1.0 / mpmat_ctrl_mut_count_nrom
        else:
            mut_count_FC = mpmat_treat_mut_count_nrom / 1.0

        if mut_count_FC > 0:
            log2FC_mut_count = math.log(mut_count_FC, 2)
        else:
            log2FC_mut_count = "NA"

        # add all info into dict
        mpmat_info_dict["ctrl_all_count"].append(ctrl_mpmat_count_dict["all_align_count"])
        mpmat_info_dict["treat_all_count"].append(treat_mpmat_count_dict["all_align_count"])

        mpmat_info_dict["ctrl_mut_count"].append(ctrl_mpmat_count_dict["region_mut_count"])
        mpmat_info_dict["treat_mut_count"].append(treat_mpmat_count_dict["region_mut_count"])

        mpmat_info_dict["ctrl_all_count.norm"].append(mpmat_ctrl_all_count_nrom)
        mpmat_info_dict["treat_all_count.norm"].append(mpmat_treat_all_count_nrom)

        mpmat_info_dict["ctrl_mut_count.norm"].append(mpmat_ctrl_mut_count_nrom)
        mpmat_info_dict["treat_mut_count.norm"].append(mpmat_treat_mut_count_nrom)

        mpmat_info_dict["log2FC_all_count"].append(log2FC_all_count)
        mpmat_info_dict["log2FC_mut_count"].append(log2FC_mut_count)

        # add pvalue into list
        # poisson_pvalue_list.append(max(mut_pvalue, all_pvalue))

        # add lambda into list 
        # treat_lambda_list.append(treat_mpmat_lambda)
        # ctrl_lambda_list.append(ctrl_mpmat_lambda)

        # treat_all_lambda_list.append(treat_mpmat_lambda_all)
        # ctrl_all_lambda_list.append(ctrl_mpmat_lambda_all)

        # run_bg_lambda_list.append(bg_lambda)

    # close file
    mpmat_file.close()
    ctrl_bam_obj.close()
    treat_bam_obj.close()

    # End for loop and fix pvalue to FDR
    FDR_qvalue = multi.multipletests(np.array(poisson_pvalue_list), alpha=0.05, method="fdr_bh", is_sorted=False)
    FDR_qvalue_vec = FDR_qvalue[1]

    # out final mpmat 
    if ARGS.output == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(ARGS.output, "wb")

    out_header_list = [
        "chr_name",
        "region_start",
        "region_end",
        "ctrl_count",
        "treat_count",
        "ctrl_mut_count",
        "treat_mut_count",
        "ctrl_count.norm",
        "treat_count.norm",
        "ctrl_mut_count.norm",
        "treat_mut_count.norm",
        "log2_FC",
        "log2_FC_mut",
        "p_value",
        "FDR"
    ]

    with open(mpmat_filename, "rb") as mpmat_raw_file:
        # out header 
        mpmat_header_str = "\t".join(out_header_list)
        out_file.write(mpmat_header_str + "\n")

        if ARGS.mpmat_header:
            header = mpmat_raw_file.readline()

        for index, line in enumerate(mpmat_raw_file):
            line_list = line.strip().split("\t")[:3]
            
            append_info_list = [
                mpmat_info_dict["ctrl_all_count"][index],
                mpmat_info_dict["treat_all_count"][index],
                mpmat_info_dict["ctrl_mut_count"][index],
                mpmat_info_dict["treat_mut_count"][index],
                mpmat_info_dict["ctrl_all_count.norm"][index],
                mpmat_info_dict["treat_all_count.norm"][index],
                mpmat_info_dict["ctrl_mut_count.norm"][index],
                mpmat_info_dict["treat_mut_count.norm"][index],
                mpmat_info_dict["log2FC_all_count"][index],
                mpmat_info_dict["log2FC_mut_count"][index],
                poisson_pvalue_list[index],
                FDR_qvalue_vec[index]
            ]

            line_list += append_info_list
            line_str = "\t".join(map(str,line_list))
            out_file.write(line_str + "\n")

    # close output file 
    if ARGS.output != "stdout":
        out_file.close()








            


# 2019-12-28

