#! /Users/meng/menghw_HD/anaconda2/bin/python
#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-02:
    2020-01-03
        Fix NM tag bug

Version-01:
    2019-12-30 
        Only keep stats function
    
    2019-12-23 
        Find significant mpmat signal

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""

Design pipeline:

Input:
    1. Input BAM
    2. reference
    3. mut_type_list
    4. mut_num_cutoff

Output:
    1. BAM json

"""
# Learning Part END-----------------------------------------------------------


import argparse
import os
import random
import string
import sys
import time
import multiprocessing
import json

import pysam 
from Bio import  SeqIO

from DetectSeqLib.FindMutRegion import *
###############################################################################
# function part 
###############################################################################
def run_BAM_count_on_chromosme(
    in_bam_filename, 
    select_chr_name, 
    ref_genome_path,
    query_mismatch_type_list = ["CT","GA"], 
    other_mismatch_type_list = "All",
    query_mut_min_cutoff = 2,
    query_mut_max_cutoff = 16,
    total_mut_max_cutoff = 16,
    other_mut_max_cutoff = 16,
    filter_align_log_state = True,
    filter_align_log_filename = "Stderr"):
    """
    INPUT:
        <in_bam_filename>

        <select_chr_name>

        <query_mismatch_type_list>
            like ["CT","GA"], and set "All" means count every type of mismatch mutations

        <other_mismatch_type_list>
            like ["GT"], and set "All" means count every type of mismatch mutations, which not in <query_mismatch_type_list>
            
        <query_mut_min_cutoff>
            lower than this, considered as "non_mut_align"
            no less than this, considered as "mut_align"
        
        <query_mut_max_cutoff>
            higher than this, considered as "query_high_mismatch"

        <total_mut_max_cutoff>
            higher than this, considered as "total_high_mismatch"
        
        <other_mut_max_cutoff>
            higher than this, considered as "other_high_mismatch"
        
            !!! The priority is  total_mut_max_cutoff > other_mut_max_cutoff > query_mut_max_cutoff > query_mut_min_cutoff
        
        <filter_align_log_state>
            True meas output align filter region
        
        <filter_align_log_filename>        
    
    RETURN:
        if everything word well return:
            dict = {
                "all_align_count" : val,
                "mut_align_count" : val,
                "non_mut_align_count" : val,
                "all_filter_count" : val,
                "total_high_mismatch_count" : val,
                "other_high_mismatch_count" : val,
                "query_high_mismatch_count" : val
            }
     
        something wrong return:
            None
    """
    # load fasta file
    ref_genome_dict = load_reference_fasta_as_dict(ref_genome_path,ref_name_list=[select_chr_name],show_load_id=True)

    # open log file 
    if filter_align_log_state:
        if filter_align_log_filename == "Stderr":
            filter_log_file = sys.stderr
        else:
            filter_log_file = open(filter_align_log_filename, "wb")
    
    # define dict 
    align_count_dict = {
        "all_align_count" : 0,
        "mut_align_count" : 0,
        "non_mut_align_count" : 0,
        "all_filter_count" : 0,
        "total_high_mismatch_count" : 0,
        "other_high_mismatch_count" : 0,
        "query_high_mismatch_count" : 0
    }
    
    # open BAM file
    in_bam_file = pysam.AlignmentFile(in_bam_filename, "rb")
    
    # run
    for align in in_bam_file.fetch(reference = select_chr_name):
        # count total
        align_count_dict["all_align_count"] +=1

        # make NM state 
        NM_state = False
        NM_tag_str = None
        
        try:
            NM_tag_str = align.get_tag("NM")
            NM_state = True
        except:
            NM_state = False
            align_count_dict["all_filter_count"] += 1
            continue

        # make sure MD state 
        MD_state = False
        MD_tag_str = None

        try:
            MD_tag_str = align.get_tag("MD")
            MD_state = True
        except:
            MD_state = False

        # init
        align_mismatch_pairs = None

        # get align mismatch info
        if MD_state:
            #  calculate mismatch info with MD tag 
            align_mismatch_pairs = get_align_mismatch_pairs(align)

        if (not MD_state) or (align_mismatch_pairs == None):
            # calculate mismatch info without MD tag 
            align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)

        # count mismatch info
        if align_mismatch_pairs != None:
            count_dict = make_mismatch_count_dict(align_mismatch_pairs)
            mut_count_res = count_mismatch_dict(count_dict, query_key_list=["CT","GA"], other_key_list="All")
            query_mut_count, other_mut_count, total_mut_count, region_mut_count = mut_count_res

            filter_reason = ""

            # filter part
            if total_mut_count >= total_mut_max_cutoff:
                align_count_dict["total_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Total High Mismatch Count"

            elif other_mut_count >= other_mut_max_cutoff:
                align_count_dict["other_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Other High Mismatch Count"

            elif query_mut_count >= query_mut_max_cutoff:
                align_count_dict["query_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Query High Mismatch Count"

            else:
                if query_mut_count < query_mut_min_cutoff:
                    # non mut 
                    align_count_dict["non_mut_align_count"] +=1
                else:
                    # mut count 
                    align_count_dict["mut_align_count"] +=1

            # filter log 
            if filter_align_log_state:
                if filter_reason != "":
                    filter_log_list = ["Filter Align", filter_reason] + mut_count_res + [align.query_name, align.reference_name, align.reference_start]
                    filter_log_str = "\t".join(map(str,filter_log_list))
                    filter_log_file.write(filter_log_str + "\n")

        else:
            # no mismatch align
            align_count_dict["non_mut_align_count"] +=1
    
    # close file 
    in_bam_file.close()
    
    if filter_align_log_filename != "Stderr":
        filter_log_file.close()
    
    # return part 
    return([select_chr_name, align_count_dict])


def merge_log_file(in_filename_list, out_filename = "stdout", out_mode="wb", rm_input_file = True):
    """
    INPUT
        <filename_list>
        
        <out_filename>
            output filename, default = stdout 
        
        <out_mode>
            set open file mode 
        
        <rm_input_file>
            default is True, means remove all in_filename_list files
    
    RETURN
        0, means everthing is okay !
    """
    
    # set output 
    if out_filename == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(out_filename, out_mode)
    
    # merge 
    for in_filename in in_filename_list:
        try:
            in_file = open(in_filename, "rb")
        except:
            continue
        
        for line in in_file:
            out_file.write(line)
        
        in_file.close()

        # remove 
        if rm_input_file:
            os.remove(in_filename)
    
    # close file 
    if out_filename != "stdout":
        out_file.close()
    
    return(0)
    

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Filter BAM high mismatches alignment")

    parser.add_argument("-i", "--input_bam",
        help="Input BAM file have to sorted and indexed",required=True)

    parser.add_argument("-r", "--reference",
        help="Genome FASTA file",required=True)

    parser.add_argument("-p", "--threads",
        help="Multiple threads number, default=1",default=1, type=int)

    parser.add_argument("-o","--out_json",
        help="Output a file contain count results as .json format, default=stdout",default="stdout")

    parser.add_argument("--query_mutation_type",
        help="Reads mutation type, which considered contain useful information, default=CT,GA",default="CT,GA")

    parser.add_argument("--query_mutation_min_cutoff",
        help="An alignment contain mutation number lower than this, considered as 'non_mut_align', or as 'mut_align', default=2",default=2, type=int)

    parser.add_argument("--query_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'query_high_mismatch', which often caused by Bismark mapping error, default=16",default=16, type=int)

    parser.add_argument("--other_mismatch_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'other_high_mismatch', default=12",default=12, type=int)

    parser.add_argument("--total_mismatch_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'total_high_mismatch', default=24",default=24, type=int)

    parser.add_argument("--filter_log",
        help="Alignment filter log file, default=./filter.log",default="./filter.log")

    

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    in_bam_filename = ARGS.input_bam

    query_mut_type_list = ARGS.query_mutation_type.split(",")
    query_mut_min_cutoff = ARGS.query_mutation_min_cutoff
    query_mut_max_cutoff = ARGS.query_mutation_max_cutoff
    other_mut_max_cutoff = ARGS.other_mismatch_max_cutoff
    total_mut_max_cutoff = ARGS.total_mismatch_max_cutoff
    
    total_threads = ARGS.threads

    filter_log_dir = os.path.dirname(ARGS.filter_log)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make full cmd
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    full_cmd_str = "python calculate-mut-stats-V01.py \n \
    --input_bam {input_bam} \n \
    --reference {reference} \n \
    --threads {threads} \n \
    --out_json {out_json} \n \
    --query_mutation_type {query_mutation_type} \n \
    --query_mutation_min_cutoff {query_mutation_min_cutoff} \n \
    --query_mutation_max_cutoff {query_mutation_max_cutoff} \n \
    --other_mismatch_max_cutoff {other_mismatch_max_cutoff} \n \
    --total_mismatch_max_cutoff {total_mismatch_max_cutoff} \n \
    --filter_log {filter_log} \n \
    ".format(
        input_bam = ARGS.input_bam,
        reference = ARGS.reference,
        threads = ARGS.threads,
        out_json = ARGS.out_json,
        query_mutation_type = ARGS.query_mutation_type,
        query_mutation_min_cutoff = ARGS.query_mutation_min_cutoff,
        query_mutation_max_cutoff = ARGS.query_mutation_max_cutoff,
        other_mismatch_max_cutoff = ARGS.other_mismatch_max_cutoff,
        total_mismatch_max_cutoff = ARGS.total_mismatch_max_cutoff,
        filter_log = ARGS.filter_log
    )

    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(full_cmd_str)
    sys.stderr.write("\n" + "-" * 80 + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check input bam file
    check_res = check_input_bam_file(in_bam_filename)

    if check_res == 1:
        raise IOError("BAM not exisit!")

    elif check_res == 2:
        raise IOError("BAM not sorted by coordinate!")

    elif check_res == 3:
        raise IOError("BAM doesn't contain index file!")

    # check other file exist
    if not os.path.exists(ARGS.reference):
        raise IOError("--reference does not exist! ")        

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load ref name
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ref_name_list = get_BAM_ref_name(in_bam_filename, ref_count_cutoff=-1)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # ctrl BAM count 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Processing BAM ... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    pool = multiprocessing.Pool(processes = total_threads)

    input_BAM_filter_log_filename_list = []
    input_BAM_count_result = []

    for index, ref_name in enumerate(ref_name_list):
        # make log filename
        input_log_basename = "input_filter_" + ref_name + "." +  "".join(random.sample(string.ascii_letters + string.digits, 16)) + ".log"
        input_log_filename = os.path.join(filter_log_dir, input_log_basename)
        input_BAM_filter_log_filename_list.append(input_log_filename)

        input_BAM_count_result.append(
            pool.apply_async(
                func = run_BAM_count_on_chromosme, 
                args = (
                    in_bam_filename,
                    ref_name,
                    ARGS.reference,
                    query_mut_type_list,
                    "All",
                    query_mut_min_cutoff,
                    query_mut_max_cutoff,
                    total_mut_max_cutoff,
                    other_mut_max_cutoff,
                    True,
                    input_log_filename,
                )
            )
        )

    pool.close()
    pool.join()

    # out pool result
    input_BAM_count_dict = {}
    for res in input_BAM_count_result:
        run_res = res.get()
        input_BAM_count_dict[run_res[0]] = run_res[1].copy()

    # merge log file
    input_merge_log = merge_log_file(
        in_filename_list = input_BAM_filter_log_filename_list, 
        out_filename = ARGS.filter_log, 
        out_mode="wb", 
        rm_input_file = True
    )

    sys.stderr.write("Processing BAM Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # store dict info 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ARGS.out_json == "stdout":
        out_filename_str = "stdout"
    else:
        out_filename_str = os.path.abspath(ARGS.out_json)

    meta_count_dict = {
        "meta_data" : {
            "create_date" : time.strftime("%Y-%m-%d %H:%M:%S",time.localtime()),
            "input_BAM" : os.path.abspath(in_bam_filename),
            "out_json" : out_filename_str,
            "reference" : os.path.abspath(ARGS.reference),
            "threads" : ARGS.threads,
            "query_mutation_type" : ARGS.query_mutation_type,
            "query_mutation_min_cutoff" : ARGS.query_mutation_min_cutoff,
            "query_mutation_max_cutoff" : ARGS.query_mutation_max_cutoff,
            "other_mismatch_max_cutoff" : ARGS.other_mismatch_max_cutoff,
            "total_mismatch_max_cutoff" : ARGS.total_mismatch_max_cutoff,
            "filter_log" : os.path.abspath(ARGS.filter_log)
        },
        "count_dict" : input_BAM_count_dict.copy()
    }

    # make final output 
    json_obj = json.dumps(meta_count_dict, encoding="utf-8", sort_keys=True, indent=4, separators=(', ', ': '))

    if ARGS.out_json == "stdout":
        sys.stdout.write(json_obj + "\n")

    else:
        with open(ARGS.out_json, "wb") as json_out_file:
            json_out_file.write(json_obj)

# 2019-12-30

