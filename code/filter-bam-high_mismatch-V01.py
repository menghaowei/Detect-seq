#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
  2019-12-17 filter BAM with high mismatches

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""

Design pipeline:

1. filter BAM file by chromosome

2. merge BAM file 

3. make BAM file index


"""
# Learning Part END-----------------------------------------------------------


import argparse
import os
import random
import string
import sys
import time
import multiprocessing
import subprocess

import pysam 
from Bio import  SeqIO


###############################################################################
# function part 
###############################################################################
def back_indel_shift(info_index_list, cur_index):
    """
    INPUT:
        <info_index_list>
            generated from align.cigartuples
        
        <cur_index>
            index related to MD tag in BAM file
    
    RETURN
        <acc_shift>
    """
    
    # parse softclip and insertion 
    if len(info_index_list) == 0:
        return 0 

    acc_shift = 0
    for info_start_index, info_len in info_index_list:
        if info_start_index >= cur_index:
            return acc_shift

        else:
            acc_shift += info_len

    return acc_shift


def get_align_mismatch_pairs(align):
    """
    INPUT
        <align>
            pysam AlignedSegment object
    
    RETURN
        <mismatch_pair_list>
            [ref_index, align_index, ref_base, align_base]
            
            When NM == 0, return None
    """
    # No mismatch 
    if align.get_tag("NM") == 0:
        return None
    
    # parse softclip, insertion and deletion
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
            cur_index += int(cur_base)
            cur_base = ""

            if  base == "^":
                i += 1 
                del_str = ""

                while (bases[i].isalpha()) and (i < len(bases)):
                    del_str += bases[i]
                    i += 1

                cur_index += len(del_str)
                del_str = ""

            elif base.isalpha():
                cur_index += 1 
                ref_base = base
                i += 1 
                
                # add into list
                fix_index = cur_index + back_indel_shift(info_index_list, cur_index)                 

                if fix_index < len(align.query_sequence):
                    mismatch_pair_list.append([cur_index-1, fix_index-1, ref_base, align.query_sequence[fix_index-1]])
                else:
                    return(None)
                
    return(mismatch_pair_list)


def get_No_MD_align_mismatch_pairs(align, genome_dict):
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
            
            When NM == 0, return None
            
    """
    # No mismatch 
    if align.get_tag("NM") == 0:
        return None
    
    mismatch_pair_list = []
    for align_idx, ref_idx in align.get_aligned_pairs():
        if (align_idx != None) and (ref_idx != None):
            align_base = align.query_sequence[align_idx]
            ref_base = genome_dict[align.reference_name][ref_idx]
            
            if align_base != ref_base:
                mismatch_pair_list.append([ 
                    ref_idx - align.reference_start,
                    align_idx,
                    align_base,
                    ref_base
                ])
                
    return(mismatch_pair_list)
                 

def make_mismatch_count_dict(align_mismatch_pairs):
    """
    RETURN
        <mismatch_count_dict>
            key like CT means C to T mismatch 
    """
    count_dict = {}
    
    for ref_idx, align_idx, ref_base, align_base in align_mismatch_pairs:
        key = ref_base + align_base
        
        if count_dict.get(key) == None:
            count_dict[key] = 1
            
        else:
            count_dict[key] += 1
        
    return(count_dict)


def query_mismatch_count(count_dict, query_key_list = "All", block_key_list=None, ignore_key=True):
    """
    HELP:
    
        <count_dict> 
            generate from Function <make_mismatch_count_dict>
        
        <query_key>
            default = "All", return all mismatch count 
            or only return sum of  select key count in query_key_list
        
        <block_key_list>
            the key in this list will not be counted
            
        <ignore_key>
            True means if key not in dict return 0 as count, 
            False means if key not in dict raise an IOError
    
    """
    mis_count = 0
    
    if query_key_list == "All":
        query_key_list = count_dict.keys()
    
    for key in query_key_list:
        if key not in block_key_list:
            count = count_dict.get(key)
            
            if count != None:
                mis_count += count 
    
    return(mis_count)


def get_BAM_ref_name(bam_filename, ref_count_cutoff = 0):
    """
    INPUT:
        <bam_filename>
            BAM file path with .bai index file in a same dir
        
    RETURN
        <ref_name_list>
            return ref name list, which contain more than ref_count_cutoff mapped reads
    
    """
    bam_file = pysam.AlignmentFile(bam_filename, "rb")
    ref_name_list = []
    
    for bam_index_info in bam_file.get_index_statistics():
        ref_name = bam_index_info[0]
        map_count = bam_index_info[1]
        unmap_count = bam_index_info[2]
        
        if map_count > ref_count_cutoff:
            if ref_name != "*":
                ref_name_list.append(ref_name)
    
    bam_file.close()
    
    return(ref_name_list)
    

def run_BAM_filter_on_chromosme(in_bam_filename, out_bam_filename, select_chr_name, block_mismatch_type_list, query_mismatch_type_list = "All", mismatch_filter_cutoff = 5):
    """
    INPUT:
        <bam_filename>

        <out_bam_filename>

        <chr_name>

        <block_mismatch>

        <query_mismatch>
    
    RETURN:
        if everything word well return:
            [total_align_count, filter_align_count, final_align_count]
                > total_align_count: total align count 
                > filter_align_count: align with too many mismatches
                > final_align_count: final out align count
            
        something wrong return:
            None
    """
    
    in_bam_file = pysam.AlignmentFile(in_bam_filename, "rb")
    out_bam_file = pysam.AlignmentFile(out_bam_filename, "wb", template = in_bam_file)

    total_align_count = 0
    filter_align_count = 0
    final_align_count = 0

    for align in in_bam_file.fetch(reference = select_chr_name):
        # count 
        total_align_count += 1

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
            align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, genome_dict)

        # count mismatch info
        if align_mismatch_pairs != None:
            count_dict = make_mismatch_count_dict(align_mismatch_pairs)
            mismatch_count = query_mismatch_count(count_dict, query_key_list="All", block_key_list=block_mismatch_type_list)

            if mismatch_count <= mismatch_filter_cutoff:
                out_bam_file.write(align)
                final_align_count += 1

            else:
                filter_align_count += 1 

        else:
            out_bam_file.write(align)
            final_align_count += 1
    
    # close file 
    in_bam_file.close()
    out_bam_file.close()

    # return 
    return ([select_chr_name, total_align_count, filter_align_count, final_align_count])
    

def check_input_bam_file(input_bam_filename):
    """
    HELP
    
        1. check exisit
        
        2. check sort state 
        
        3. check .bai index file 
    
    RETURN
        0 alright
        
        Not 0:
            1 not exist 
            2 not sorted 
            3 index not exist 
    """
    # check exist 
    if not os.path.exists(input_bam_filename):
        return 1 
    
    # check sort 
    bam_file = pysam.AlignmentFile(input_bam_filename, "rb")
    if bam_file.header.get("HD").get("SO") != "coordinate":
        return 2
    
    # check bai index 
    bai_filename_1 = input_bam_filename + ".bai"
    bai_filename_2 = os.path.splitext(input_bam_filename)[0] + ".bai"
    if (not os.path.exists(bai_filename_1)) and (not os.path.exists(bai_filename_2)):
        return 3 
    
    return 0

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Filter BAM high mismatches alignment")

    parser.add_argument("-i", "--input_bam",
        help="Input BAM file, have to sorted by coordinate and indexed",required=True)

    parser.add_argument("-r", "--reference",
        help="Genome FASTA file",required=True)

    parser.add_argument("-o", "--output_bam",
        help="Output BAM file",required=True)

    parser.add_argument("-p", "--threads",
        help="Multiple threads number, default=1",default=1, type=int)

    parser.add_argument("--mismatch_cutoff",
        help="Reads in BAM file with more than <mismatch_cutoff> mismatch will be discarded. default=5",default=5, type=int)

    parser.add_argument("--block_mismatch_type",
        help="Reads with those mismach types will be not counted, default=CT,GA",default="CT,GA")

    parser.add_argument("--samtools_path",
        help="samtools path, default considered can run at current PATH environment",default="samtools")


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    input_bam_filename = ARGS.input_bam
    output_bam_filename = ARGS.output_bam

    block_mismatch_type_list = ARGS.block_mismatch_type.split(",")
    mismatch_cutoff = ARGS.mismatch_cutoff
    total_threads = ARGS.threads

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # check input bam file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    check_res = check_input_bam_file(input_bam_filename)

    if check_res == 1:
        raise IOError("BAM not exisit!")

    elif check_res == 2:
        raise IOError("BAM not sorted by coordinate!")

    elif check_res == 3:
        raise IOError("BAM doesn't contain index file!")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make temp files list
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    ref_name_list = get_BAM_ref_name(input_bam_filename, ref_count_cutoff=0)
    temp_filename_list = []

    for ref_name in ref_name_list:
        temp_file_basename = "temp_" + ref_name + "." +  "".join(random.sample(string.ascii_letters + string.digits, 16)) + ".bam"
        temp_filename = os.path.join(os.path.dirname(output_bam_filename), temp_file_basename)
        temp_filename_list.append(temp_filename)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load reference fasta 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    sys.stderr.write("Load reference... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    genome_dict = {}
    genome_fa =  SeqIO.parse(handle= ARGS.reference,format="fasta")

    for ref in genome_fa:
        if ref.id in ref_name_list:
            sys.stderr.write("Loading...\t" + ref.id + "\n")
            genome_dict[ref.id] = ref.seq.upper()

        else:
            sys.stderr.write("Skip...\t" + ref.id + "\n")

    sys.stderr.write("Load reference... Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # run filter 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Filter BAM ... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    pool = multiprocessing.Pool(processes = total_threads)

    filter_result = []

    for index, ref_name in enumerate(ref_name_list):
        out_temp_bam_filename = temp_filename_list[index]
        
        run_result = (
            pool.apply_async(
                func = run_BAM_filter_on_chromosme, 
                args = (
                    input_bam_filename,
                    out_temp_bam_filename,
                    ref_name,
                    block_mismatch_type_list,
                    "All",
                    mismatch_cutoff,
                )
            )
        )

        filter_result.append(run_result)

    pool.close()
    pool.join()

    sys.stderr.write("Filter BAM Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # out log
    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("Filter Info:\n")
    sys.stderr.write("FI\tchr_name\ttotal_align_count\tfilter_align_count\tfinal_align_count\n")

    for res in filter_result:
        log_result = res.get()
        log_str = "\t".join(map(str,log_result))
        sys.stderr.write("FI\t" + log_str + "\n")

    sys.stderr.write("-" * 80 + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # merge files
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Merge temp file... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    temp_bam_str = " ".join(temp_filename_list)

    merge_cmd = "{samtools_path} cat -o {output} {temp_bam}".format(
        samtools_path = ARGS.samtools_path,
        output = output_bam_filename,
        temp_bam = temp_bam_str
    )
    
    # Log
    sys.stderr.write("Merge BAM command:\n%s\n" % merge_cmd)

    # run bedtools 
    merge_runcode = subprocess.call(merge_cmd, shell=True)

    # return 
    if merge_runcode != 0:
        sys.stderr.write("Merge temp file... Error! Unsuccessfully!\t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    else:
        sys.stderr.write("Merge temp file... Done!\t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))
         
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # rm temp files
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    for temp_filename in temp_filename_list:
        os.remove(temp_filename)

# 2019-12-17