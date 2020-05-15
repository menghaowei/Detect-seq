#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
    2020-04-02
        Filter raw mpmat file with or without header


E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
None

"""
# Learning Part END-----------------------------------------------------------


import argparse
import time
import os
import sys

import pysam 
from Bio import  SeqIO

import numpy as np
import pandas as pd


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

    
def get_region_high_mismatch_state(chr_name, region_start, region_end, bam_obj, genome_dict, 
                                   CT_GA_mismatch_cutoff = 16, 
                                   other_mismatch_cutoff = 12, 
                                   total_mismatch_cutoff = 24,
                                   CT_GA_domiant_cutoff = 3,
                                   high_mis_reads_ratio_cutoff = 0.2):
    """
    INPUT:
        <chr_name>, <region_start>, <region_end>
            set region coordinate information
        
        <bam_obj>
            pysam.AlignmentFile object, BAM file have to be sorted
        
        <genome_dict>
            A dict, key is chr_name, value is chromosome seq with capital alphabet
        
        <CT_GA_mismatch_cutoff>
            Alignment read with CT or GA mismatch larger than this cutoff will be counted as 'high_mismatch_read'
        
        <other_mismatch_cutoff>
            Other type mutation
        
        <total_mismatch_cutoff>
            Total mismatch cutoff 
        
        <CT_GA_domiant_cutoff>
            other_mismatch_cutoff - CT_GA_mismatch_cutoff should smaller than a cutoff 
        
        <high_mis_reads_ratio_cutoff>
            if 'high_mismatch_read' count take more than this cutoff, consider the region as a high mismatch region.
    
    RETURN:
        [True/False, region_total_reads, high_mismatch_reads, normal_reads]
    """
    
    region_start = int(region_start)
    region_end = int(region_end)
    
    if (region_end - region_start) >= 100:
        extend_length = 0
    else:
        extend_length = 50

    # fix region 
    get_region_start = region_start - extend_length
    if get_region_start < 1:
        get_region_start = 1

    # get bam region
    bam_seg_obj = bam_obj.fetch(
        contig = chr_name,
        start = get_region_start,
        end = region_end + extend_length
    )

    # get align initial number
    total_align_count = 0
    high_mis_align_count = 0
    normal_align_count = 0
    region_high_mis_state = False

    for align in bam_seg_obj:
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
            other_mismatch_count = query_mismatch_count(count_dict, query_key_list="All", block_key_list=["CT","GA"])
            CT_GA_mismatch_count = query_mismatch_count(count_dict, query_key_list=["CT","GA"], block_key_list=[])

        else:
            other_mismatch_count = 0
            CT_GA_mismatch_count = 0

        # find high mismatch region 
        if CT_GA_mismatch_count >= CT_GA_mismatch_cutoff:
            high_mis_align_count += 1
            continue

        if other_mismatch_count >= other_mismatch_cutoff:
            high_mis_align_count += 1
            continue

        if (CT_GA_mismatch_count + other_mismatch_count) >= total_mismatch_cutoff:
            high_mis_align_count += 1
            continue

        if (other_mismatch_count - CT_GA_mismatch_count) >= CT_GA_domiant_cutoff:
            high_mis_align_count += 1
            continue

        normal_align_count += 1


    if high_mis_align_count / 1.0 / total_align_count <= 0.2:
        if high_mis_align_count >= 3:
            region_high_mis_state = True 

    # return 
    return([region_high_mis_state, total_align_count, high_mis_align_count, normal_align_count])    



###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="A tool to find significant enrichment mpmat region")

    parser.add_argument("-i", "--mpmat_table",
        help=".mpmat table, which can be generated from <pmat-merge.py> code",required=True)

    parser.add_argument("-o", "--output",
        help="Output file, default=stdout",default="stdout")

    parser.add_argument("-b", "--input_BAM",
        help="Control BAM file",required=True)

    parser.add_argument("-r", "--reference",
        help="Genome FASTA file",required=True)

    parser.add_argument("--CT_GA_mutation_max_cutoff",
        help="Alignment read with CT or GA mismatch larger than this cutoff will be counted as 'high_mismatch_read', default=16",default=16, type=int)

    parser.add_argument("--other_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'high_mismatch_read', default=12",default=12, type=int)

    parser.add_argument("--total_mutation_max_cutoff",
        help="An alignment contain mutation number higher than this, considered as 'high_mismatch_read', default=24",default=24, type=int)

    parser.add_argument("--CT_GA_domiant_cutoff",
        help="<other_mismatch_cutoff> - <CT_GA_mismatch_cutoff> should smaller than a cutoff, default=3",default=3, type=int)

    parser.add_argument("--high_mis_reads_ratio_cutoff",
        help="If 'high_mismatch_read' count take more than this cutoff, consider the region as a high mismatch region. default=0.2",default=0.2, type=float)

    parser.add_argument("--input_mpmat_header",
        help="If contain header in input mpmat file, default=False",default=False, type=bool)

    parser.add_argument("--out_filter_line",
        help="Output high mismatch line, default=False",default=False, type=bool)


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    mpmat_filename = ARGS.mpmat_table
    input_bam_filename = ARGS.input_BAM
    
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # make full cmd
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    full_cmd_str = "python filter-mpmat-high-mismatch.py \n \
    --mpmat_table {mpmat_table} \n \
    --input_BAM {input_BAM} \n \
    --reference {reference} \n \
    --output {output} \n \
    --CT_GA_mutation_max_cutoff {CT_GA_mutation_max_cutoff} \n \
    --other_mutation_max_cutoff {other_mutation_max_cutoff} \n \
    --total_mutation_max_cutoff {total_mutation_max_cutoff} \n \
    --CT_GA_domiant_cutoff {CT_GA_domiant_cutoff} \n \
    --high_mis_reads_ratio_cutoff {high_mis_reads_ratio_cutoff} \n \
    --input_mpmat_header {input_mpmat_header} \n \
    --out_filter_line {out_filter_line} \
    ".format(
        mpmat_table = ARGS.mpmat_table,
        input_BAM = ARGS.input_BAM,
        reference = ARGS.reference,
        output = ARGS.output,
        CT_GA_mutation_max_cutoff = ARGS.CT_GA_mutation_max_cutoff,
        other_mutation_max_cutoff = ARGS.other_mutation_max_cutoff,
        total_mutation_max_cutoff = ARGS.total_mutation_max_cutoff,
        CT_GA_domiant_cutoff = ARGS.CT_GA_domiant_cutoff,
        high_mis_reads_ratio_cutoff = ARGS.high_mis_reads_ratio_cutoff,
        input_mpmat_header = ARGS.input_mpmat_header,
        out_filter_line = ARGS.out_filter_line
    )

    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(full_cmd_str)
    sys.stderr.write("\n" + "-" * 80 + "\n")

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

    input_bam = pysam.AlignmentFile(input_bam_filename, "rb")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load genome as dict
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    sys.stderr.write("Load reference... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    genome_dict = {}
    
    genome_fa =  SeqIO.parse(handle= ARGS.reference, format="fasta")

    for ref in genome_fa:
        sys.stderr.write("Loading...\t" + ref.id + "\n")
        genome_dict[ref.id] = ref.seq.upper()

    sys.stderr.write("Load reference... Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # define list
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    header_list = [
        "chr_name",
        "region_start",
        "region_end",
        "region_site_num",
        "region_mut_site_num",
        "region_SNP_mut_num",
        "region_site_index",
        "mut_base_num",
        "cover_base_num",
        "mut_ratio",
        "SNP_ann",
        "tandem_info",
        "Pass_info"
    ]

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load mpmat 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ARGS.input_mpmat_header:
        mpmat_df = pd.read_csv(mpmat_filename,sep="\t")
    
    else:
        mpmat_df = pd.read_csv(mpmat_filename,sep="\t", header=None)

    # fix column name 
    column_rename_dict = {}
    for index, new_name in enumerate(header_list[:-1]):
        column_rename_dict[mpmat_df.columns[index]]= new_name

    mpmat_df.rename(columns = column_rename_dict,inplace=True)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # set output 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ARGS.output == "stdout":
        output_file = sys.stdout
    else:
        output_file = open(ARGS.output,"wb")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # filter mpmat file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    region_state_list = []
    filter_line_num = 0

    for index, row in mpmat_df.iterrows():
        if index % 100 == 0:
            sys.stderr.write("Processed %s \t %s \n" % (index,time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        region_state = get_region_high_mismatch_state(
            chr_name = row["chr_name"], 
            region_start = row["region_start"],
            region_end = row["region_end"], 
            bam_obj= input_bam, 
            genome_dict= genome_dict,
            CT_GA_mismatch_cutoff = ARGS.CT_GA_mutation_max_cutoff,
            other_mismatch_cutoff = ARGS.other_mutation_max_cutoff,
            total_mismatch_cutoff = ARGS.total_mutation_max_cutoff,
            CT_GA_domiant_cutoff = ARGS.CT_GA_domiant_cutoff,
            high_mis_reads_ratio_cutoff = ARGS.high_mis_reads_ratio_cutoff
        )
        
        region_state_list.append(region_state)

        # output row 
        if region_state[0] == True:
            if ARGS.out_filter_line == False:
                output_file.write("\t".join(map(str,list(row))) + "\n")

        elif region_state[0] == False:
            if ARGS.out_filter_line == True:
                output_file.write("\t".join(map(str,list(row))) + "\n")

            filter_line_num += 1

    # output log
    sys.stderr.write("Total input line: %s \n" % len(region_state_list))
    sys.stderr.write("Filter line: %s \n" % filter_line_num)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # close file 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    input_bam.close()

    if ARGS.output != "stdout":
        output_file.close()








            


# 2020-04-02

