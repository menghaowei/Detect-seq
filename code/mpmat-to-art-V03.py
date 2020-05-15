#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""

Author: MENG Howard

Version-03:
    2020-03-31
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
LEARNING_PART = \
"""

There will be a lot of parameters to set,

but default paramters can handle output figure in most cases.

"""


PIPELINE = \
"""

INPUT
    <mpmat>
    
    <fasta>
    
    <alignment settings>
    
    <id_prefix>
    
OUTPUT
    <alignment_result_table>
        format like:
            chr_name
            region_start
            region_end
            region_id
            region_site_index
            SNP_ann
            Pass_info
            align_direction
            PAM_type
            PAM_seq
            align_total_match
            align_total_mismatch
            align_total_gap
            align_score
            seed_5bp_mismatch
            seed_5bp_gap
            head_5bp_mismatch
            head_5bp_gap
            head_8bp_mismatch
            head_8bp_gap
            align_target_seq
            align_info_state
            align_query_seq
            align_chr_name
            align_chr_start
            align_chr_en
            ... ...
            
DESIGN
    0. use bedtools getfasta to get fasta temp files
    
    1. read all mpmat and fa file 
    
    2. all fa file are Watson strand seq
    
    3. alignment and back alignment info (fwd / rev)
    
    4. calculate output info (.art) file

"""
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import sys
import os 
import time
import subprocess

import pandas as pd
import numpy as np
import string
import random
from Bio import SeqIO
import Bio

np.set_printoptions(suppress=True)

###############################################################################
# function part 
###############################################################################
def find_seq_PAM_index(query_seq):
    """
    <INPUT>
        query_seq
    
    <HELP>
        e.g. query_seq = "AAAGAGAG"
        return index list = [1,3,5], which means search NAG, NGG at the same time
    """
    pam_index_list = []
    
    for index, base in enumerate(query_seq[:-2]):
        if query_seq[index+1:index+3] == "AG":
            pam_index_list.append((index, "NAG"))
            
        elif query_seq[index+1:index+3] == "GG":
            pam_index_list.append((index, "NGG"))
            
    return(pam_index_list)


def analysis_align_obj(alignment, reverse_state = False):
    """
    INPUT:
        <alignment obj>
        
    OUTPUT:
        <info> 
            1. match count 
            2. mismatch count 
            3. gap count 
            4. alignment.score
            
    HELP:
        2019-11-15 fix-1
            The gap count should be the num of gap contain in sgRNA alignment region.
            
            e.g.
            
            AGTGGTAAGAAGAAGACGAGACATAATGAG
            ------||||||||||||||X|----||||
            ------AAGAAGAAGACGAGCC----TGAG
            
            gap count should be 4, rather than 10.
        
        2019-11-15 fix-2
            add return info, start_index, end_index, now the retrun list will be 
            
            return_list = [
                match_count,
                mismatch_count,
                gap_count,
                alignment.score,
                start_index,
                end_index
            ]
            
            The <start_index> and <end_index> are index related to sgRNA alignment string
        
    """
    
    # define params 
    match_count = 0
    mismatch_count = 0
    gap_count = 0

    alignment_list =  str(alignment).split("\n")
    query_length = len(alignment.query)    
    
    if reverse_state:
        target_list = alignment_list[0][::-1]
        info_list = alignment_list[1][::-1]
        query_list = alignment_list[2][::-1]

    else:
        target_list = alignment_list[0]
        info_list = alignment_list[1]
        query_list = alignment_list[2]        
        
    count_state = False
    count_query_base = 0
    
    ref_align_start_index = 0
    ref_align_end_index = len(alignment.target) - 1
    
    # counting 
    for index, info_base in enumerate(info_list):
        if not count_state:
            if query_list[index] != "-":
                count_state = True
                ref_align_start_index = index
                
        if not count_state:
            continue
        
        else:
            if info_base == "|":
                match_count += 1
                count_query_base += 1

            elif info_base == "X":
                mismatch_count += 1
                count_query_base += 1

            elif info_base == "-":
                gap_count += 1
                if query_list[index] != "-":
                    count_query_base += 1 
            
        if count_query_base >= query_length:
            ref_align_end_index = index
            break
                
    return_list = [
        match_count,
        mismatch_count,
        gap_count,
        alignment.score,
        ref_align_start_index,
        ref_align_end_index
    ]
    
    return(return_list)


def sign_value(x):
    """
    HELP
        sign function
    """
    if x < 0:
        return(-1)
    elif x == 0:
        return(0)
    else:
        return(1)    


def cmp_align_list(align_a, align_b):
    """
    INPUT
        like [17, 3, 0, 73.0, 1, 22 'AAGAAGAAGACGAGTCTGCA', '||||||||||||||X|||XX', 'AAGAAGAAGACGAGCCTGAG']
    
    HELP
        compare function for align_list 
    """
    align_alphabet = {"-":0, "X":1, "|":2}
    sort_index_list = [3, 2, 1, 0]
    sort_rev_state_list = [True, False, False, True]    
    
    for order_index, align_index in enumerate(sort_index_list):
        if (align_a[align_index]  - align_b[align_index]) == 0:
                continue
        else:
            if sort_rev_state_list[order_index]:
                return( -1 * sign_value(align_a[align_index]  - align_b[align_index]))
            else:
                return( sign_value(align_a[align_index]  - align_b[align_index]))
    
    for index, char_a in enumerate(align_a[7]):
        if index <= (len(align_b[7]) - 1):
            value_a = align_alphabet[char_a]
            value_b = align_alphabet[align_b[7][index]]

            if value_a > value_b:
                return 1
            elif value_a < value_b:
                return -1 

    return 0


def run_sgRNA_alignment(align_ref_seq, align_sgRNA, sgRNA_aligner, extend_len=3):
    """
    INPUT
        <align_ref_seq>
            
        <align_sgRNA> 
            sgRNA seq without PAM
            
        <possible_sgRNA_region>
    
    RETURN
        <final_align_res_list>
    """
    align_sgRNA_rev = align_sgRNA[::-1]
    
    # find all PAM
    PAM_info_list = find_seq_PAM_index(align_ref_seq)
    if len(PAM_info_list) == 0:
        return([])
    
    # forward alignment 
    final_align_res_list = []

    for PAM_start_idx, PAM_type in PAM_info_list:

        # filter 5' end PAM 
        if (PAM_start_idx - extend_len) < len(align_sgRNA):
            continue

        # select PAM and sgRNA in the possible region 
        region_seq_start = PAM_start_idx-len(align_sgRNA) - extend_len
        region_seq_end = PAM_start_idx + 3 

        # alignment part 
        region_seq = align_ref_seq[region_seq_start : PAM_start_idx]
        align_res = sgRNA_aligner.align(region_seq[::-1], align_sgRNA_rev)

        # parse alignment
        ## if contain multiple alignment result with the same score, keep the best one;
        ## sort reason score -> gap -> mismatch -> match 
        align_res_list = []
        for align in align_res:
            align_analysis_res = analysis_align_obj(align, reverse_state=True)
            align_info_list = [ x[::-1] for x in str(align).strip().split("\n")]
            align_analysis_res += align_info_list
            align_analysis_res += [PAM_start_idx, PAM_type]
            align_res_list.append(align_analysis_res)

        if len(align_res_list) == 1:
            final_align_res_list.append(align_res_list[0])

        else:
            align_res_list.sort(cmp=cmp_align_list)
            final_align_res_list.append(align_res_list[0])

    # sort final alignment 
    final_align_res_list.sort(cmp=cmp_align_list)
    
    return(final_align_res_list)


def run_no_PAM_sgRNA_alignment_no_chop(align_ref_seq, align_sgRNA_full, no_PAM_sgRNA_aligner):
    """
    INPUT
        <align_ref_seq>
            
        <align_sgRNA> 
            sgRNA seq without PAM
            
        <no_PAM_sgRNA_aligner>
            An obj from BioPython pairwise alignment
    
    RETURN
        <final_align_res_list>
    """

    # alignment part 
    align_res = no_PAM_sgRNA_aligner.align(align_ref_seq, align_sgRNA_full)

    # parse alignment
    ## if contain multiple alignment result with the same score, keep the best one;
    ## sort reason score -> gap -> mismatch -> match 
    align_res_list = []
    for align in align_res:
        align_analysis_res_temp = analysis_align_obj(align, reverse_state=False)
        align_analysis_res = align_analysis_res_temp[:]
        align_analysis_res += str(align).strip().split("\n")

        PAM_start_index = align_analysis_res[5] - 2
        PAM_type = align_ref_seq[PAM_start_index : PAM_start_index+3]

        align_analysis_res += [PAM_start_index, PAM_type]
        align_res_list.append(align_analysis_res)

    if len(align_res_list) == 1:
        return (align_res_list)

    else:
        align_res_list.sort(cmp=cmp_align_list)
        return (align_res_list)


def get_mut_region_strand_direction(site_index_list_str):
    """
    INPUT
        <site_index_list_str> 
            A str like 'chr1_143265617_CT,chr1_143265637_CT,chr1_143265643_CT,chr1_143265655_CT,chr1_143265663_CT,chr1_143265666_CT'
    
    RETRUN
        count CT and GA site number then return CT / GA mutation region tyep
    """
    
    site_index_list = site_index_list_str.split(",")
    CT_count = GA_count = 0 
    
    for site_index in site_index_list:
        if "CT" in site_index:
            CT_count += 1
        
        elif "GA" in site_index:
            GA_count += 1
    
    if CT_count >= GA_count:
        return("CT")
    
    else:
        return("GA")


def back_align_res_stats(align_direction, align_res_list, align_ref_seq, align_sgRNA_seq_full):
    """
    INPUT
        <align_direction>
            PAM_fwd, PAM_rev, no_PAM_fwd, no_PAM_rev
        
        <align_res_list>
            format like 
                [20, 0, 0, 100.0, 0, 19, 'CCTGAGTCCGAGCAGAAGAAGAA', '---||||||||||||||||||||', '---GAGTCCGAGCAGAAGAAGAA', 58, 'NGG']
            
            values are
                match_num
                mismatch_num
                gap_num
                alignment_score
                align_start_index
                align_end_index
                align_ref_seq
                align_state
                align_query_seq
                PAM_index
                PAM_info
            
            Warning:
            
                When <align_direction> is no_PAM_fwd, no_PAM_rev, the align_query_seq is full seq of sgRNA, 
            
                Which contain NGG/NAG PAM info.
                
        RETURN:

            align_direction
            PAM_type
            PAM_seq
            align_total_match
            align_total_mismatch
            align_total_gap
            align_score
            seed_5bp_mismatch
            seed_5bp_gap
            head_5bp_mismatch
            head_5bp_gap
            head_8bp_mismatch
            head_8bp_gap
            align_target_seq
            align_info_state
            align_query_seq
            
    """

    align_dict = {}
    
    # initialize
    align_dict["align_direction"] = align_direction
    align_dict["align_total_match"] = align_res_list[0]
    align_dict["align_total_mismatch"] = align_res_list[1]
    align_dict["align_total_gap"] = align_res_list[2]
    align_dict["align_score"] = align_res_list[3]   

    align_dict["align_target_seq"] = align_res_list[6][align_res_list[4]: align_res_list[5]+1]
    align_dict["align_info_state"] = align_res_list[7][align_res_list[4]: align_res_list[5]+1]
    align_dict["align_query_seq"] = align_res_list[8][align_res_list[4]: align_res_list[5]+1]

    align_dict["seed_5bp_mismatch"] = 0
    align_dict["seed_5bp_gap"] = 0

    align_dict["head_5bp_mismatch"] = 0
    align_dict["head_5bp_gap"] = 0

    align_dict["head_8bp_mismatch"] = 0
    align_dict["head_8bp_gap"] = 0

    # set ref_seq
    align_ref_seq_rc = Bio.Seq.reverse_complement(align_ref_seq_rc)
    align_sgRNA_seq = align_sgRNA_seq_full[-3:]

    count_index = 1
    for index in xrange(align_res_list[4], align_res_list[5]+1):
        align_info_char = align_res_list[7][index]

        if align_info_char == "X":
            if count_index <= 5:
                align_dict["head_5bp_mismatch"] += 1 

            if count_index <= 8:
                align_dict["head_8bp_mismatch"] += 1 

            if (len(align_sgRNA_seq) - count_index + 1) <= 5:
                align_dict["seed_5bp_mismatch"] += 1 

        elif align_info_char == "-":
            if count_index <= 5:
                align_dict["head_5bp_gap"] += 1 

            if count_index <= 8:
                align_dict["head_8bp_gap"] += 1 

            if 0 <= (len(align_sgRNA_seq) - count_index + 1) <= 5:
                align_dict["seed_5bp_gap"] += 1 

        count_index +=1 


    if (final_align_direction == "PAM_fwd") or (final_align_direction == "PAM_rev"):
        align_dict["PAM_type"] = align_res_list[10]

        if final_align_direction == "PAM_fwd":
            align_dict["PAM_seq"] = align_ref_seq[align_res_list[9]:align_res_list[9]+3]

        elif final_align_direction == "PAM_rev":
            align_dict["PAM_seq"] = align_ref_seq_rc[align_res_list[9]:align_res_list[9]+3]

        # add PAM info 
        align_dict["align_target_seq"] += align_dict["PAM_seq"]
        align_dict["align_query_seq"] += align_sgRNA_seq_full[-3:]

        for char_idx, target_char in enumerate(align_dict["PAM_seq"]):
            if target_char == align_sgRNA_seq_full[len(align_sgRNA_seq)+char_idx:len(align_sgRNA_seq)+char_idx+1]:
                align_dict["align_info_state"] += "|"
            else:
                align_dict["align_info_state"] += "X"

        
    elif (final_align_direction == "No_PAM_fwd") or (final_align_direction == "No_PAM_rev"):
        align_dict["PAM_type"] = "NoPAM"
        align_dict["PAM_seq"] = align_res_list[10]

    # return part 
    return_list = [
        align_dict["align_direction"],
        align_dict["PAM_type"],
        align_dict["PAM_seq"],
        align_dict["align_total_match"],
        align_dict["align_total_mismatch"],
        align_dict["align_total_gap"],
        align_dict["align_score"],
        align_dict["seed_5bp_mismatch"],
        align_dict["seed_5bp_gap"],
        align_dict["head_5bp_mismatch"],
        align_dict["head_5bp_gap"],
        align_dict["head_8bp_mismatch"],
        align_dict["head_8bp_gap"],
        align_dict["align_target_seq"],
        align_dict["align_info_state"],
        align_dict["align_query_seq"]
    ]

    return(return_list)


def make_temp_FASTA_file(mpmat_table_df, extend_size, sample_prefix, out_dir, reference_fa, bedtools_path):
    """
    INPUT:
        <sample_prefix>
        <out_dir>
        <reference_fa>
        <bedtools_path>

    RETURN:
        If successed return
            [0, temp_FASTA_filename]

        If unsuccessed return
            [1, None]
    """

    # extend region start, end
    region_fix_start = mpmat_table_df["region_start"] - extend_size
    region_fix_end = mpmat_table_df["region_end"] + extend_size

    # add region id
    region_id = []
    for index in mpmat_table_df.index:
        region_id_list = [
            sample_prefix,
            mpmat_table_df["chr_name"][index], 
            mpmat_table_df["region_start"][index], 
            mpmat_table_df["region_site_num"][index]
        ]
        
        if "CT" in mpmat_table_df["region_site_index"][index]:
            region_id_list.append("CT")
        else:
            region_id_list.append("GA")
            
        region_id_list.append(index + 1)
        region_id.append( "_".join( map(str, region_id_list)) )

    mpmat_table_df["region_id"] = region_id

    # make out bed 
    region_bed_table = pd.DataFrame()
    region_bed_table["extend_chr_name"] = mpmat_table_df["chr_name"]
    region_bed_table["extend_chr_start"] = region_fix_start
    region_bed_table["extend_chr_end"] = region_fix_end
    region_bed_table["id"] = region_id

    # make temp BED file 
    temp_bed_out_dir = out_dir
    temp_bed_basename =  "temp_" + sample_prefix + "_"+ "".join(random.sample(string.ascii_letters + string.digits, 16)) + ".bed"
    temp_bed_filename = os.path.join(temp_bed_out_dir, temp_bed_basename)
    region_bed_table.to_csv(temp_bed_filename, sep="\t", header=False, index=False)

    # Log
    sys.stderr.write("Output temp bed with filename: %s\n" % temp_bed_filename)

    # use bedtools to make fasta file
    out_fasta_filename = temp_bed_filename + ".fa"
    genome_ref_fasta = ARGS.reference
    tool_bedtools_path = ARGS.bedtools_path

    tool_bedtools_cmd_fmt = "{bedtools} getfasta -fi {ref_fasta}  -fo {out_fasta} -bed {input_bed} -name+"
    tool_bedtools_cmd = tool_bedtools_cmd_fmt.format(
        bedtools = tool_bedtools_path,
        ref_fasta = genome_ref_fasta,
        out_fasta = out_fasta_filename,
        input_bed = temp_bed_filename
    )

    # Log
    sys.stderr.write("bedtools command:\n%s\n" % tool_bedtools_cmd)

    # run bedtools 
    get_fasta_runcode = subprocess.call(tool_bedtools_cmd, shell=True)

    # rm temp file
    os.remove(temp_bed_filename)

    # return 
    if get_fasta_runcode != 0:
        os.remove(out_fasta_filename)
        sys.stderr.write("Error! Unsuccessfully get mpmat region FASTA seq!\n")
        return([1, None ,None])

    else:
        sys.stderr.write("Successfully get mpmat region FASTA seq!\n")
        return([0, out_fasta_filename, region_bed_table])


def align_position_on_chrom(final_align, final_align_direction, ref_seq_chr_name, ref_seq_chr_start, ref_seq_chr_end, ref_seq_len):
    """
    INPUT
        <final_align>
        
        <final_align_direction>

        <ref_seq_chr_start>
            chromsome position 

        <ref_seq_chr_end>
            chromsome position 
            
    RETURN
        a list like:
        [ <align_pos_on_chr_name>, <align_pos_on_chr_start>, <align_pos_on_chr_end> ]
        
    """
    if final_align_direction == "PAM_fwd":
        end_index = final_align[9] + 3 
        start_index = final_align[9] - len(final_align[6]) + final_align[6].count("-") + final_align[4]

        align_pos_on_chr_name = ref_seq_chr_name
        align_pos_on_chr_start = ref_seq_chr_start + start_index
        align_pos_on_chr_end = ref_seq_chr_start + end_index - 1 

    elif final_align_direction == "PAM_rev":        
        end_index = final_align[9] + 3 
        start_index = final_align[9] - len(final_align[6]) + final_align[6].count("-") + final_align[4]

        align_pos_on_chr_name = ref_seq_chr_name
        align_pos_on_chr_start = ref_seq_chr_start + ref_seq_len - end_index
        align_pos_on_chr_end = ref_seq_chr_start + ref_seq_len - start_index - 1 

    elif final_align_direction == "No_PAM_fwd":
        start_index = final_align[4]
        end_index = final_align[5]        

        align_pos_on_chr_name = ref_seq_chr_name
        align_pos_on_chr_start = ref_seq_chr_start + start_index
        align_pos_on_chr_end = ref_seq_chr_start + end_index
        
    elif final_align_direction == "No_PAM_rev":
        start_index = final_align[4]
        end_index = final_align[5]        

        align_pos_on_chr_name = ref_seq_chr_name
        align_pos_on_chr_start = ref_seq_chr_start + ref_seq_len - end_index - 1 
        align_pos_on_chr_end = ref_seq_chr_start + ref_seq_len - start_index - 1 
        
    return([align_pos_on_chr_name, align_pos_on_chr_start, align_pos_on_chr_end])
    
###############################################################################
# default sgRNA
###############################################################################

default_sgRNA_dict = {
    "VEGFA": "GACCCCCTCCACCCCGCCTCCGG", 
    "EMX1": "GAGTCCGAGCAGAAGAAGAAGGG", 
    "HEK3": "GGCCCAGACTGAGCACGTGATGG", 
    "HEK4":"GGCACTGCGGCTGGAGGTGGGGG", 
    "RNF2":"GTCATCTTAGTCATTACCTGAGG"
}   


###############################################################################
# main part
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="From .mpmat to alignment result table (.art) file")

    parser.add_argument("-i", "--mpmat_table",
        help=".mpmat table file, generated from <pmat-merge.py> or <mpmat-select.py>",required=True)

    parser.add_argument("--sgRNA",
        help="sgRNA sequence with PAM (NGG/NAG) sequence",required=True)

    parser.add_argument("-o","--out_alignment_result",
        help="Output alignment result table (.art) filename, default=Stdout",default="Stdout")

    parser.add_argument("-r","--reference",
        help="Reference genome fasta file",required=True)

    parser.add_argument("-s","--sample_name",
        help="Sample name of this .mpmat, default=run_mpmat",default="run_mpmat")

    parser.add_argument("--align_dist_to_signal",
        help="If distance too far between mutation signal and alignment, consider as there is no appropriate alignment, default=20",default="20")

    parser.add_argument("--bedtools_path",
        help="Software <bedtools> PATH, default considered <bedtools> is already included in current PATH environment.",default="bedtools")

    parser.add_argument("--align_settings",
        help="Set <align_match_score>  <align_mismatch_score> <align_gap_open_score> <align_gap_extension_score>, default=5,-4,-24,-8",default="5,-4,-24,-8")

    parser.add_argument("--align_min_score",
        help="If alignment score lower than this, consider as no appropriate alignment, default=15",default="15")

    parser.add_argument("--input_header",
        help="If .mpmat file contain header, default=False",default="False")

    parser.add_argument("--input_sep",
        help=r"default=\t",default=r"\t")

    parser.add_argument("--more_colname",
        help="More info you want to include in output table, so you can write the column names like col1,col2... default=None",default="None")

    ARGS = parser.parse_args() 

    # ------------------------------------------------------------------------>>>>>>>
    # load sgRNA
    # ------------------------------------------------------------------------>>>>>>> 
    # set sgRNA info
    if ARGS.sgRNA in default_sgRNA_dict.keys():
        sgRNA_seq_full = default_sgRNA_dict[ARGS.sgRNA]
    else:
        sgRNA_seq_full = str(ARGS.sgRNA).upper()
    
    sgRNA_seq = sgRNA_seq_full[:-3]
    sgRNA_PAM_info = sgRNA_seq_full[-3:]
    
    # ------------------------------------------------------------------------>>>>>>>
    # load mpmat table 
    # ------------------------------------------------------------------------>>>>>>> 
    input_header_state = eval(ARGS.input_header)

    input_header_list = [
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

    # load table and set column name
    if not input_header_state:
        mpmat_table = pd.read_csv(ARGS.mpmat_table, sep=ARGS.input_sep, header=None)
        mpmat_table.columns= input_header_list[:len(mpmat_table.columns)]
    else:
        mpmat_table = pd.read_csv(ARGS.mpmat_table, sep=ARGS.input_sep)
        # print mpmat_table[:10]
        # mpmat_table.columns= input_header_list[:len(mpmat_table.columns)]

    # change data type 
    mpmat_table[["region_start"]].astype(int)
    mpmat_table[["region_end"]].astype(int)

    # ------------------------------------------------------------------------>>>>>>>
    # get fasta file
    # ------------------------------------------------------------------------>>>>>>> 
    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("From mpmat file to generate FASTA file...\n")

    signal_to_align_dist_cutoff = int(ARGS.align_dist_to_signal)
    extend_size = len(sgRNA_seq_full) + signal_to_align_dist_cutoff
    sample_prefix = ARGS.sample_name

    # get out dir
    if ARGS.out_alignment_result == "Stdout":
        temp_bed_out_dir = os.path.abspath(".")
    else:
        temp_bed_out_dir = os.path.dirname(ARGS.out_alignment_result)  

    # run bedtools
    make_fa_state, fasta_filename, region_bed_table = make_temp_FASTA_file(
        mpmat_table, 
        extend_size, 
        sample_prefix, 
        temp_bed_out_dir, 
        ARGS.reference, 
        ARGS.bedtools_path)
    
    if make_fa_state != 0:
        raise IOError("Make FASTA file unsuccessed!")
    
    # ------------------------------------------------------------------------>>>>>>>
    # merge fasta seq into mpmat table 
    # ------------------------------------------------------------------------>>>>>>> 
    # add seq info into pd.df
    input_fasta = SeqIO.parse(fasta_filename, format="fasta")

    region_seq_list = []
    for record in input_fasta:
        region_seq_list.append(str(record.seq).upper())

    input_fasta.close()
        
    # fix region_bed_table
    region_bed_table["extend_chr_start"] +=1
    region_bed_table["seq"] = region_seq_list

    # add all info and make a merge table 
    mpmat_table_merge = pd.concat( [ mpmat_table, region_bed_table[["extend_chr_start", "extend_chr_end", "seq"]]], axis=1)

    # remove FASTA 
    os.remove(fasta_filename)

    # ------------------------------------------------------------------------>>>>>>>
    # make alinger
    # ------------------------------------------------------------------------>>>>>>> 
    align_match, align_mismatch, align_gap_open, align_gap_extension = map(int,str(ARGS.align_settings).split(","))

    sgRNA_aligner = Bio.Align.PairwiseAligner()
    sgRNA_aligner.match = align_match
    sgRNA_aligner.mismatch = align_mismatch
    sgRNA_aligner.open_gap_score = align_gap_open
    sgRNA_aligner.extend_gap_score = align_gap_extension
    sgRNA_aligner.query_left_gap_score = align_gap_open
    sgRNA_aligner.query_right_gap_score = 0
    sgRNA_aligner.mode = "global"

    no_PAM_sgRNA_aligner = Bio.Align.PairwiseAligner()
    no_PAM_sgRNA_aligner.match = align_match
    no_PAM_sgRNA_aligner.mismatch = align_mismatch
    no_PAM_sgRNA_aligner.open_gap_score = align_gap_open
    no_PAM_sgRNA_aligner.extend_gap_score = align_gap_extension
    no_PAM_sgRNA_aligner.query_left_gap_score = 0
    no_PAM_sgRNA_aligner.query_right_gap_score = 0
    no_PAM_sgRNA_aligner.mode = "global"

    # Log
    sys.stderr.write("-" * 80 + "\n")
    sys.stderr.write("sgRNA Aligner:\n%s" % str(sgRNA_aligner))

    sys.stderr.write("-" * 80 + "\n")

    sys.stderr.write("No sgRNA Aligner:\n%s" % str(no_PAM_sgRNA_aligner))
    sys.stderr.write("-" * 80 + "\n")

    # ------------------------------------------------------------------------>>>>>>>
    # open output file 
    # ------------------------------------------------------------------------>>>>>>> 
    if ARGS.out_alignment_result == "Stdout":
        output_file = sys.stdout
    else:
        output_file = open(ARGS.out_alignment_result,"w")

    # ------------------------------------------------------------------------>>>>>>>
    # set output colname
    # ------------------------------------------------------------------------>>>>>>> 
    add_colname_list = None
    add_colname_state = False
    if ARGS.more_colname != "None":
        add_colname_list = ARGS.more_colname.split(",")
        add_colname_state = all([ x in mpmat_table_merge.columns for x in add_colname_list])

    # ------------------------------------------------------------------------>>>>>>>
    # output header
    # ------------------------------------------------------------------------>>>>>>> 
    out_header_list = [
        "chr_name",
        "region_start",
        "region_end",
        "region_id",
        "region_site_index",
        "SNP_ann",
        "Pass_info",
        "align_direction",
        "PAM_type",
        "PAM_seq",
        "align_total_match",
        "align_total_mismatch",
        "align_total_gap",
        "align_score",
        "seed_5bp_mismatch",
        "seed_5bp_gap",
        "head_5bp_mismatch",
        "head_5bp_gap",
        "head_8bp_mismatch",
        "head_8bp_gap",
        "align_target_seq",
        "align_info_state",
        "align_query_seq",
        "align_chr_name",
        "align_chr_start",
        "align_chr_end"
    ]

    if add_colname_state:
        out_header_list += add_colname_list

    output_file.write("\t".join(out_header_list) + "\n")
    # ------------------------------------------------------------------------>>>>>>>
    # run alignment
    # ------------------------------------------------------------------------>>>>>>> 
    DNA_rev_cmp_dict = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "-":"-" }
    align_min_score = int(ARGS.align_min_score)
    mut_region_tolerance_score = 5

    sgRNA_align_res_list = []

    for index in mpmat_table_merge.index:
        
        # ------------------------------->>>>>>
        # prepare ref_seq to align
        # ------------------------------->>>>>>
        ref_seq = mpmat_table_merge["seq"][index]
        ref_seq_BioPy = Bio.Seq.Seq(ref_seq, Bio.Alphabet.IUPAC.unambiguous_dna)
        ref_seq_rc = str(ref_seq_BioPy.reverse_complement())
        mut_region_type = get_mut_region_strand_direction( mpmat_table_merge["region_site_index"][index] )
        
        # ------------------------------->>>>>>
        # PAM alignment
        # ------------------------------->>>>>>    
        # PAM fwd alignment 
        final_align_fwd = run_sgRNA_alignment(ref_seq, sgRNA_seq, sgRNA_aligner, 3)
        
        # PAM rev alignment 
        final_align_rev = run_sgRNA_alignment(ref_seq_rc, sgRNA_seq, sgRNA_aligner, 3)
        
        # get fwd and rev state 
        align_state_fwd = False
        align_state_rev = False
        final_align_direction = None
        
        if len(final_align_fwd) > 0:
            if final_align_fwd[0][3] >= align_min_score:
                align_state_fwd = True
                
        if len(final_align_rev) > 0:
            if final_align_rev[0][3] >= align_min_score:
                align_state_rev = True            
        
        if align_state_fwd and align_state_rev:
            if mut_region_type == "CT":
                if final_align_fwd[0][3] + mut_region_tolerance_score >= final_align_rev[0][3]:
                    final_align_direction = "PAM_fwd"
                else:
                    final_align_direction = "PAM_rev"
            
            elif mut_region_type == "GA":
                if final_align_rev[0][3] + mut_region_tolerance_score >= final_align_fwd[0][3]:
                    final_align_direction = "PAM_rev"
                else:
                    final_align_direction = "PAM_fwd"
                
        elif align_state_fwd and (not align_state_rev):
            final_align_direction = "PAM_fwd"
        
        elif (not align_state_fwd) and (align_state_rev):
            final_align_direction = "PAM_rev"
        
        else:
            final_align_direction = "No_PAM"
        
        # ------------------------------->>>>>>
        # run no PAM alignment 
        # ------------------------------->>>>>>   
        if final_align_direction == "No_PAM":
            # fwd 
            final_align_fwd_no_PAM = run_no_PAM_sgRNA_alignment_no_chop(ref_seq, sgRNA_seq_full, no_PAM_sgRNA_aligner)
            
            # rev 
            final_align_rev_no_PAM = run_no_PAM_sgRNA_alignment_no_chop(ref_seq_rc, sgRNA_seq_full, no_PAM_sgRNA_aligner )
            
            if final_align_fwd_no_PAM[0][3] > final_align_rev_no_PAM[0][3]:
                final_align_direction = "No_PAM_fwd"
                
            else:
                final_align_direction = "No_PAM_rev"
            
        # ------------------------------->>>>>>
        # get final alignment 
        # ------------------------------->>>>>>
        if final_align_direction == "PAM_fwd":
            final_align = final_align_fwd[0]
        
        elif final_align_direction == "PAM_rev":
            final_align = final_align_rev[0]
        
        elif final_align_direction == "No_PAM_fwd":
            final_align = final_align_fwd_no_PAM[0]
        
        elif final_align_direction == "No_PAM_rev":
            final_align = final_align_rev_no_PAM[0]
        
        # ------------------------------->>>>>>
        # analysis and make table 
        # ------------------------------->>>>>>
        # get alignment chr position
        ref_seq_chr_name, ref_seq_chr_start, ref_seq_chr_end = region_bed_table.loc[index, ["extend_chr_name","extend_chr_start","extend_chr_end"]]
        align_pos_chr_res = align_position_on_chrom(final_align, final_align_direction, ref_seq_chr_name, ref_seq_chr_start, ref_seq_chr_end, len(ref_seq))    

        # get alignment stats
        align_res_stats = back_align_res_stats(align_direction=final_align_direction, align_res_list=final_align, align_ref_seq=ref_seq, align_sgRNA_seq_full=sgRNA_seq_full)

        # make output list 
        out_list = []
        out_list += list(mpmat_table_merge.loc[index,["chr_name", "region_start", "region_end", "region_id","region_site_index","SNP_ann","Pass_info"]])
        out_list += align_res_stats
        out_list += align_pos_chr_res

        # add column
        if add_colname_state:
            out_list += list(mpmat_table_merge.loc[index,add_colname_list])

        output_file.write("\t".join(map(str,out_list)) + "\n")


    # close file
    output_file.close()

# 2019-12-3 by MENG Howard
