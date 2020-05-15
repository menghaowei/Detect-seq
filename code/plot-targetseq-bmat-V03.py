#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
  2019-11-13 
        TargetSeq bmat file plot

Version-02:
  2019-12-01 
        1. Add sgRNA sequence
        2. Define plot region 
        3. Re-write sgRNA alignment methods

Version-03:
  2020-03-22
        1. fix get rev dict bug

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

1. try to find all PAM and align sgRNA related to the PAM position

2. if no appropriate alignment, try to run alignment with reverse complementary seq

3. get the plot region

4. make plot paramters and data

5. plot

"""
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import os 
import string
import time

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')

from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection

import Bio
from Bio import Align
from Bio import SeqIO

np.set_printoptions(suppress=True)


###############################################################################
# function part 
###############################################################################
def hex_to_rgb(hex_color):
    """
    INPUT
        <hex_color> 
            Color format  like #FFFFAA
    
    RETURN
        <rgb_color>    
            RGB tuple like (255, 255, 170)
    """
    value = hex_color.lstrip('#')
    value_len = len(value)
    return( tuple(int(value[index: index + 2], 16) for index in range(0, value_len, 2)) )


def rgb_to_hex(rgb_color):
    """
    INPUT
        <rgb_color>
            Color format like (255, 255, 170)
    
    RETURN
        <hex_color>
            Color format  like #FFFFAA
    """
    rgb_color = ('#%02x%02x%02x' % rgb_color).upper()
    return(rgb_color)


def make_color_list(low_color_RGB, high_color_RGB, length_out = 20, back_format="Hex"):
    """
    INPUT
        <low_color_RGB> <high_color_RGB>
            Format like (210, 179, 150), tuple, list, or np.array
        
        <back_format>
            Hex OR RGB
        
    RETURN
        <color_list>
    """
    low_color = np.array(low_color_RGB)
    high_color = np.array(high_color_RGB)
    
    color_list = []
    for index in range(0, length_out + 1):
        rgb_color = low_color + (high_color - low_color) // length_out * index
        if back_format == "Hex":
            color_list.append(rgb_to_hex(tuple(rgb_color)))
        else:
            color_list.append(tuple(rgb_color))
            
    return(color_list)


def map_color(value_vec, breaks, color_list):
    """
    INPUT:
        <value_vec>
            np.array or a list of values
            
        <breaks>
            A sorted value list, which can split all num into len(color_list) intervals. 
            e.g. [0.01, 0.1, 0.5, 1] make all real num into 5 intervals, (-Inf,0.01], (0.01,0.1], (0.1, 0.5],  (0.5, 1], (1, +Inf] 
        
        <color_list>
            A hex-format color list, which have to match with breaks
    
    RETURN
        <value_color_vec>
            A list map the value_vec with breaks 
    """
    value_idx_list = []
    
    for value in value_vec:
        match_state = False
        for index, break_value in enumerate(breaks):
            if value <= break_value:
                value_idx_list.append(index)
                match_state = True
                break
        
        if not match_state:
            value_idx_list.append(index+1)
    
    return( tuple(color_list[col_idx] for col_idx in value_idx_list))


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
        PAM_type = ref_seq[PAM_start_index : PAM_start_index+3]

        align_analysis_res += [PAM_start_index, PAM_type]
        align_res_list.append(align_analysis_res)

    if len(align_res_list) == 1:
        return (align_res_list)

    else:
        align_res_list.sort(cmp=cmp_align_list)
        return (align_res_list)    


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

    parser = argparse.ArgumentParser(description="convert mpileup file to info file")

    parser.add_argument("-i", "--input_bmat",
        help=".bmat file of Target-Seq data",required=True)

    parser.add_argument("--sgRNA",
        help="sgRNA sequence with PAM (NGG/NAG) sequence",default=None)

    parser.add_argument("-o","--out_figure",
        help="Output figure filename",required=True)

    parser.add_argument("--out_figure_format",
        help="Support 'pdf' and 'png' result, default=pdf",default="pdf")

    parser.add_argument("--out_figure_dpi",
        help="Out figure dpi, default=100",default="100")

    parser.add_argument("--show_indel",
        help="If show indel info in the out figure, default=True",default="True")

    parser.add_argument("--show_index",
        help="If show index info in the out figure, default=True",default="True")

    parser.add_argument("--box_border",
        help="If show box border in the out figure, default=False",default="False")

    parser.add_argument("--box_space",
        help="Set space size between two boxes, default=0.03",default="0.03")

    parser.add_argument("--min_color",
        help="min color to plot heatmap with RGB format, default=250,239,230",default="250,239,230")

    parser.add_argument("--max_color",
        help="max color to plot heatmap with RGB format, default=154,104,57",default="166,117,71")

    parser.add_argument("--min_ratio",
        help="Lower than this ratio plot as white color, default=0.001",default="0.001")

    parser.add_argument("--max_ratio",
        help="Higher than this ratio plot as max color, default=0.99",default="0.99")    

    parser.add_argument("--region_extend_length",
        help="From the middle site to extend <region_extend_length> bp at both side, default=25",default="25")

    parser.add_argument("--align_settings",
        help="Set <align_match_score>  <align_mismatch_score> <align_gap_open_score> <align_gap_extension_score>, default=5,-4,-24,-8",default="5,-4,-24,-8")

    parser.add_argument("--align_min_score",
        help="If alignment score lower than this, consider as no appropriate alignment, default=15",default="15")

    ARGS = parser.parse_args() 

    ###############################################################################
    # read parameters and load file 
    ###############################################################################
    # set sgRNA info
    if ARGS.sgRNA in default_sgRNA_dict.keys():
        sgRNA_full_length = default_sgRNA_dict[ARGS.sgRNA]
    else:
        sgRNA_full_length = str(ARGS.sgRNA).upper()
    
    sgRNA_seq = sgRNA_full_length[:-3]
    sgRNA_PAM_info = sgRNA_full_length[-3:]

    # set plot region
    region_extend_length = int(ARGS.region_extend_length)

    # set min align score
    align_min_score = float(ARGS.align_min_score)

    # ---------------------------------------------------------------->>>>>
    # load .bmat file
    # ---------------------------------------------------------------->>>>>
    bmat_table = pd.read_csv(ARGS.input_bmat,sep="\t")
    ref_seq =  "".join(bmat_table.ref_base)

    # ---------------------------------------------------------------->>>>>
    # set alignment 
    # ---------------------------------------------------------------->>>>>
    align_match, align_mismatch, align_gap_open, align_gap_extension = map(int,str(ARGS.align_settings).split(","))

    sgRNA_aligner = Align.PairwiseAligner()
    sgRNA_aligner.match = align_match
    sgRNA_aligner.mismatch = align_mismatch
    sgRNA_aligner.open_gap_score = align_gap_open
    sgRNA_aligner.extend_gap_score = align_gap_extension
    sgRNA_aligner.query_left_gap_score = align_gap_open
    sgRNA_aligner.query_right_gap_score = 0
    sgRNA_aligner.mode = "global"

    print("-" * 80)
    print(str(sgRNA_aligner))
    print("-" * 80)

    no_PAM_sgRNA_aligner = Align.PairwiseAligner()
    no_PAM_sgRNA_aligner.match = align_match
    no_PAM_sgRNA_aligner.mismatch = align_mismatch
    no_PAM_sgRNA_aligner.open_gap_score = align_gap_open
    no_PAM_sgRNA_aligner.extend_gap_score = align_gap_extension
    no_PAM_sgRNA_aligner.query_left_gap_score = 0
    no_PAM_sgRNA_aligner.query_right_gap_score = 0
    no_PAM_sgRNA_aligner.mode = "global"

    print("-" * 80)
    print(str(no_PAM_sgRNA_aligner))
    print("-" * 80)

    # ---------------------------------------------------------------->>>>>
    # alignment 
    # ---------------------------------------------------------------->>>>>
    sgRNA_align_extend_len = 3

    # make ref seq
    ref_seq_BioPy = Bio.Seq.Seq(ref_seq, Bio.Alphabet.IUPAC.unambiguous_dna)
    ref_seq_rc = str(ref_seq_BioPy.reverse_complement())

    # PAM fwd alignment 
    final_align_fwd = run_sgRNA_alignment(ref_seq, sgRNA_seq, sgRNA_aligner, sgRNA_align_extend_len)

    # PAM rev alignment 
    final_align_rev = run_sgRNA_alignment(ref_seq_rc, sgRNA_seq, sgRNA_aligner, sgRNA_align_extend_len)

     # fwd alignment 
    print("Forward best alignment:")
    print(final_align_fwd[0][6][final_align_fwd[0][4] : final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][7][final_align_fwd[0][4] : final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][8][final_align_fwd[0][4] : final_align_fwd[0][5] + 1])
    print(final_align_fwd[0][3])
    print("-" * 80)

    # rev alignment 
    print("Reverse best alignment:")
    print(final_align_rev[0][6][final_align_rev[0][4] : final_align_rev[0][5] + 1])
    print(final_align_rev[0][7][final_align_rev[0][4] : final_align_rev[0][5] + 1])
    print(final_align_rev[0][8][final_align_rev[0][4] : final_align_rev[0][5] + 1])
    print(final_align_rev[0][3])
    print("-" * 80)

    # make alignment info 
    sgRNA_align = [""] * len(ref_seq)
    sgRNA_align_insert = [""] * len(ref_seq)

    final_align_direction = None
    DNA_rev_cmp_dict = { "A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "-":"-" }

    # define align direction and final align res
    if final_align_fwd[0][3] >= final_align_rev[0][3]:
        if final_align_fwd[0][3] >= align_min_score:
            final_align_direction = "Forward PAM Alignment"
            final_align = final_align_fwd[0]
            
    elif final_align_rev[0][3] > final_align_fwd[0][3]:
        if final_align_rev[0][3] >= align_min_score:
            final_align_direction = "Reverse PAM Alignment"
            final_align = final_align_rev[0]

    if final_align_direction == None:
        final_align_direction = "No PAM Alignment"
        final_align = run_no_PAM_sgRNA_alignment_no_chop(ref_seq, sgRNA_full_length, no_PAM_sgRNA_aligner)[0]

    # make sgRNA alignment info 
    final_align_info  = final_align[7][final_align[4] : final_align[5] + 1]
    final_align_ref = final_align[6][final_align[4] : final_align[5] + 1]
    final_align_sgRNA = final_align[8][final_align[4] : final_align[5] + 1]

    final_align_ref_gap_count = final_align_ref.count("-")
    final_align_sgRNA_gap_count = final_align_sgRNA.count("-")

    ref_align_gap_count = 0 
    ref_del_str = ""

    if (final_align_direction == "Forward PAM Alignment") or (final_align_direction == "Reverse PAM Alignment") :
        sgRNA_start = final_align[9] - len(sgRNA_seq) - final_align_sgRNA_gap_count + final_align_ref_gap_count
        
        for align_index, align_ref in enumerate(final_align_ref):
            if align_ref != "-":
                sgRNA_align[sgRNA_start + align_index - ref_align_gap_count] = ref_del_str + final_align_sgRNA[align_index]
                ref_del_str = ""
            else:
                ref_align_gap_count += 1 
                ref_del_str += final_align_sgRNA[align_index]            
                sgRNA_align_insert[sgRNA_start + align_index] = final_align_sgRNA[align_index]
        
        # add PAM info 
        PAM_ref_start_index = sgRNA_start + align_index + 1 - ref_align_gap_count
        sgRNA_align[PAM_ref_start_index : PAM_ref_start_index + 3] = sgRNA_full_length[-3:]

        if final_align_direction == "Reverse PAM Alignment":
            sgRNA_align = sgRNA_align[::-1]
            sgRNA_align_insert = sgRNA_align_insert[::-1]

    elif final_align_direction == "No PAM Alignment":
        sgRNA_start = final_align[4] 
        
        for align_index, align_ref in enumerate(final_align_ref):
            if align_ref != "-":
                sgRNA_align[sgRNA_start + align_index - ref_align_gap_count] = final_align_sgRNA[align_index]
            else:
                ref_align_gap_count += 1 
                sgRNA_align_insert[sgRNA_start + align_index] = final_align_sgRNA[align_index]

    # set possible_sgRNA_region
    if final_align_direction == None:
        possible_sgRNA_region = [len(ref_seq) // 2  - region_extend_length , len(ref_seq) // 2  + region_extend_length]

    else:
        if final_align_direction == "No PAM Alignment": 
            sgRNA_align_start = final_align[4] 
            sgRNA_align_end = final_align[5]

        else:
            if final_align_direction == "Forward PAM Alignment":
                print 1
                sgRNA_align_end = final_align[9] + 3 
                sgRNA_align_start = sgRNA_align_end - len(sgRNA_full_length)        
            
            elif final_align_direction == "Reverse PAM Alignment":
                print 2
                sgRNA_align_start = len(ref_seq) - (final_align[9] + 3 )
                sgRNA_align_end = sgRNA_align_start + len(sgRNA_full_length)
        
        possible_sgRNA_region_start = max(sgRNA_align_start - region_extend_length, 0)
        possible_sgRNA_region_end = min(sgRNA_align_end + region_extend_length, len(ref_seq) - 1 )
        possible_sgRNA_region = [possible_sgRNA_region_start, possible_sgRNA_region_end]
        
    # select bmat_table 
    bmat_table_select = bmat_table[possible_sgRNA_region[0] : possible_sgRNA_region[1]]   

    # ---------------------------------------------------------------->>>>>
    # make plot 
    # ---------------------------------------------------------------->>>>>

    # --------------------------------------------------->>>>>
    # set color
    # --------------------------------------------------->>>>>
    # show indel
    indel_plot_state = eval(ARGS.show_indel)
    index_plot_state = eval(ARGS.show_index)
    box_border_plot_state = eval(ARGS.box_border)

    # set panel size 
    panel_box_width = 0.4
    panel_box_heigth = 0.4
    panel_space = 0.05
    panel_box_space = float(ARGS.box_space)

    # color part 
    base_color_dict = {"A":"#04E3E3", "T":"#F9B874", "C":"#B9E76B", "G":"#F53798","N":"#AAAAAA","-":"#AAAAAA"}

    # make color breaks 
    color_break_num = 20
    break_step = 1.0 / color_break_num
    min_color_value = float(ARGS.min_ratio)
    max_color_value = float(ARGS.max_ratio)
    color_break = np.round(np.arange(min_color_value,max_color_value,break_step),5)

    # make color list 
    # low_color = (210, 179, 150)
    # low_color = (240, 225, 212)
    # low_color = (250, 239, 230)
    # high_color = (154, 104, 57)
    # low_color = (217, 231, 245)
    # high_color = (9, 42, 96)

    low_color = tuple( map(int, ARGS.min_color.split(",")) )
    high_color = tuple( map(int, ARGS.max_color.split(",")) )

    try:
        color_list = make_color_list(low_color, high_color, len(color_break)-1,"Hex")
        color_list = ["#FFFFFF"] + color_list
    except:
        print low_color, high_color
        print color_break

    # get plot info 
    total_box_count = possible_sgRNA_region[1] - possible_sgRNA_region[0]

    # calculate base info and fix zero
    base_sum_count = bmat_table_select[["A","G","C","T"]].apply(lambda x: x.sum(), axis=1)
    total_sum_count = bmat_table_select[["A","G","C","T","del_count","insert_count"]].apply(lambda x: x.sum(), axis=1)
    base_sum_count[np.where(base_sum_count==0)[0]] = 1
    total_sum_count[np.where(total_sum_count==0)[0]] = 1

    # make plot size 
    plot_data_list = None
    panel_space_coef = None
    panel_height_coef = None

    if indel_plot_state:
        panel_height_coef = [0.5, 0.9, 0.9] + [0.5] * 6 + [0.5] * 6
        panel_space_coef = [1,1,1] + [0.3] * 3 + [1,0.3,1] + [0.3] * 3 + [1,0.3]
        plot_data_list = [
            ["Index",bmat_table_select.chr_index],
            ["On-target",sgRNA_align[possible_sgRNA_region[0] : possible_sgRNA_region[1]]],
            ["Ref",bmat_table_select.ref_base],
            ["A",np.array(bmat_table_select["A"])],
            ["G",np.array(bmat_table_select["G"])],
            ["C",np.array(bmat_table_select["C"])],
            ["T",np.array(bmat_table_select["T"])], 
            ["Del",np.array(bmat_table_select["del_count"])],
            ["Ins",np.array(bmat_table_select["insert_count"])],
            ["A.ratio",list(bmat_table_select["A"] / base_sum_count )],
            ["G.ratio",list(bmat_table_select["G"] / base_sum_count )],
            ["C.ratio",list(bmat_table_select["C"] / base_sum_count )],
            ["T.ratio",list(bmat_table_select["T"] / base_sum_count )], 
            ["Del.ratio",list(bmat_table_select["del_count"] / total_sum_count )],
            ["Ins.ratio",list(bmat_table_select["insert_count"] / total_sum_count)]
        ]
    else:
        panel_height_coef = [0.5, 0.9, 0.9] + [0.5] * 4 + [0.5] * 4
        panel_space_coef = [1,1,1] + [0.3] * 3 + [1] + [0.3] * 3
        plot_data_list = [
            ["Index",bmat_table_select.chr_index],
            ["On-target",sgRNA_align[possible_sgRNA_region[0] : possible_sgRNA_region[1]]],
            ["Ref",bmat_table_select.ref_base],
            ["A",np.array(bmat_table_select["A"])],
            ["G",np.array(bmat_table_select["G"])],
            ["C",np.array(bmat_table_select["C"])],
            ["T",np.array(bmat_table_select["T"])], 
            ["A.ratio",list(bmat_table_select["A"] / base_sum_count )],
            ["G.ratio",list(bmat_table_select["G"] / base_sum_count )],
            ["C.ratio",list(bmat_table_select["C"] / base_sum_count )],
            ["T.ratio",list(bmat_table_select["T"] / base_sum_count )]
        ]

    # get box and space info
    box_height_list = np.array(panel_height_coef) * panel_box_heigth
    panel_space_list = np.array(panel_space_coef) * panel_space
        
    # calculate figure total width and height
    figure_width = total_box_count * panel_box_width  + (total_box_count - 1) * panel_box_space + panel_box_width * 2
    figure_height = sum(box_height_list) + sum(panel_space_list)

    # make all box_x 
    box_x_vec = np.arange(0, figure_width+panel_box_width, panel_box_width + panel_box_space)
    box_x_vec = box_x_vec[:(len(ref_seq) + 1)]

    # make box border 
    if box_border_plot_state:
        box_edgecolor = "#AAAAAA"
        box_linestyle = "-"
        box_linewidth = 2
    else:
        box_edgecolor = "#FFFFFF"
        box_linestyle = "None"
        box_linewidth = 0

    # make box_y initialize
    current_y = 0

    # ---------------------------------------------------------------->>>>>>>>
    # plot region
    # ---------------------------------------------------------------->>>>>>>>
    # set new figure 
    fig = plt.figure(figsize=(figure_width*1.1, figure_height*1.1))
    ax = fig.add_subplot(111,aspect="equal")
    plt.xlim([0,figure_width])
    plt.ylim([-figure_height,0])
    plt.axis("off")

    # make plot 
    text_list = []
    patches = []

    for panel_index in range(len(panel_height_coef)):
        
        # panel name 
        panel_name = plot_data_list[panel_index][0]
        panel_name_x = box_x_vec[0]
        panel_name_y = current_y - box_height_list[panel_index] * 0.5
        text_list.append((panel_name_x, panel_name_y, panel_name, 10))
        
        # plot panel box 
        if panel_name == "Index":      
            # don't draw box, only add text
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                box_x = box_x_vec[index+1]
                text_list.append((box_x + panel_box_width * 0.5, current_y - box_height_list[panel_index] * 0.5, str(box_value), 10))
               
            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])
            
        elif panel_name in ["On-target", "Ref"]:
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                if box_value == "":
                    box_fill = False
                    box_color = "#FFFFFF"

                else:
                    if "Reverse" in final_align_direction:
                        if panel_name == "Ref":
                            box_value =  "".join([DNA_rev_cmp_dict.get(x) for x in box_value])
                        else:
                            pass
                        
                        box_color = base_color_dict.get(box_value[0])
                        
                    else:
                        box_color = base_color_dict.get(box_value[-1])
                        
                    if not box_color:
                        box_fill = False
                        box_color = "#FFFFFF"
                    else:
                        box_fill = True

                box_x = box_x_vec[index+1]

                patches.append(Rectangle( 
                    xy = (box_x, current_y - box_height_list[panel_index]), 
                    width = panel_box_width,  
                    height = box_height_list[panel_index], 
                    fill = box_fill,
                    alpha = 1,
                    linestyle = box_linestyle,
                    linewidth = box_linewidth,
                    edgecolor = box_edgecolor,
                    facecolor= box_color)
                )

                # text 
                text_list.append((box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index] , str(box_value), 16))

            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])   

        elif (panel_name in ["A","G","C","T","Del","Ins"]):
            if panel_name in ["Del","Ins"]:
                box_ratio = plot_data_list[panel_index][1] / total_sum_count
            else:
                box_ratio = plot_data_list[panel_index][1] / base_sum_count
                
            box_color_list = map_color(box_ratio, color_break, color_list)
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                box_color = box_color_list[index]
                box_x = box_x_vec[index+1]
                
                patches.append(Rectangle( 
                    xy = (box_x, current_y - box_height_list[panel_index]), 
                    width = panel_box_width,  
                    height = box_height_list[panel_index], 
                    fill = True,
                    alpha = 1,
                    linestyle = box_linestyle,
                    linewidth = box_linewidth,
                    edgecolor = box_edgecolor,
                    facecolor = box_color)
                ) 
            
                # text 
                text_list.append((box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index], str(box_value), 6))
                
            # make next panel_y
            current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])        

        else:
            box_color_list = map_color(plot_data_list[panel_index][1], color_break, color_list)
            for index, box_value in enumerate(plot_data_list[panel_index][1]):
                box_color = box_color_list[index]
                box_x = box_x_vec[index+1]

                patches.append(Rectangle( 
                    xy = (box_x, current_y-box_height_list[panel_index]), 
                    width = panel_box_width,  
                    height = box_height_list[panel_index], 
                    fill = True,
                    alpha = 1,
                    linestyle = box_linestyle,
                    linewidth = box_linewidth,
                    edgecolor = box_edgecolor,
                    facecolor = box_color)
                ) 

                # text 
                text_list.append((box_x + 0.5 * panel_box_width, current_y - 0.5 * box_height_list[panel_index], round(box_value*100,4), 6))

            # make next panel_y
            if panel_index < len(panel_space_list):
                current_y = current_y - (box_height_list[panel_index] + panel_space_list[panel_index])        
               
    # plot box 
    ax.add_collection(PatchCollection(patches, match_original=True))

    # add text 
    for text_x, text_y, text_info, text_fontsize in text_list:
         plt.text(
            x = text_x,
            y = text_y,
            s = text_info,
            horizontalalignment='center',
            verticalalignment='center',
            fontsize = text_fontsize,
            fontname="Arial"
         )

    # output plot
    if str(ARGS.out_figure_format) == "PNG":
        fig.savefig(fname = ARGS.out_figure, bbox_inches='tight', dpi = int(ARGS.out_figure_dpi), format = "png")

    else:
        fig.savefig(fname = ARGS.out_figure, bbox_inches='tight', dpi = int(ARGS.out_figure_dpi), format = "pdf")


    









   





# 2019-12-01 by MENG Howard



