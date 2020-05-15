#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-02:
  2019-09-11 tolerance option fix 

Version-01:
  2019-08-08 mpmat select with different cutoff

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
output format like:

# column 1~3 basic info
<chr_name> <region_start> <region_end> 

# annotation info
<region_site_num> <region_PASS_mut_site_num> <region_site_index.list> <to_base_num.list> <cover_base_num.list> <mut_ratio.list> <SNP_info.list> <tandem_info> <PASS_info>

TAB separate
"""
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import sys
import time
import string

###############################################################################
# function part 
###############################################################################


###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="mpmat file selection")

    parser.add_argument("-i", "--Input",
        help="Input mpmat file",required=True)

    parser.add_argument("-o", "--Output",
        help="Output mpmat select file",default="Stdout")

    parser.add_argument("-f", "--FromBase",
        help="Ref base, accept A,G,C,T default=C",default="C")

    parser.add_argument("-t", "--ToBase",
        help="Mut base, accept A,G,C,T default=T",default="T")

    parser.add_argument("-m", "--SiteMutNum",
        help="Site-cutoff, mutation reads number, default=1",default="1")

    parser.add_argument("-c", "--SiteCoverNum",
        help="Site-cutoff, total cover reads number, default=5",default="5")

    parser.add_argument("-r", "--SiteMutRatio",
        help="Site-cutoff, site mutation ratio, default=0.1",default="0.1")        

    parser.add_argument("--RegionPassNum",
        help="Region-cutoff region should contain no-less than Pass site number, default=3",default="3")

    parser.add_argument("--RegionToleranceNum",
        help="Region-cutoff region should tolerace site number within a reproted region, default=1, if set as False, don't consider this paramter",default="1")

    parser.add_argument("--OutHeader",
        help="If contain header line in output file, default=True",default="True")

    parser.add_argument("--InHeader",
        help="If contain header line in input file, default=True",default="True")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    input_file_path = ARGS.Input
    output_file_path = ARGS.Output

    from_base = ARGS.FromBase
    to_base = ARGS.ToBase

    site_cutoff_coverage = int(ARGS.SiteCoverNum)
    site_cutoff_mut_num = int(ARGS.SiteMutNum)
    site_cutoff_mut_ratio = float(ARGS.SiteMutRatio)

    if ARGS.RegionToleranceNum == "False":
        region_cutoff_tandem_not_pass_tolerance = 10000
    else:
        region_cutoff_tandem_not_pass_tolerance = int(ARGS.RegionToleranceNum)

    region_cutoff_pass_site_num = int(ARGS.RegionPassNum)

    if_output_header = eval(ARGS.OutHeader)
    if_input_header = eval(ARGS.InHeader)
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open input file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ".gz" == input_file_path[-3:]:
        in_mpmat_file = gzip.open(input_file_path,"r")
    else:
        in_mpmat_file = open(input_file_path,"r")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open output file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if output_file_path == "Stdout":
        output_file = sys.stdout
    elif ".gz" == output_file_path[-3:]:
        output_file = gzip.open(output_file_path,"w")
    else:
        output_file = open(output_file_path,"w")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # header output
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # output header 
    if if_output_header:
        header = [
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
            "pass_info"
        ]
        output_file.write("\t".join(header) + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # header output
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if if_input_header:
        header = in_mpmat_file.readline()

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # main part
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    for mpmat_line in in_mpmat_file:
        mpmat_line_list = mpmat_line.strip().split("\t")
        
        # initialize the line 
        region_pass_start_idx = None
        region_pass_end_idx = None
        region_pass_count = 0
        
        region_pass_state_list = []
        site_index_list = mpmat_line_list[6].split(",")
        mut_num_list = map(int,mpmat_line_list[7].split(","))
        cover_num_list = map(int,mpmat_line_list[8].split(","))
        mut_ratio_list = map(float,mpmat_line_list[9].split(","))
        SNP_ann_list = mpmat_line_list[10].split(",")
        tandem_state_list = mpmat_line_list[11].split(",")
        
        # test SNP, and annotate site cutoff
        for index in range(int(mpmat_line_list[3])):
            # SNP 是否有注释
            if SNP_ann_list[index] == "True":
                region_pass_state_list.append("Filter")
                continue
            
            if site_index_list[index][-2:] != (from_base + to_base):
                region_pass_state_list.append("Filter")
                continue
            
            if mut_num_list[index] < site_cutoff_mut_num:
                region_pass_state_list.append("Filter")
                continue

            if cover_num_list[index] < site_cutoff_coverage:
                region_pass_state_list.append("Filter")
                continue

            if mut_ratio_list[index] < site_cutoff_mut_ratio:
                region_pass_state_list.append("Filter")
                continue

            # count part
            region_pass_state_list.append("Pass")

            if region_pass_start_idx == None:
                region_pass_start_idx = index
                
            region_pass_end_idx = index

        # test region cutoff
        if (region_pass_end_idx != None) and (region_pass_start_idx != None):
            select_region_site_count = region_pass_end_idx - region_pass_start_idx + 1
            select_region_pass_count = 0 
            select_region_filter_count = 0 
            select_region_SNP_count = 0

            for index in range(region_pass_start_idx, region_pass_end_idx+1):
                if region_pass_state_list[index] == "Pass":
                    select_region_pass_count += 1
                else:
                    select_region_filter_count += 1

                if SNP_ann_list[index] == "True":
                    select_region_SNP_count += 1

            if (select_region_pass_count >= region_cutoff_pass_site_num) and (select_region_filter_count <= region_cutoff_tandem_not_pass_tolerance):

                region_new_start = site_index_list[region_pass_start_idx].split("_")[1]
                region_new_end = site_index_list[region_pass_end_idx].split("_")[1]
                region_site_count = region_pass_end_idx - region_pass_start_idx + 1

                out_list = [
                    mpmat_line_list[0],
                    region_new_start,
                    region_new_end,
                    select_region_site_count,
                    select_region_pass_count,
                    select_region_SNP_count,
                    ",".join(site_index_list[region_pass_start_idx : (region_pass_end_idx + 1)]),
                    ",".join(map(str,mut_num_list[region_pass_start_idx : (region_pass_end_idx + 1)])),
                    ",".join(map(str,cover_num_list[region_pass_start_idx : (region_pass_end_idx + 1)])),
                    ",".join(map(str,mut_ratio_list[region_pass_start_idx : (region_pass_end_idx + 1)])),
                    ",".join(map(str,SNP_ann_list[region_pass_start_idx : (region_pass_end_idx + 1)])),
                    ",".join(tandem_state_list[region_pass_start_idx : (region_pass_end_idx + 1)]),
                    ",".join(region_pass_state_list[region_pass_start_idx : (region_pass_end_idx + 1)])
                ]

                out_str = "\t".join(map(str,out_list))
                output_file.write(out_str + "\n")

    

in_mpmat_file.close()
output_file.close()








# 2019-08-08