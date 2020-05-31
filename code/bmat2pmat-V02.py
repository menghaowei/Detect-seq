#! /Users/meng/menghw_HD/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
  2019-06-10 对生成的bmat格式进行转换，目标为类似bed格式的pmat

Version-02:
  2019-08-04 修复bug：
             1. 直接从bmat格式转成pmat格式，而不需要增加一列形成类似bed的format
             2. 更名bmatbed2pmat -> bmat2pmat

E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
output format like:

<chr_name> <chr_index> <site_index> <A> <G> <T> <C> <type> <ref_base> <mut_base> <ref_num> <mut_num> <cover_num> <mut_ratio>
chr1    54795   54795   chr1_54795_TC   1   0   2   6   TC  T   C   6   2   8   0.25

TAB separate
"""
# Learning Part END-----------------------------------------------------------

import argparse
import gzip
import sys

###############################################################################
# function part 
###############################################################################
def parse_line(Line,inlike_bed=False):
    """
    FUN:
        从bmat格式中解析数据

    RETURN:
        [ref_base, ref_num, mut_base, mut_num]
    """
    Line_list = Line.strip().split("\t")
    
    if inlike_bed:
        # set intial 
        ref_base = Line_list[3]
        ref_num = 0
        mut_base = None
        mut_num = 0
        
        # get ref_num and mut_num
        base_dict = {
            "A":int(Line_list[4]),
            "G":int(Line_list[5]),
            "C":int(Line_list[6]),
            "T":int(Line_list[7])
        }

    else:
        # set intial 
        ref_base = Line_list[2]
        ref_num = 0
        mut_base = None
        mut_num = 0
        
        # get ref_num and mut_num
        base_dict = {
            "A":int(Line_list[3]),
            "G":int(Line_list[4]),
            "C":int(Line_list[5]),
            "T":int(Line_list[6])
        }
    
    for base_key in base_dict.keys():        
        if ref_base == base_key:
            ref_num = base_dict[base_key]
        
        else:
            if base_dict[base_key] >  mut_num:
                mut_base = base_key
                mut_num =  base_dict[base_key]
    
    return([ref_base, ref_num, mut_base, mut_num])

###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="convert bmat file to pmat file")

    parser.add_argument("-i", "--Input",
        help="Input bmat file",required=True)

    parser.add_argument("-o", "--Output",
        help="Output BED format file",default="Stdout")

    parser.add_argument("-c", "--CoverNumCutoff",
        help="Site coverage number cutoff default=0",default="0")

    parser.add_argument("-m", "--MutNumCutoff",
        help="Site mutation number cutoff default=0",default="0")

    parser.add_argument("-r", "--MutRatioCutoff",
        help="Site mutation ratio cutoff default=0",default="0")

    parser.add_argument("-t", "--MutType",
        help="Select mutation type, ALL means no selection, can set like CT, default=ALL",default="ALL")

    parser.add_argument("--InHeader",
        help="If contain header line in input file, default=True",default="True")

    parser.add_argument("--InLikeBED",
        help="If input bmat file looks like bed file, default=False",default="False")

    parser.add_argument("--OutHeader",
        help="If contain header line in output file, default=True",default="True")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    INPUT_FILE_PATH = ARGS.Input
    OUTPUT_FILE_PATH = ARGS.Output

    IF_IN_HEADER = eval(ARGS.InHeader)
    IF_OUT_HEADER = eval(ARGS.OutHeader)
    IF_LIKE_BED = eval(ARGS.InLikeBED)

    BASE_MUT_NUM_CUTOFF = int(ARGS.MutNumCutoff)
    BASE_MUT_RATIO_CUTOFF = float(ARGS.MutRatioCutoff)
    BASE_COVER_NUM_CUTOFF = int(ARGS.CoverNumCutoff)

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open input file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ".gz" in INPUT_FILE_PATH:
        INPUT_FILE = gzip.open(INPUT_FILE_PATH,"r")
    else:
        INPUT_FILE = open(INPUT_FILE_PATH,"r")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open output file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if OUTPUT_FILE_PATH == "Stdout":
        OUTPUT_FILE = sys.stdout
    elif ".gz" in OUTPUT_FILE_PATH:
        OUTPUT_FILE = gzip.open(OUTPUT_FILE_PATH,"w")
    else:
        OUTPUT_FILE = open(OUTPUT_FILE_PATH,"w")


    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # main part
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # output header 
    if IF_OUT_HEADER:
        header = [
            "chr_name",
            "chr_index",
            "site_index",
            "A",
            "G",
            "C",
            "T",
            "mut_type",
            "ref_base",
            "mut_base",
            "ref_num",
            "mut_num",
            "cover_num",
            "mut_ratio"
        ]
        OUTPUT_FILE.write("\t".join(header) + "\n")

    # parse bmat file
    if IF_IN_HEADER:
        data = INPUT_FILE.readline()

    # run part
    for line in INPUT_FILE:
        line_list = line.strip().split("\t")
        ref_base, ref_num, mut_base, mut_num = parse_line(line, inlike_bed = IF_LIKE_BED)
        cover_num = ref_num + mut_num

        if cover_num > 0:
            mut_ratio = round(mut_num / 1.0 / cover_num, 5)
        else:
            mut_ratio = 0

        if mut_base == None:
            mut_base = "."
        
        if (cover_num >= BASE_COVER_NUM_CUTOFF) and (mut_num >= BASE_MUT_NUM_CUTOFF) and (mut_ratio >= BASE_MUT_RATIO_CUTOFF):

            # make type 
            mut_type = ref_base + mut_base    
            
            # site index 
            site_index = "%s_%s_%s" % (line_list[0], line_list[1], mut_type)
            
            if IF_LIKE_BED:
                # make output list 
                out_list = line_list[:3] + [site_index] + line_list[4:8] + [mut_type, ref_base, mut_base, ref_num, mut_num, cover_num, mut_ratio]
            else:
                out_list = line_list[:2] + [line_list[1],site_index] + line_list[3:7] + [mut_type, ref_base, mut_base, ref_num, mut_num, cover_num, mut_ratio]
            
            out_str = "\t".join(map(str,out_list))

            if ARGS.MutType == "ALL":
                OUTPUT_FILE.write(out_str + "\n")
            elif ARGS.MutType == mut_type:
                OUTPUT_FILE.write(out_str + "\n")


INPUT_FILE.close()
OUTPUT_FILE.close()









# 2019-06-10