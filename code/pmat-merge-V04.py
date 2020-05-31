#! /Users/meng/menghw_HD/anaconda2/bin/python
# _*_ coding: UTF-8 _*_


# Version information START --------------------------------------------------
VERSION_INFO = \
"""
Author: MENG Howard

Version-01:
  2019-08-06 对bmat2pmat的输出结果进行merge

Version-02:
  2019-08-08 修改输出的列及header
             修改 back_trandem_site_index 函数可能会造成ref index溢出的问题

Version-03:
  2019-09-12 增加SNP注释信息的支持格式      

Version-04:
  2019-09-19 

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
<region_site_num> <region_conatin_mut_site_num(rm SNP)> <region_SNP_num> <region_site_index.list> <to_base_num.list> <cover_base_num.list> <mut_ratio.list> <SNP_info.list> <tandem_info>

TAB separate
"""
# Learning Part END-----------------------------------------------------------

from Bio import SeqIO
import argparse
import gzip
import sys
import time
import string

###############################################################################
# function part 
###############################################################################
def back_trandem_site_index(chr_name, chr_start, base_type, ref_dict, distance_size_cutoff = 150, back_site_num_cutoff = None):
    """
    <INPUT> 
        chr_name: 
            chr1
        
        chr_start: 
            1-based index, same with BED and VCF files
        
        base_type: 
            the site base info, accept A,T,C,G
        
        ref_dict:
            genome dict 
        
        distance_size_cutoff: 
            if region larger than <distance_size_cutoff> bp, the function with return list
        
        back_site_num_cutoff:
            if return list number 
        
    <RETURN>
        a list of site index, format like: 
        ["chr1_10100_C",...]
    
    """
    site_index_list = []
    chr_index = int(chr_start)

    if ref_dict[chr_name][chr_index -1] == base_type:
        site_index = "%s_%s_%s" % (chr_name, chr_index, base_type)
        site_index_list.append(site_index)
    else:
        return(None)

    chr_end_index = chr_index - 1 + distance_size_cutoff
    if chr_end_index > len(ref_dict[chr_name]):
        chr_end_index = len(ref_dict[chr_name])

    for run_chr_index in xrange(chr_index, chr_end_index):
        run_base = ref_dict[chr_name][run_chr_index]
        if run_base == base_type:
            site_index = "%s_%s_%s" % (chr_name, run_chr_index + 1, base_type)
            site_index_list.append(site_index)

            if back_site_num_cutoff != None:
                if len(site_index_list) >= back_site_num_cutoff:
                    return(site_index_list)

    return(site_index_list)


def load_VCF_as_dict(snp_vcf_file_obj, input_dict):
    """
    <HELP>
        load vcf file input a dict, and return dict
    """
    
    for line in snp_vcf_file_obj:
        if line[0] != "#":
            line_list = line.strip().split("\t")
            vcf_index = "%s_%s_%s%s" % (line_list[0], line_list[1], line_list[3], line_list[4])
            input_dict[vcf_index] = 0

    return(input_dict)


def load_BED_as_dict(snp_vcf_file_obj, input_dict):
    """
    <HELP>
        load vcf file input a dict, and return dict
    """
    
    for line in snp_vcf_file_obj:
        if line[0] != "#":
            line_list = line.strip().split("\t")
            site_index = line_list[3]
            input_dict[site_index] = 0

    return(input_dict)


def get_region_start_line(pmat_obj, from_base, to_base, pmat_temp_line = None):
    """
    <HELP>
        help to get region start site
    """
    # consider temp 
    if pmat_temp_line != None:
        line_list = pmat_temp_line.strip().split("\t")
        if (int(line_list[12]) > 0) and (line_list[9] == from_base) and (line_list[10] == to_base) :
            return(True,pmat_temp_line) 
    
    # initial the line info
    line = pmat_obj.readline()
    while line:
        line_list = line.strip().split("\t")
        if (int(line_list[12]) > 0) and (line_list[9] == from_base) and (line_list[10] == to_base) :
            return(True,line) 
        else:
            line = in_pmat_file.readline()
    
    return(False,line)    
  

def report_region(out_region_dict, to_base = "T"):
    """
    <HELP>
        make output str from region_dict
        report region from mut to mut
        
        header_order:
        # column 1~3 basic info
        <chr_name> <region_start> <region_end> 

        # annotation info
        <region_site_num> <region_conatin_mut_site_num(rm SNP)> <region_SNP_num> <region_site_index.list> <to_base_num.list> <cover_base_num.list> <mut_ratio.list> <SNP_info.list> <tandem_info>
    """
    out_str_list = []

    # site selection
    region_start_mut_index = None
    region_end_mut_index = None
    
    # get start index
    for index in  range(len(out_region_dict["mut_count_list"])):
        if (int(out_region_dict["mut_count_list"][index]) > 0) and (out_region_dict["site_idx_list"][index][-1] == to_base): 

            if not out_region_dict["SNP_ann_list"][index]:
                region_start_mut_index = index
                break

    # get end index
    for index in range(len(out_region_dict["mut_count_list"]) -1 , -1 , -1):
        if (int(out_region_dict["mut_count_list"][index]) > 0) and (out_region_dict["site_idx_list"][index][-1] == to_base): 
            
            if not out_region_dict["SNP_ann_list"][index]:
                region_end_mut_index = index
                break

    if (region_start_mut_index == None) or (region_end_mut_index == None):
        return(None)
    
    else:
        # column 1~3 basic info
        out_str_list.append(out_region_dict["chr_name"])
        out_str_list.append(out_region_dict["site_idx_list"][region_start_mut_index].split("_")[1])
        out_str_list.append(out_region_dict["site_idx_list"][region_end_mut_index].split("_")[1])

        # region site number 
        out_str_list.append(region_end_mut_index - region_start_mut_index + 1)
        
        # get mutation count
        mut_signal_site_num = 0 
        snp_site_num = 0
        
        for index in range(region_start_mut_index, region_end_mut_index + 1):
            if out_region_dict["SNP_ann_list"][index]:
                snp_site_num += 1
                continue
            
            if (int(out_region_dict["mut_count_list"][index]) > 0) and (out_region_dict["site_idx_list"][index][-1] == to_base):
                mut_signal_site_num += 1

        out_str_list.append(mut_signal_site_num)
        out_str_list.append(snp_site_num)

        # annotation info 
        out_str_list.append(",".join(map(str,out_region_dict["site_idx_list"][region_start_mut_index : region_end_mut_index+1])))
        out_str_list.append(",".join(map(str,out_region_dict["mut_count_list"][region_start_mut_index : region_end_mut_index+1])))
        out_str_list.append(",".join(map(str,out_region_dict["cover_count_list"][region_start_mut_index : region_end_mut_index+1])))
        out_str_list.append(",".join(map(str,out_region_dict["mut_ratio_list"][region_start_mut_index : region_end_mut_index+1])))
        out_str_list.append(",".join(map(str,out_region_dict["SNP_ann_list"][region_start_mut_index : region_end_mut_index+1]))) 

        # make output str
        out_str_list.append(",".join(map(str,out_region_dict["tandem_state"][region_start_mut_index : region_end_mut_index+1])))  
        out_str = "\t".join(map(str,out_str_list))

        return (out_str)


###############################################################################
# read parameters
###############################################################################
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="merge pmat file")

    parser.add_argument("-i", "--Input",
        help="Input bmat file",required=True)

    parser.add_argument("-o", "--Output",
        help="Output BED format file",default="Stdout")

    parser.add_argument("-f", "--FromBase",
        help="Ref base, accept A,G,C,T default=C",default="C")

    parser.add_argument("-t", "--ToBase",
        help="Mut base, accept A,G,C,T default=T",default="T")

    parser.add_argument("-r","--reference",
        help="Reference genome fasta file",required=True)

    parser.add_argument("-d", "--MaxSiteDistance",
        help="Max distance between two sites in one region, default=50",default="50")

    parser.add_argument("-D", "--MaxRegionDistance",
        help="Max length of a mutation region, default=100",default="100")

    parser.add_argument("--NoMutNumCutoff",
        help="The number of site without mutation --ToBase signal in a mutation region, default=2",default="2")

    parser.add_argument("--OmitTandemNumCutoff",
        help="The omit tande site cutoff, default=2",default="2")

    parser.add_argument("--SNP",
        help="SNP file with vcf or bed format, if use multiple file, use ',' to separate, default=None",default="None")

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

    if_output_header = eval(ARGS.OutHeader)

    from_base = ARGS.FromBase
    to_base = ARGS.ToBase

    bisite_distance_cutoff = int(ARGS.MaxSiteDistance)
    region_distance_cutoff = int(ARGS.MaxRegionDistance)
    omit_tandem_num_cutoff = int(ARGS.OmitTandemNumCutoff)
    region_no_mut_num_cutoff = int(ARGS.NoMutNumCutoff)

    genome_fa_file_path = ARGS.reference
    SNP_vcf_file_path = ARGS.SNP

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # open input file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ".gz" == input_file_path[-3:]:
        in_pmat_file = gzip.open(input_file_path,"r")
    else:
        in_pmat_file = open(input_file_path,"r")

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
    # load SNP VCF file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # log
    SNP_vcf_dict = {}

    if SNP_vcf_file_path != "None":
        # log
        sys.stderr.write("Loading SNP file...\t%s\n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

        SNP_filename_list = SNP_vcf_file_path.split(",")

        for SNP_filename in SNP_filename_list:
            # log
            sys.stderr.write("%s...\t%s\n" % (SNP_filename,time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

            if ".gz" == SNP_filename[-3:]:
                SNP_vcf_file_obj = gzip.open(SNP_filename,"r")
                SNP_filename_base = os.path.splitext(SNP_filename)[0]
            else:
                SNP_vcf_file_obj = open(SNP_filename,"r")
                SNP_filename_base = SNP_filename

            if SNP_filename_base[-3:].upper() == "VCF":
                SNP_vcf_dict = load_VCF_as_dict(SNP_vcf_file_obj, SNP_vcf_dict)

            elif SNP_filename_base[-3:].upper() == "BED":
                SNP_vcf_dict = load_BED_as_dict(SNP_vcf_file_obj, SNP_vcf_dict)

        sys.stderr.write("Done!\t%s\n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    else:
        sys.stderr.write("No SNP VCF file provided.\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load genome 
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    genome_dict = {}
    genome_fa =  SeqIO.parse(handle=genome_fa_file_path,format="fasta")

    sys.stderr.write("Loading genome... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    for ref in genome_fa:
        sys.stderr.write("Loading...\t" + ref.id + "\n")
        genome_dict[ref.id] = ref.seq.upper()

    sys.stderr.write("Loading genome... Done! \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

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
            "tandem_info"
        ]
        output_file.write("\t".join(header) + "\n")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # header output
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ARGS.InHeader == "True":
        header = in_pmat_file.readline()

    # set initial var    
    region_dict = {}
    vcf_temp_line = None
    region_detect_state = False
    tandem_site_index_list = []
    tandem_order_index = 0
    pmat_line = in_pmat_file.readline()

    # log
    sys.stderr.write("Start to merge pmat file... \t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))

    while pmat_line:
        if not region_detect_state:
            # initial the line info
            region_query_res = get_region_start_line(in_pmat_file, from_base, to_base, pmat_line)
            pmat_line = region_query_res[1]

            if region_query_res[0]:
                pmat_line_list = pmat_line.strip().split("\t")

                region_dict = {
                    "chr_name" : pmat_line_list[0],
                    "chr_start" : pmat_line_list[1],
                    "chr_end" : None,
                    "site_idx_list" : [ pmat_line_list[3] ],
                    "mut_count_list" : [ pmat_line_list[12] ],
                    "cover_count_list" : [ pmat_line_list[13] ],
                    "mut_ratio_list" : [ pmat_line_list[14] ],
                    "SNP_ann_list" : [ ],
                    "miss_tandem_num" : 0,
                    "no_mut_num": 0,
                    "mut_signal_site_num" : 1,
                    "tandem_state": [ 0 ]
                }

                # get snp query res
                if SNP_vcf_dict:
                    query_snp_res = SNP_vcf_dict.get(pmat_line_list[3])
                    if query_snp_res == None:
                        region_dict["SNP_ann_list"].append(False)
                    else:
                        region_dict["SNP_ann_list"].append(True)
                else:
                    region_dict["SNP_ann_list"].append(False)

                # get tandem mut site 
                tandem_site_index_list = back_trandem_site_index(
                    chr_name = region_dict["chr_name"],
                    chr_start = region_dict["chr_start"],
                    base_type = from_base,
                    ref_dict = genome_dict,
                    distance_size_cutoff = region_distance_cutoff,
                    back_site_num_cutoff = 20
                )

                # 如果返回的tandem site 数目很少，则直接输出最终的结果
                if tandem_site_index_list != None:
                    if len(tandem_site_index_list) == 1:
                        out_str = report_region(region_dict, to_base)
                        if out_str:
                            output_file.write(out_str + "\n")

                        # set initial var    
                        region_dict = {}
                        region_detect_state = False
                        tandem_site_index_list = []
                        tandem_order_index = 0 
                        pmat_line = in_pmat_file.readline()

                    else:
                        region_detect_state = True
                        tandem_order_index = 1
                        pmat_line = in_pmat_file.readline()

        else:
            pmat_line_list = pmat_line.strip().split("\t")

            # 判断是否为相同的染色体
            if pmat_line_list[0] != region_dict["chr_name"]:
                out_str = report_region(region_dict, to_base)
                if out_str:
                    output_file.write(out_str + "\n")

                # set initial var    
                region_dict = {}
                region_detect_state = False
                tandem_site_index_list = []
                tandem_order_index = 0 

            else:
                # 判断是否符合-d 以及 -D的距离cutoff
                if int(pmat_line_list[1]) - int(region_dict["site_idx_list"][-1].split("_")[1]) >= bisite_distance_cutoff:
                    out_str = report_region(region_dict, to_base)
                    if out_str:
                        output_file.write(out_str + "\n")

                    # set initial var    
                    region_dict = {}
                    region_detect_state = False
                    tandem_site_index_list = []
                    tandem_order_index = 0 

                elif (int(region_dict["site_idx_list"][-1].split("_")[1]) - int(region_dict["site_idx_list"][0].split("_")[1])) >= region_distance_cutoff:
                    out_str = report_region(region_dict, to_base)
                    if out_str:
                        output_file.write(out_str + "\n")

                    # set initial var    
                    region_dict = {}
                    region_detect_state = False
                    tandem_site_index_list = []
                    tandem_order_index = 0 

                elif tandem_order_index >= len(tandem_site_index_list):
                    out_str = report_region(region_dict, to_base)
                    if out_str:
                        output_file.write(out_str + "\n")

                    # set initial var    
                    region_dict = {}
                    region_detect_state = False
                    tandem_site_index_list = []
                    tandem_order_index = 0 

                else:
                    if pmat_line_list[9] == from_base:
                        cur_tandem_site_index_split = tandem_site_index_list[tandem_order_index].split("_")
                        
                        # 判断是否为tandem info
                        if pmat_line_list[1] == cur_tandem_site_index_split[1]:

                            #判断是否带有mut 信息
                            if int(pmat_line_list[12]) == 0 :
                                region_dict["no_mut_num"] += 1
                            elif pmat_line_list[10] != to_base:
                                region_dict["no_mut_num"] += 1

                            if region_dict["no_mut_num"] >= region_no_mut_num_cutoff:
                                out_str = report_region(region_dict, to_base)
                                if out_str:
                                    output_file.write(out_str + "\n")
                                
                                # set initial var    
                                region_dict = {}
                                region_detect_state = False
                                tandem_site_index_list = []
                                tandem_order_index = 0
                                pmat_line = in_pmat_file.readline()

                            else:
                                if pmat_line_list[10] == to_base:
                                    # add info into region_dict
                                    #### get snp query res
                                    if SNP_vcf_dict:
                                        query_snp_res = SNP_vcf_dict.get(pmat_line_list[3])
                                        if query_snp_res == None:
                                            region_dict["SNP_ann_list"].append(False)
                                        else:
                                            region_dict["SNP_ann_list"].append(True)
                                    else:
                                        region_dict["SNP_ann_list"].append(False)

                                    #### add info into region_dict
                                    region_dict["site_idx_list"].append(pmat_line_list[3])
                                    region_dict["mut_count_list"].append(pmat_line_list[12])
                                    region_dict["cover_count_list"].append(pmat_line_list[13])
                                    region_dict["mut_ratio_list"].append(pmat_line_list[14])
                                    region_dict["tandem_state"].append(0)

                                    if int(pmat_line_list[12]) > 0:
                                        region_dict["mut_signal_site_num"] += 1

                                    pmat_line = in_pmat_file.readline()
                                    tandem_order_index += 1
                                    
                                else:
                                    # 具有非C-T mutation 的信息
                                    if SNP_vcf_dict:
                                        query_snp_res = SNP_vcf_dict.get(pmat_line_list[3])
                                        if query_snp_res == None:
                                            region_dict["SNP_ann_list"].append(False)
                                        else:
                                            region_dict["SNP_ann_list"].append(True)
                                    else:
                                        region_dict["SNP_ann_list"].append(False)
                                        
                                    #### add info into region_dict
                                    region_dict["site_idx_list"].append(pmat_line_list[3])
                                    region_dict["mut_count_list"].append(0)
                                    region_dict["cover_count_list"].append(pmat_line_list[13])
                                    region_dict["mut_ratio_list"].append(0)
                                    region_dict["tandem_state"].append(0)

                                    pmat_line = in_pmat_file.readline()
                                    tandem_order_index += 1
                                    
                        elif int(pmat_line_list[1]) > int(cur_tandem_site_index_split[1]):
                            region_dict["miss_tandem_num"] += 1

                            if region_dict["miss_tandem_num"] > omit_tandem_num_cutoff:
                                out_str = report_region(region_dict, to_base)
                                if out_str:
                                    output_file.write(out_str + "\n")
                                # set initial var    
                                region_dict = {}
                                region_detect_state = False
                                tandem_site_index_list = []
                                tandem_order_index = 0
                            else:
                                # 没有测到的tandem info
                                region_dict["site_idx_list"].append(tandem_site_index_list[tandem_order_index] + ".")
                                region_dict["mut_count_list"].append(0)
                                region_dict["cover_count_list"].append(0)
                                region_dict["mut_ratio_list"].append(0)
                                region_dict["no_mut_num"] += 1
                                region_dict["tandem_state"].append(1)

                                # get snp query res
                                region_dict["SNP_ann_list"].append(False)

                                tandem_order_index += 1

                    else:
                        pmat_line = in_pmat_file.readline()

    sys.stderr.write("Merge done!\t %s \n" % (time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())))    
    

in_pmat_file.close()
output_file.close()








# 2019-09-12