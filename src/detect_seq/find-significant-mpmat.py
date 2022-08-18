#! /gpfs/user/menghaowei/anaconda2/bin/python
# _*_ coding: UTF-8 _*_

import argparse

from DetectSeqLib.CalculateBackground import *
from DetectSeqLib.CheckAndLoadFiles import *
from DetectSeqLib.OutputAndClearTemp import *
from DetectSeqLib.RegionStatsTest import *

# Version information START --------------------------------------------------
VERSION_INFO = """
Author: MENG Howard

Version-01:
    2019-12-30
        Input ctrl json and treat json 

    2019-12-29
        Fix output part 

    2019-12-28 
        Find significant mpmat signal

Version-02:
    2020-01-27
        Fix some bugs

Version-03:
    2020-08-21
        Fix ABE support

Version-04:
    2020-10-23
        Update code to a new version, this code could analyze CBE, ABE, GBE, PE, DdCBE data...
        
        NEW:
            1. support ctrl hard filter, binomial test to block background sites
            2. support multi-threads calculation
            3. support reads tandem element info
                  
Version-05:
    2020-11-09
        Update code to a new version, this code could analyze CBE, ABE, GBE, PE, DdCBE data...
        
        Fix small bugs and can run from mpmat file with block-info
        
Version-06
    2021-02-13
        Update new method, compare directly without sequencing depth normalization
    
Version-07
    2022-01-06
        Add filter by mpmat SNV/SNP info

Version-08
    2022-07-30
        1. add block function
        2. output with highest signal site
        3. output with site index and block info
        
E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = """
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
    9. added info mpmat file
"""


# Learning Part END-----------------------------------------------------------

############################################################################################
# FUN Part
############################################################################################
def _log_cmd_str(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <full_cmd_str>
            str, record full command str info
    """
    full_cmd_str = """python find-significant-mpmat.py
    --mpmat_table {mpmat_table}
    --output {output}
    --ctrl_BAM {ctrl_BAM}
    --treat_BAM {treat_BAM}
    --reference {reference}
    --thread {thread}
    --query_mutation_type {query_mutation_type}
    --verbose {verbose}
    --keep_temp_file {keep_temp_file}
    --mpmat_filter_info_col_index {mpmat_filter_info_col_index}
    --mpmat_block_info_col_index {mpmat_block_info_col_index}
    --region_block_mut_num_cutoff {region_block_mut_num_cutoff}
    --query_mut_min_cutoff {query_mut_min_cutoff}
    --query_mut_max_cutoff {query_mut_max_cutoff}
    --total_mut_max_cutoff {total_mut_max_cutoff}
    --other_mut_max_cutoff {other_mut_max_cutoff}
    --seq_reads_length {seq_reads_length}
    --scale_reads_count {scale_reads_count}
    --lambda_method {lambda_method}
    --poisson_method {poisson_method}
    """.format(
        mpmat_table=args.mpmat_table,
        output=args.output,
        ctrl_BAM=args.ctrl_BAM,
        treat_BAM=args.treat_BAM,
        reference=args.reference,
        thread=args.thread,
        query_mutation_type=args.query_mutation_type,
        verbose=args.verbose,
        keep_temp_file=args.keep_temp_file,
        mpmat_filter_info_col_index=args.mpmat_filter_info_col_index,
        mpmat_block_info_col_index=args.mpmat_block_info_col_index,
        region_block_mut_num_cutoff=args.region_block_mut_num_cutoff,
        query_mut_min_cutoff=args.query_mut_min_cutoff,
        query_mut_max_cutoff=args.query_mut_max_cutoff,
        total_mut_max_cutoff=args.total_mut_max_cutoff,
        other_mut_max_cutoff=args.other_mut_max_cutoff,
        seq_reads_length=args.seq_reads_length,
        scale_reads_count=args.scale_reads_count,
        lambda_method=args.lambda_method,
        poisson_method=args.poisson_method
    )

    return full_cmd_str


def _check_file_exist(args):
    """
    INPUT:
        <args>
            obj, argparse obj

    RETURN:
        <check_all_state>
            bool, True means all files are exist, False means at least one of files can't pass file check step.

        <check_exist_dict>
            dict, each item contain 3 elements:
                1.input filename
                2.check state
                3. check reason
    """
    # init list
    check_all_state = True
    check_exist_dict = {}

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 1. check ctrl bam file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ctrl_check_res = check_input_bam_file(os.path.abspath(args.ctrl_BAM))

    if ctrl_check_res == 1:
        check_exist_dict["ctrl_BAM"] = [args.ctrl_BAM, False, "BAM not exist!"]
        check_all_state = False

    elif ctrl_check_res == 2:
        check_exist_dict["ctrl_BAM"] = [args.ctrl_BAM, False, "BAM not sorted by coordinate!"]
        check_all_state = False

    elif ctrl_check_res == 3:
        check_exist_dict["ctrl_BAM"] = [args.ctrl_BAM, False, "BAM doesn't contain index file!"]
        check_all_state = False

    else:
        check_exist_dict["ctrl_BAM"] = [args.ctrl_BAM, True, None]

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 2. check PD bam file
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    treat_check_res = check_input_bam_file(os.path.abspath(args.treat_BAM))

    if treat_check_res == 1:
        check_exist_dict["treat_BAM"] = [args.treat_BAM, False, "BAM not exist!"]
        check_all_state = False

    elif treat_check_res == 2:
        check_exist_dict["treat_BAM"] = [args.treat_BAM, False, "BAM not sorted by coordinate!"]
        check_all_state = False

    elif treat_check_res == 3:
        check_exist_dict["treat_BAM"] = [args.treat_BAM, False, "BAM doesn't contain index file!"]
        check_all_state = False

    else:
        check_exist_dict["treat_BAM"] = [args.treat_BAM, True, None]

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # 3. check other file exist
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # mpmat_table
    if not os.path.exists(os.path.abspath(ARGS.mpmat_table)):
        check_exist_dict["mpmat_table"] = [args.mpmat_table, False, "--mpmat_table does not exist!"]
        check_all_state = False
    else:
        check_exist_dict["mpmat_table"] = [args.mpmat_table, True, None]

    # reference
    if not os.path.exists(os.path.abspath(ARGS.reference)):
        check_exist_dict["reference"] = [args.reference, False, "--reference does not exist!"]
        check_all_state = False
    else:
        check_exist_dict["reference"] = [args.reference, True, None]

    # make log
    logging.basicConfig(level=10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    for file_type in ["mpmat_table", "ctrl_BAM", "treat_BAM", "reference"]:
        if check_exist_dict[file_type][1]:
            logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
            logging.info("Yes!")
        else:
            logging.info("Check exist \n \t--%s %s" % (file_type, check_exist_dict[file_type][0]))
            logging.info("No! Reason: %s" % (check_exist_dict[file_type][2]))

        sys.stderr.write("-" * 80 + "\n")

    # return part
    return check_all_state, check_exist_dict


header_list = [
    "chr_name",
    "region_start",
    "region_end",
    "mpmat_index",
    "region_site_num",
    "region_block_site_num",
    "region_mut_site_num",
    "region_site_index",
    "region_block_state",
    "region_highest_site_index",
    "region_highest_site_mut_num",
    "region_highest_site_cover_num",
    "region_highest_site_mut_ratio",
    "ctrl_count",
    "treat_count",
    "ctrl_mut_count",
    "treat_mut_count",
    "ctrl_count.norm",
    "treat_count.norm",
    "ctrl_mut_count.norm",
    "treat_mut_count.norm",
    "count_info",
    "log2_FC",
    "log2_FC_mut",
    "test_state",
    "p_value",
    "FDR"
]

# ---------------------------------------------------------------------------->>>>>>>>>>
#  main part
# ---------------------------------------------------------------------------->>>>>>>>>>
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="A tool to find significant enrichment mpmat region. "
                                                 "You can set a lot of parameters with this program, "
                                                 "usually default parameters can work well.")

    # ==========================================================================================>>>>>
    # Input and output params
    # ==========================================================================================>>>>>
    parser.add_argument("-i", "--mpmat_table",
                        help=".mpmat table, which can be generated from <pmat-merge.py> code",
                        required=True)

    parser.add_argument("-o", "--output",
                        help="Output Poisson test result, default=poisson_output.tsv", default="./poisson_output.tsv",
                        type=str)

    parser.add_argument("-c", "--ctrl_BAM",
                        help="Control BAM file", required=True)

    parser.add_argument("-t", "--treat_BAM",
                        help="Treatment BAM file", required=True)

    parser.add_argument("-r", "--reference",
                        help="Genome FASTA file", required=True)

    parser.add_argument("-p", "--thread",
                        help="Number of threads used to process data", default=1, type=int)

    parser.add_argument("--query_mutation_type",
                        help="Query mutation type, which will be considered as mutation signals, default=CT,GA",
                        default="CT,GA")

    parser.add_argument("--verbose",
                        help="Larger number means out more log info, can be 0,1,2,3 default=3",
                        default=3, type=int)

    parser.add_argument("--keep_temp_file",
                        help="If keep temp files, default=False",
                        default="False")

    # ==========================================================================================>>>>>
    # mpmat file parsing
    # ==========================================================================================>>>>>
    parser.add_argument("--mpmat_filter_info_col_index",
                        help="Column index for filter info, -1 means no such info. Default=-1. Usually can set as 13",
                        default=-1, type=int)

    parser.add_argument("--mpmat_block_info_col_index",
                        help="Column index for region block info, -1 means no such info. Default=-1.",
                        default=-1, type=int)

    # ==========================================================================================>>>>>
    # mpmat region block filter
    # ==========================================================================================>>>>>
    parser.add_argument("--region_block_mut_num_cutoff",
                        help="Site filter cutoff, if a site has a mutation signal >= this cutoff in ctrl sample, "
                             "the site will be blocked in the downstream analysis. Default=2",
                        default=2, type=int)

    # ==========================================================================================>>>>>
    # Poisson test params
    # ==========================================================================================>>>>>
    parser.add_argument("--query_mut_min_cutoff",
                        help="An alignment contains query mutation count lower than this cutoff in mpmat region "
                             "will considered as 'non_mut_align', default=1",
                        default=1, type=int)

    parser.add_argument("--query_mut_max_cutoff",
                        help="An alignment contain mutation number higher than this, considered as "
                             "'query_high_mismatch', which often caused by Bismark mapping error, default=16",
                        default=16, type=int)

    parser.add_argument("--total_mut_max_cutoff",
                        help="An alignment contain mutation number higher than this, considered as "
                             "'total_high_mismatch', default=20",
                        default=20, type=int)

    parser.add_argument("--other_mut_max_cutoff",
                        help="An alignment contain mutation number higher than this, considered as "
                             "'other_high_mismatch', default=12",
                        default=12, type=int)

    parser.add_argument("--seq_reads_length",
                        help="Sequencing reads length, default=150",
                        default=150, type=int)

    parser.add_argument("--scale_reads_count",
                        help="Scaled final output region signal, default=1000000 Usually, default=1e6, means CPM",
                        default=int(1e6), type=int)

    parser.add_argument("--lambda_method",
                        help="Can be set as 'ctrl_max', 'treat_max', 'max', 'raw', default=ctrl_max",
                        default="ctrl_max")

    parser.add_argument("--poisson_method",
                        help="Can be set as 'mutation' OR 'all', default=mutation. "
                             "'mutation' means only use mutation alignments to run Poisson test,"
                             "'all' means use all alignments to run Poisson, which similar to MACS2",
                        default="mutation")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # * load the parameters
    # ---------------------------------------------------------------------------->>>>>>>>>>
    ARGS = parser.parse_args()

    # output full cmd
    sys.stderr.write("\n" + "-" * 80 + "\n")
    sys.stderr.write(_log_cmd_str(args=ARGS))

    # check file exist
    sys.stderr.write("\n" + "-" * 80 + "\n")
    check_res = _check_file_exist(args=ARGS)
    if not check_res[0]:
        raise IOError("Please make sure each file is exist!")

    # check output file dir
    if ARGS.output == "stdout":
        temp_dir = os.getcwd()
    else:
        temp_dir = os.path.dirname(ARGS.output)

    # fix ARGS
    query_mutation_type = ARGS.query_mutation_type.split(",")

    # log format
    logging.basicConfig(level=(4 - ARGS.verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # if ARGS.mpmat_block_info_col_index != 14:
    #     logging.warning("mpmat file block col index not set as default!")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part I split mpmat
    # ---------------------------------------------------------------------------->>>>>>>>>>
    # split
    sys.stderr.write("\n" + "-" * 80 + "\n")
    split_mpmat_filename_dict = split_mpmat_by_chr_name(
        input_mpmat_filename=ARGS.mpmat_table,
        temp_dir=temp_dir,
        force_temp_dir=True)

    # meta dict
    meta_dict = {"filename": {}}
    meta_dict["filename"]["mpmat_region_AddBlockInfo"] = ARGS.mpmat_table
    meta_dict["filename"]["mpmat_region_AddBlockInfo_split"] = split_mpmat_filename_dict.copy()

    # make chr order
    ref_order_dict = {}
    ref_order_list = []

    for index, chr_name in enumerate(meta_dict["filename"]["mpmat_region_AddBlockInfo_split"]["chr_name_order"]):
        ref_order_dict[chr_name] = index
        ref_order_list.append(chr_name)

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part II calculation of normalization scale factor
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    logging.info("Calculating ctrl sample Poisson background...")

    ctrl_poisson_bg_meta_dict = multi_thread_chr_bam_mut_count(
            input_bam_filename=ARGS.ctrl_BAM,
            select_chr_name_list=ref_order_list,
            ref_genome_filename=ARGS.reference,
            thread_num=ARGS.thread,
            query_mut_type_list=query_mutation_type,
            query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
            query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
            total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
            other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
            log_verbose=ARGS.verbose)

    sys.stderr.write("-" * 80 + "\n")
    logging.info("Calculating treat sample Poisson background...")

    treat_poisson_bg_meta_dict = multi_thread_chr_bam_mut_count(
            input_bam_filename=ARGS.treat_BAM,
            select_chr_name_list=ref_order_list,
            ref_genome_filename=ARGS.reference,
            thread_num=ARGS.thread,
            query_mut_type_list=query_mutation_type,
            query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
            query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
            total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
            other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
            log_verbose=ARGS.verbose)

    # genome effective length
    sys.stderr.write("-" * 80 + "\n")
    genome_base_count_dict = calculate_effective_genome(ref_genome_filename=ARGS.reference, log_verbose=ARGS.verbose)

    # make meta data dict
    meta_data_dict = {"ctrl": count_dict_sum_up(ctrl_poisson_bg_meta_dict["count_dict"]),
                      "treat": count_dict_sum_up(treat_poisson_bg_meta_dict["count_dict"])}

    # make scale dict
    scale_factor_dict, normalize_scale_factor_dict, genome_bg_dict = back_all_normalization_scale_dict(
        meta_data_dict=meta_data_dict,
        genome_base_count_dict=genome_base_count_dict,
        seq_reads_length=ARGS.seq_reads_length,
        norm_scale_reads_count=int(ARGS.scale_reads_count))

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part III Poisson test on mpmat region
    # ---------------------------------------------------------------------------->>>>>>>>>>
    poisson_out_split_dict = multi_run_mpmat_poisson_test(
        mpmat_block_split_dict=meta_dict["filename"]["mpmat_region_AddBlockInfo_split"],
        ctrl_bam_filename=ARGS.ctrl_BAM,
        treat_bam_filename=ARGS.treat_BAM,
        ref_genome_fa_filename=ARGS.reference,
        scale_factor_dict=scale_factor_dict,
        normalize_scale_factor_dict=normalize_scale_factor_dict,
        genome_bg_dict=genome_bg_dict,
        lambda_bg_method=ARGS.lambda_method,
        poisson_method=ARGS.poisson_method,
        log_verbose=ARGS.verbose,
        thread=ARGS.thread,
        region_block_mut_num_cutoff=ARGS.region_block_mut_num_cutoff,
        reads_query_mut_min_cutoff=ARGS.query_mut_min_cutoff,
        reads_query_mut_max_cutoff=ARGS.query_mut_max_cutoff,
        reads_total_mut_max_cutoff=ARGS.total_mut_max_cutoff,
        reads_other_mut_max_cutoff=ARGS.other_mut_max_cutoff,
        mpmat_filter_info_col_index=ARGS.mpmat_filter_info_col_index,
        mpmat_block_info_col_index=ARGS.mpmat_block_info_col_index
    )

    meta_dict["filename"]["Poisson_out_split"] = poisson_out_split_dict.copy()
    meta_dict["filename"]["Poisson_out_merge_NoFDR"] = ARGS.output + ".NoFDR"

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part VI output Poisson test result
    # ---------------------------------------------------------------------------->>>>>>>>>>
    sys.stderr.write("-" * 80 + "\n")
    merge_state, pval_list_str = merge_split_files(split_file_dict=meta_dict["filename"]["Poisson_out_split"],
                                                   key_order_list=ref_order_list,
                                                   out_filename=ARGS.output + ".NoFDR",
                                                   header_list=None,
                                                   in_sep="\t", out_sep="\t",
                                                   log_verbose=ARGS.verbose,
                                                   return_col_index=-1)

    sys.stderr.write("-" * 80 + "\n")
    logging.info("Calculating FDR value and make final output table...")

    # change pval into float
    pval_list = []
    for pval_str in pval_list_str:
        if pval_str != "NA":
            pval_list.append(eval(pval_str))
        else:
            pval_list.append("NA")

    # make qval with BH method
    qval_list = make_qvalue_with_BH_method(pval_list)

    # make final output
    final_out_file = open(ARGS.output, "w")
    final_out_file.write("\t".join(header_list) + "\n")

    with open(meta_dict["filename"]["Poisson_out_merge_NoFDR"], "r") as no_fdr_file:
        for index, line in enumerate(no_fdr_file):
            line_list = line.strip().split("\t")
            line_list.append(str(qval_list[index]))
            final_out_file.write("\t".join(line_list) + "\n")

    final_out_file.close()

    logging.info("Done!")

    # ---------------------------------------------------------------------------->>>>>>>>>>
    # Part VII remove temp files
    # ---------------------------------------------------------------------------->>>>>>>>>>
    if ARGS.keep_temp_file == "False":
        sys.stderr.write("-" * 80 + "\n")
        logging.info("Set <keep_temp_file> as False... Removing temp files...")

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing mpmat split files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["mpmat_region_AddBlockInfo_split"],
                                 log_verbose=ARGS.verbose)

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing Poisson test split files...")
        clear_temp_files_by_dict(temp_file_dict=meta_dict["filename"]["Poisson_out_split"],
                                 log_verbose=ARGS.verbose)

        sys.stderr.write("-" * 80 + "\n")
        logging.info("Removing Poisson test raw file...")
        try:
            os.remove(meta_dict["filename"]["Poisson_out_merge_NoFDR"])
        except Warning as w:
            logging.warning("Removing file error: \n\t%s" % meta_dict["filename"]["Poisson_out_merge_NoFDR"])

    logging.info("Everything done!")
