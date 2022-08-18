#! /Volumes/data_HD/miniconda3/bin/python
# _*_ coding: UTF-8 _*_

import argparse
import pysam
import pysamstats
import pandas as pd
import random
import string
import os
import multiprocessing
import logging
import sys

# Version information START --------------------------------------------------
VERSION_INFO = """
Author: MENG Howard

Version-01:
    2022-07-24
        From bam to count info like .pmat / .bmat file
        
E-Mail: menghaowei@pku.edu.cn
"""
# Version information END ----------------------------------------------------


# Learning Part START --------------------------------------------------------
LEARNING_PART = """
output pmat format like

<chr_name> <chr_index> <site_index> <A> <G> <T> <C> <type> <ref_base> <mut_base> <ref_num> <mut_num> <cover_num> <mut_ratio>
chr1    54795   54795   chr1_54795_TC   1   0   2   6   TC  T   C   6   2   8   0.25

TAB separate
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
    pass

    return full_cmd_str


def get_bmat_info(row_info, bed_like_format=True):
    """
    INPUT:
        <row_info>
            pd.Series

        <bed_list_format>
            bool, Ture means the first 3 columns of output is chr_name, index, index; while False means chr_name, index

    RETURN:
        <bmat_format_list>
            list,
                "chr_name",  "chr_index", "ref_base", "A", "G", "C", "T", "del_count", "insert_count",
                "ambiguous_count", "deletion",  "insertion", "ambiguous",  "mut_num"

                The pysamstats lib doesn't support load INDEL detials, so we set INDEL info as '.'
    """

    # return part
    if bed_like_format:
        return_list = (
            row_info["chrom"],
            row_info["pos"],
            row_info["pos"],
            row_info["ref"],
            row_info["A_pp"],
            row_info["G_pp"],
            row_info["T_pp"],
            row_info["C_pp"],
            row_info["deletions_pp"],
            row_info["insertions_pp"],
            0,
            ".",
            ".",
            ".",
            row_info["mismatches_pp"]
        )
    else:
        return_list = (
            row_info["chrom"],
            row_info["pos"],
            row_info["ref"],
            row_info["A_pp"],
            row_info["G_pp"],
            row_info["T_pp"],
            row_info["C_pp"],
            row_info["deletions_pp"],
            row_info["insertions_pp"],
            0,
            ".",
            ".",
            ".",
            row_info["mismatches_pp"]
        )

    return return_list


def get_pmat_info(row_info,
                  base_list=("A_pp", "T_pp", "C_pp", "G_pp", "N_pp"),
                  bed_like_format=True):
    """
    INPUT:
        <row_info>
            pd.Series

        <base_list>
            list, selected row name,

        <bed_list_format>
            bool, Ture means the first 3 columns of output is chr_name, index, index; while False means chr_name, index

    RETURN:
        <pmat_format_list>
            list,
                 <chr_name> <chr_index>  <chr_index>  <site_index>   <A> <G> <T> <C> <type> <ref_base> <mut_base>
                 <ref_num> <mut_num> <cover_num> <mut_ratio>
                   chr1       54795          54795   chr1_54795_TC   1   0   2   6   TC  T   C   6   2   8   0.25

            ! Note, mut_num only count the largest mutation type,
              this number could not equals to the mismatch number in bmat file !
    """

    # count mut
    ref_num = mut_num = 0
    ref_base = mut_base = row_info["ref"]

    for base_key in base_list:
        if row_info.get(base_key) is not None:
            if base_key[0] == ref_base:
                ref_num = row_info[base_key]

            else:
                if row_info[base_key] > mut_num:
                    mut_base = base_key[0]
                    mut_num = row_info[base_key]

    # count total
    cover_num = mut_num + ref_num

    if cover_num > 0:
        mut_ratio = round(mut_num / 1.0 / cover_num, 5)
    else:
        mut_ratio = 0

    # mut type
    if mut_base == ref_base:
        mut_type = ref_base + "."
    else:
        mut_type = ref_base + mut_base

    # get chr info
    chr_name = row_info["chrom"]
    chr_index = row_info["pos"]
    site_index = "%s_%s_%s" % (chr_name, chr_index, mut_type)

    # return part
    if bed_like_format:
        return_list = (
            chr_name,
            chr_index,
            chr_index,
            site_index,
            row_info["A_pp"],
            row_info["G_pp"],
            row_info["T_pp"],
            row_info["C_pp"],
            mut_type,
            ref_base,
            mut_base,
            ref_num,
            mut_num,
            cover_num,
            mut_ratio
        )
    else:
        return_list = (
            chr_name,
            chr_index,
            site_index,
            row_info["A_pp"],
            row_info["G_pp"],
            row_info["T_pp"],
            row_info["C_pp"],
            mut_type,
            ref_base,
            mut_base,
            ref_num,
            mut_num,
            cover_num,
            mut_ratio
        )

    return return_list


def get_BAM_ref_length(bam_filename):
    """
    INPUT:
        <bam_filename>
            BAM file path with .bai index file in a same dir

    RETURN
        <ref_name_list>
            return ref name list, which contain more than ref_count_cutoff mapped reads

    """
    bam_file = pysam.AlignmentFile(bam_filename, "r")
    ref_dict = {
        "key_order": []
    }

    chr_name_list = bam_file.references
    chr_len_list = bam_file.lengths

    for index, chr_name in enumerate(chr_name_list):
        ref_dict["key_order"].append(chr_name)
        ref_dict[chr_name] = chr_len_list[index]

    bam_file.close()

    return ref_dict


def split_genome_into_bin_list(ref_genome_len_dict, exclude_chr_list=(), binsize=20000000):
    """
    INPUT:
        <ref_genome_len_dict>
            dict, generate from FUN <get_BAM_ref_length>

        <exclude_chr_list>
            list or tuple, like ["chr1", "chr2"], item in this list will be excluded in output

        <binsize>
            int, set genome bin size

    RETURN:
        <ref_genome_bin_list>
            list, each row contain chr_name, start_index, end_index, region_index
    """

    ref_genome_bin_list = []

    for chr_name in ref_genome_len_dict["key_order"]:
        chr_len = ref_genome_len_dict[chr_name]

        if chr_name in exclude_chr_list:
            continue

        for start_index in range(0, chr_len, binsize):

            if (start_index + binsize) > (chr_len + 1):
                end_index = chr_len + 1
            else:
                end_index = start_index + binsize

            if start_index == 0:
                start_index = 1

            region_index = "%s_%s_%s" % (chr_name, start_index, end_index)

            ref_genome_bin_list.append([
                chr_name,
                start_index,
                end_index,
                region_index
            ])

    return ref_genome_bin_list


def generate_temp_filename(genome_bin_list, temp_dir="."):
    """
    INPUT:
        <bin_list>
            list, generated from FUN <split_genome_into_bin_list>


    RETURN:
        <temp_file_dict>
            dict, key is region_index, value is a filename

    """

    temp_dir_abs = os.path.abspath(temp_dir)

    file_dict = {}

    for run_chr_name, run_region_start, run_region_end, region_id in genome_bin_list:
        temp_file_basename = "temp__" + region_id + "__" + "".join(
            random.sample(string.ascii_letters + string.digits, 16)
        )

        temp_file_name = os.path.join(temp_dir_abs, temp_file_basename)

        file_dict[region_id] = temp_file_name

    return file_dict


def convert_BAM_to_count_info(
        bam_filename,
        ref_fa_filename,
        region_chr_name,
        query_region_start,
        query_region_end,
        out_filename,
        out_format="pmat",
        bed_like_format=True,
        mapq_cutoff=20,
        base_cutoff=20,
        max_base_depth=8000,
        every_pos_base=False,
        cover_num_cutoff=0,
        mut_num_cutoff=0,
        mut_ratio_cutoff=0,
        mut_type="ALL"
):
    """
    INPUT:
        <bam_filename>
            str, filename

        <ref_fa_filename>
            str, ref FASTA filename

        <region_chr_name>  <region_start>  <region_end>
            str,int,int like chr1 100000000 120000000

        <out_filename>
            str, output filename

        <out_format>
            str, bmat / pmat

        <bed_like_format>
            bool, True / False

        <mapq_cutoff>
            int, 20

        <base_cutoff>
            int, 20

        <every_pos_base>
            bool, True means out all genome bases, even if there is no alignment.

    RETURN
        <out_state>
            int, 0, everything is okay

    """

    # open bam file and ref fasta
    load_bam_obj = pysam.AlignmentFile(bam_filename, "r")
    load_ref_fa = pysam.Fastafile(ref_fa_filename)

    temp_count_res = pysamstats.load_variation(alignmentfile=load_bam_obj,
                                               fafile=load_ref_fa,
                                               chrom=region_chr_name,
                                               start=query_region_start,
                                               end=query_region_end,
                                               one_based=True,
                                               truncate=True,
                                               max_depth=max_base_depth,
                                               min_mapq=mapq_cutoff,
                                               min_baseq=base_cutoff,
                                               pad=every_pos_base)

    temp_count_df = pd.DataFrame(temp_count_res)

    # decoding with UTF-8
    temp_count_df['chrom'] = temp_count_df['chrom'].apply(lambda s: s.decode('utf-8'))
    temp_count_df['ref'] = temp_count_df['ref'].apply(lambda s: s.decode('utf-8'))

    if out_format == "bmat":
        temp_out_df = temp_count_df.apply(
            lambda x: get_bmat_info(x, bed_like_format=bed_like_format),
            axis=1).apply(pd.Series)

    elif out_format == "pmat":
        header_not_bed = [
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

        header_like_bed = [
            "chr_name",
            "chr_index",
            "chr_index2",
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

        temp_out_df = temp_count_df.apply(lambda x: get_pmat_info(x, bed_like_format=bed_like_format), axis=1)

        if bed_like_format:
            temp_out_df = temp_out_df.apply(pd.Series, index=header_like_bed)
        else:
            temp_out_df = temp_out_df.apply(pd.Series, index=header_not_bed)

        # select mut_type
        if mut_type != "ALL":
            temp_out_df = temp_out_df.loc[temp_out_df["mut_type"] == mut_type]

        # filter cutoff
        if (mut_num_cutoff > 0) or (mut_ratio_cutoff > 0) or (cover_num_cutoff > 0):
            temp_out_df = temp_out_df.loc[
                (temp_out_df["cover_num"] >= cover_num_cutoff) &
                (temp_out_df["mut_num"] >= mut_num_cutoff) &
                (temp_out_df["mut_ratio"] >= mut_ratio_cutoff)
            ]

    else:
        raise IOError("out_format is wrong! FUN: convert_BAM_to_count_info")

    # output
    temp_out_df.to_csv(path_or_buf=out_filename, sep="\t", header=False, index=False, encoding="utf-8")

    # delete
    del temp_out_df

    # close file
    load_bam_obj.close()
    load_ref_fa.close()

    return 0


def merge_split_files(input_bin_list,
                      temp_file_dict,
                      out_filename="stdout",
                      header_list=None,
                      in_sep="\t",
                      out_sep="\t",
                      log_verbose=3):
    """
    INPUT:
        <input_bin_list>
            list, typical format like
                 [
                     ['chr1', 10000000, 10010000, 'chr1_10000000_10010000'],
                     ['chr1', 10010000, 10020000, 'chr1_10010000_10020000'],
                     ['chr1', 10020000, 10030000, 'chr1_10020000_10030000'],
                     ['chr1', 10030000, 10040000, 'chr1_10030000_10040000'],
                     ['chr1', 10040000, 10050000, 'chr1_10040000_10050000'],
                     ['chr1', 10050000, 10060000, 'chr1_10050000_10060000'],
                     ['chr1', 10060000, 10070000, 'chr1_10060000_10070000'],
                     ['chr1', 10070000, 10080000, 'chr1_10070000_10080000'],
                     ['chr1', 10080000, 10090000, 'chr1_10080000_10090000'],
                     ['chr1', 10090000, 10100000, 'chr1_10090000_10100000']
                 ]

        <temp_file_dict>
            dict, typical format like
            {
                'chr1_10000000_10010000': '/Users/meng/data_HD/01.temp/temp__chr1_10000000_10010000__91jBQNPoKeq48rlE',
                'chr1_10010000_10020000': '/Users/meng/data_HD/01.temp/temp__chr1_10010000_10020000__j5Dvtzn90EF76wHR',
                'chr1_10020000_10030000': '/Users/meng/data_HD/01.temp/temp__chr1_10020000_10030000__OzqXUTkbxjWDeKBr',
                'chr1_10030000_10040000': '/Users/meng/data_HD/01.temp/temp__chr1_10030000_10040000__xLcQqUHh5j2Tds8M',
                'chr1_10040000_10050000': '/Users/meng/data_HD/01.temp/temp__chr1_10040000_10050000__FXa9uKQ5bYcsnV63',
                'chr1_10050000_10060000': '/Users/meng/data_HD/01.temp/temp__chr1_10050000_10060000__sch51e9E7mnINgOJ',
                'chr1_10060000_10070000': '/Users/meng/data_HD/01.temp/temp__chr1_10060000_10070000__xqfC61UPlAy93bJs',
                'chr1_10070000_10080000': '/Users/meng/data_HD/01.temp/temp__chr1_10070000_10080000__a4xdNpiJn7LyI2Kf',
                'chr1_10080000_10090000': '/Users/meng/data_HD/01.temp/temp__chr1_10080000_10090000__zoWim4DVGuyTCbh6',
                'chr1_10090000_10100000': '/Users/meng/data_HD/01.temp/temp__chr1_10090000_10100000__Y17mqoMAxR9sDHEV'
            }


        <out_filename>
            str

        <header_list>
            list, length have to match with split files

        <sep>
            str, default '\t', can set by custom

    RETURN:
        0, everything is okay.

    INFO:
        Final date 2022-07-26
    """

    # log
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # open output file
    if out_filename == "stdout":
        out_file = sys.stdout
    else:
        out_file = open(out_filename, "w")

    # header output
    if header_list is not None:
        out_file.write(out_sep.join(map(str, header_list)) + "\n")

    # merge files
    logging.info("Starting to merge temp files...")

    for run_info in input_bin_list:
        run_region_index = run_info[3]
        run_temp_filename = temp_file_dict.get(run_region_index)

        if run_temp_filename is None:
            logging.error("Can not find temp file for %s" % run_region_index)
            raise IOError("merge temp file error!")

        # merge step
        logging.debug("Merging files, processing on \n\t%s" % run_temp_filename)

        with open(run_temp_filename, "r") as run_file:
            for line in run_file:
                if in_sep != out_sep:
                    line_list = line.strip().split(in_sep)
                    out_file.write(out_sep.join(line_list) + "\n")

                else:
                    out_file.write(line)

    # log
    logging.info("Merging files, done!")

    # close file
    if out_filename != "stdout":
        out_file.close()

    return 0


def clear_temp_files_by_dict(temp_file_dict, log_verbose=3):
    """
    INPUT:
        <temp_file_dict>
            dict, format like:
            {
                'chr1_10000000_10010000': '/Users/meng/data_HD/01.temp/temp__chr1_10000000_10010000__91jBQNPoKeq48rlE',
                'chr1_10010000_10020000': '/Users/meng/data_HD/01.temp/temp__chr1_10010000_10020000__j5Dvtzn90EF76wHR',
                'chr1_10020000_10030000': '/Users/meng/data_HD/01.temp/temp__chr1_10020000_10030000__OzqXUTkbxjWDeKBr',
                'chr1_10030000_10040000': '/Users/meng/data_HD/01.temp/temp__chr1_10030000_10040000__xLcQqUHh5j2Tds8M',
                'chr1_10040000_10050000': '/Users/meng/data_HD/01.temp/temp__chr1_10040000_10050000__FXa9uKQ5bYcsnV63',
                'chr1_10050000_10060000': '/Users/meng/data_HD/01.temp/temp__chr1_10050000_10060000__sch51e9E7mnINgOJ',
                'chr1_10060000_10070000': '/Users/meng/data_HD/01.temp/temp__chr1_10060000_10070000__xqfC61UPlAy93bJs',
                'chr1_10070000_10080000': '/Users/meng/data_HD/01.temp/temp__chr1_10070000_10080000__a4xdNpiJn7LyI2Kf',
                'chr1_10080000_10090000': '/Users/meng/data_HD/01.temp/temp__chr1_10080000_10090000__zoWim4DVGuyTCbh6',
                'chr1_10090000_10100000': '/Users/meng/data_HD/01.temp/temp__chr1_10090000_10100000__Y17mqoMAxR9sDHEV'
            }

    RETURN:
        0, everything is okay.

    INFO:
        Final date 2022-07-26
    """

    # log
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    run_state = 0

    for run_key in temp_file_dict.keys():
        if type(temp_file_dict[run_key]) is str:
            if os.path.exists(temp_file_dict[run_key]):
                if os.path.isfile(temp_file_dict[run_key]):
                    try:
                        os.remove(temp_file_dict[run_key])
                    except:
                        logging.error("Removing error on \n\t%s" % temp_file_dict[run_key])
                        run_state = 1

    return run_state


# define header info
header_pmat_list = [
    "chr_name", "chr_index", "site_index", "A", "G", "C", "T",
    "mut_type", "ref_base", "mut_base", "ref_num", "mut_num", "cover_num", "mut_ratio"
]

header_pmat_bed_list = [
    "chr_name", "chr_index", "chr_index2", "site_index", "A", "G", "C", "T",
    "mut_type", "ref_base", "mut_base", "ref_num", "mut_num", "cover_num", "mut_ratio"
]

header_bmat_list = [
    "chr_name", "chr_index", "ref_base", "A", "G", "C", "T",
    "del_count", "insert_count", "ambiguous_count", "deletion",  "insertion", "ambiguous",  "mut_num"
]

header_bmat_bed_list = [
    "chr_name", "chr_index", "ref_base", "A", "G", "C", "T",
    "del_count", "insert_count", "ambiguous_count", "deletion",  "insertion", "ambiguous",  "mut_num"
]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="convert bam file to bmat / pmat file")

    parser.add_argument("-i", "--input_bam",
                        help="Input bam file, sorted and indexed", required=True)

    parser.add_argument("-o", "--output",
                        help="Output filename, default=out.pmat", default="out.pmat")

    parser.add_argument("-r", "--reference",
                        help="Genome FASTA file", required=True)

    parser.add_argument("-p", "--threads",
                        help="multiple threads number, default=1", default=1, type=int)

    parser.add_argument("--out_format",
                        help="bmat or pmat, default=pmat", default="pmat")

    parser.add_argument("--bed_like_format",
                        help="The first 3 cols is chr_name, chr_index, chr_index, default=True", default="True")

    parser.add_argument("--block_size",
                        help="Genome block size, larger block size means more memory requirement. "
                             "default=100000 this value needs ~700MB memory for each thread.",
                        default=100000, type=int)

    # ==========================================================================================>>>>>
    # count part
    # ==========================================================================================>>>>>
    parser.add_argument("--mapq_cutoff",
                        help="MAPQ >= this cutoff will be kept. Default=20", default=20, type=int)

    parser.add_argument("--base_cutoff",
                        help="Sequencing base quality cutoff. Default=20", default=20, type=int)

    parser.add_argument("--max_depth",
                        help="Max depth for each position. Default=8000", default=8000, type=int)

    parser.add_argument("--every_position",
                        help="Set True means output every position in whole genome. Default=False",
                        default="False")

    # ==========================================================================================>>>>>
    # out filter part
    # ==========================================================================================>>>>>
    parser.add_argument("--cover_num_cutoff",
                        help="Only for pmat, site coverage number cutoff default=0", default=0, type=int)

    parser.add_argument("--mut_num_cutoff",
                        help="Only for pmat, site mutation number cutoff default=0", default=0, type=int)

    parser.add_argument("--mut_ratio_cutoff",
                        help="Only for pmat, site mutation ratio cutoff default=0", default=0.0, type=float)

    parser.add_argument("--mut_type",
                        help="Only for pmat, select mutation type, ALL means no selection, "
                             "can set like CT, default output all info",
                        default="ALL")

    parser.add_argument("--out_header",
                        help="If contain header line in output file, default=True", default="True")

    parser.add_argument("--keep_temp_file",
                        help="If keep temp files, default=False",
                        default="False")

    parser.add_argument("--verbose",
                        help="Larger number means out more log info, can be 0,1,2,3 default=3",
                        default=3, type=int)

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # load args
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    if ARGS.every_position == "True":
        args_every_position = True
    else:
        args_every_position = False

    if ARGS.out_header == "True":
        args_out_header = True
    else:
        args_out_header = False

    if ARGS.bed_like_format == "True":
        args_bed_like_format = True
    else:
        args_bed_like_format = False

    if ARGS.keep_temp_file == "True":
        args_keep_temp_file = True
    else:
        args_keep_temp_file = False

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # split genome
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    ref_genome_length_dict = get_BAM_ref_length(bam_filename=ARGS.input_bam)

    genome_bin_list = split_genome_into_bin_list(ref_genome_length_dict, binsize=ARGS.block_size)

    # genome_bin_list = [
    #      ['chr1', 10000000, 10010000, 'chr1_10000000_10010000'],
    #      ['chr1', 10010000, 10020000, 'chr1_10010000_10020000'],
    #      ['chr1', 10020000, 10030000, 'chr1_10020000_10030000'],
    #      ['chr1', 10030000, 10040000, 'chr1_10030000_10040000'],
    #      ['chr1', 10040000, 10050000, 'chr1_10040000_10050000'],
    #      ['chr1', 10050000, 10060000, 'chr1_10050000_10060000'],
    #      ['chr1', 10060000, 10070000, 'chr1_10060000_10070000'],
    #      ['chr1', 10070000, 10080000, 'chr1_10070000_10080000'],
    #      ['chr1', 10080000, 10090000, 'chr1_10080000_10090000'],
    #      ['chr1', 10090000, 10100000, 'chr1_10090000_10100000']
    # ]

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # make temp file
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    if ARGS.output == "out.pmat":
        temp_dir = os.getcwd()
    else:
        temp_dir = os.path.dirname(ARGS.output)

    temp_filename_dict = generate_temp_filename(genome_bin_list, temp_dir=temp_dir)
    # temp_filename_dict = {
    #     'chr1_10000000_10010000': '/Users/meng/data_HD/01.temp/temp__chr1_10000000_10010000__91jBQNPoKeq48rlE',
    #     'chr1_10010000_10020000': '/Users/meng/data_HD/01.temp/temp__chr1_10010000_10020000__j5Dvtzn90EF76wHR',
    #     'chr1_10020000_10030000': '/Users/meng/data_HD/01.temp/temp__chr1_10020000_10030000__OzqXUTkbxjWDeKBr',
    #     'chr1_10030000_10040000': '/Users/meng/data_HD/01.temp/temp__chr1_10030000_10040000__xLcQqUHh5j2Tds8M',
    #     'chr1_10040000_10050000': '/Users/meng/data_HD/01.temp/temp__chr1_10040000_10050000__FXa9uKQ5bYcsnV63',
    #     'chr1_10050000_10060000': '/Users/meng/data_HD/01.temp/temp__chr1_10050000_10060000__sch51e9E7mnINgOJ',
    #     'chr1_10060000_10070000': '/Users/meng/data_HD/01.temp/temp__chr1_10060000_10070000__xqfC61UPlAy93bJs',
    #     'chr1_10070000_10080000': '/Users/meng/data_HD/01.temp/temp__chr1_10070000_10080000__a4xdNpiJn7LyI2Kf',
    #     'chr1_10080000_10090000': '/Users/meng/data_HD/01.temp/temp__chr1_10080000_10090000__zoWim4DVGuyTCbh6',
    #     'chr1_10090000_10100000': '/Users/meng/data_HD/01.temp/temp__chr1_10090000_10100000__Y17mqoMAxR9sDHEV'
    # }

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # multi-threads running part
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # set log info
    logging.basicConfig(level=(4 - ARGS.verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    logging.info("Starting to count BAM info...")

    pool = multiprocessing.Pool(processes=ARGS.threads)

    # record info
    region_count_result = []

    for chr_name, region_start, region_end, region_index in genome_bin_list:

        temp_filename = temp_filename_dict.get(region_index)

        if temp_filename is None:
            raise IOError("Error about <temp_filename_dict>!")

        region_count_result.append(
            pool.apply_async(
                func=convert_BAM_to_count_info,
                args=(
                    ARGS.input_bam,
                    ARGS.reference,
                    chr_name,
                    region_start,
                    region_end,
                    temp_filename,
                    ARGS.out_format,
                    args_bed_like_format,
                    ARGS.mapq_cutoff,
                    ARGS.base_cutoff,
                    ARGS.max_depth,
                    args_every_position,
                    ARGS.cover_num_cutoff,
                    ARGS.mut_num_cutoff,
                    ARGS.mut_ratio_cutoff,
                    ARGS.mut_type,
                )
            )
        )

    pool.close()
    pool.join()

    # out pool result
    region_count_result_state = []
    for index, res in enumerate(region_count_result):
        run_res = res.get()
        region_count_result_state.append(run_res)

    print(region_count_result_state)

    logging.info("temp files preparation step done!")

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # check counting state and merge temp files
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.info("Try to merge temp files...")

    out_header_list = None

    if args_out_header:
        if ARGS.out_format == "pmat":
            if args_bed_like_format:
                out_header_list = header_pmat_bed_list
            else:
                out_header_list = header_pmat_list

        elif ARGS.out_format == "bmat":
            if args_bed_like_format:
                out_header_list = header_bmat_bed_list
            else:
                out_header_list = header_bmat_list

    merge_state = merge_split_files(input_bin_list=genome_bin_list,
                                    temp_file_dict=temp_filename_dict,
                                    out_filename=ARGS.output,
                                    header_list=out_header_list,
                                    in_sep="\t",
                                    out_sep="\t",
                                    log_verbose=3)

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # remove temp file
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    if not args_keep_temp_file:
        rm_temp_state = clear_temp_files_by_dict(temp_filename_dict, log_verbose=ARGS.verbose)

    logging.info("Everything is done!")
    # Final edit date: 2022-07-25


