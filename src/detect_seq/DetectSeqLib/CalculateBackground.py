# _*_ coding: UTF-8 _*_

import logging
import sys
import os
import gzip
import random
import string
import json
import pysam
import time

from Bio import SeqIO
import multiprocessing

from DetectSeqLib.CheckAndLoadFiles import load_reference_fasta_as_dict
from DetectSeqLib.RegionStatsTest import get_align_mismatch_pairs, get_No_MD_align_mismatch_pairs

# Version information START ----------------------------------------------------
VERSION_INFO = \
    """

Author: MENG Howard

Version-01:
  2020-11-04 function part for Detect-seq significant test

E-Mail: meng_howard@126.com

"""
# Version information END ------------------------------------------------------

# Function List START ----------------------------------------------------------
COPY_RIGHT = \
"""
Belong to MENG Haowei, YiLab @ Peking University
"""


# Function List END ------------------------------------------------------------


#################################################################################
# FUN
#################################################################################
def count_effective_genome(genome_fa_filename, genome_json_out=None, log_verbose=3):
    """
    INPUT:
        <genome_fa_filename>
            str, Genome FASTA filename
        
        <genome_json_out>
            str, Output filename with JSON format, default=stdout means print info to screen.
        
        <log_verbose>
            int, log output info level, bigger means more log info
    
    OUTPUT:
        Count each chromosome A,T,C,G,N count, and out info with JSON format.
        
        {
            "chr1": {
                "A": 67070277,
                "C": 48055043,
                "G": 48111528,
                "N": 18475410,
                "T": 67244164
            }
        }
    """

    # --------------------------------------------------->>>>>>
    # log setting 
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    logging.info("-" * 80)
    logging.info("Counting reference genome background info...")
    logging.info("Processing FASTA file...")

    # --------------------------------------------------->>>>>>
    # load genome file
    # --------------------------------------------------->>>>>>
    genome_fa = SeqIO.parse(handle=genome_fa_filename, format="fasta")

    genome_base_count_dict = {
        "total": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}
    }

    for ref in genome_fa:
        logging.debug("Counting with %s" % ref.id)

        chr_seq = ref.seq.upper()

        if genome_base_count_dict.get(ref.id) is None:
            genome_base_count_dict[ref.id] = {}

        # add into dict 
        genome_base_count_dict[ref.id]["A"] = chr_seq.count("A")
        genome_base_count_dict[ref.id]["T"] = chr_seq.count("T")
        genome_base_count_dict[ref.id]["C"] = chr_seq.count("C")
        genome_base_count_dict[ref.id]["G"] = chr_seq.count("G")
        genome_base_count_dict[ref.id]["N"] = chr_seq.count("N")

        # add to total
        genome_base_count_dict["total"]["A"] += genome_base_count_dict[ref.id]["A"]
        genome_base_count_dict["total"]["T"] += genome_base_count_dict[ref.id]["T"]
        genome_base_count_dict["total"]["C"] += genome_base_count_dict[ref.id]["C"]
        genome_base_count_dict["total"]["G"] += genome_base_count_dict[ref.id]["G"]
        genome_base_count_dict["total"]["N"] += genome_base_count_dict[ref.id]["N"]

    # log
    logging.info("Processing FASTA file. Done!")
    logging.info("-" * 80)

    # --------------------------------------------------->>>>>>
    # output json
    # --------------------------------------------------->>>>>>
    genome_base_count_dict["reference"] = os.path.abspath(genome_fa_filename)

    if genome_json_out == "stdout":
        out_file = sys.stdout

    elif genome_json_out is None:
        out_file = None

    else:
        out_file = open(genome_json_out, "w")

    if out_file is not None:
        json_obj = json.dumps(genome_base_count_dict, encoding="utf-8", sort_keys=True, indent=4,
                              separators=(', ', ': '))
        out_file.write(json_obj)

    # close file 
    if (genome_json_out != "stdout") and (out_file is not None):
        out_file.close()

    genome_fa.close()

    # --------------------------------------------------->>>>>>
    # return part
    # --------------------------------------------------->>>>>>
    return genome_base_count_dict


#################################################################################
# FUN
#################################################################################
def count_split_chrom_mut_bg(bmat_filename,
                             ref_name,
                             query_mut_type_list=["CT", "GA"],
                             json_out_filename=None,
                             log_verbose=3):
    """
    INPUT
        <bmat_filename>
            str, .bmat file
        
        <ref_name>
            str, ref_name like chr1
    
    OUTPUT
        <json_out_filename>
            str, default=stdout, output calculation reault in JSON format.
        
        <log_out_filename>
            str, default=stderr, output run log.
    
    RETURN
        dict, count_dict, with the same info to <json_out_filename>
        
        return dict like:

        {
            "A": 523999, 
            "AC": 136, 
            "AG": 553, 
            "AT": 145, 
            "C": 566244, 
            "CA": 593, 
            "CG": 190, 
            "CT": 3845, 
            "G": 540117, 
            "GA": 3405, 
            "GC": 202, 
            "GT": 551, 
            "T": 494954, 
            "TA": 168, 
            "TC": 777, 
            "TG": 154, 
            "query_mut_base_count": 7250, 
            "query_mut_bg_pval": 0.006553014793543879, 
            "query_total_base_count": 1106361, 
            "total_base_count": 2125314, 
            "total_mut_base_count": 10719, 
            "total_mut_bg_pval": 0.005043490044294632
        }

    """
    # --------------------------------------------------->>>>>>
    # log setting 
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ------------------------------------------------------>>>>>>    
    # open file
    # ------------------------------------------------------>>>>>>
    # open json out 
    if json_out_filename == "stdout":
        json_out = sys.stdout

    elif json_out_filename is None:
        json_out = None

    else:
        json_out = open(json_out_filename, "w")

    # open input bmat file
    if (bmat_filename[-3:] == ".gz") or (bmat_filename[-5:] == ".gzip"):
        input_file = gzip.open(bmat_filename, "r")
    else:
        input_file = open(bmat_filename, "r")

    # log
    logging.info("-" * 80)
    logging.info("Start to process bmat file. Set ref_name as: %s" % ref_name)

    # ------------------------------------------------------>>>>>>
    # check header 
    # ------------------------------------------------------>>>>>>    
    header = input_file.readline()

    if not ("chr_name" in header):
        input_file.seek(0)
        logging.warning("Input .bmat file does not contain header line.")

    # ------------------------------------------------------>>>>>>    
    # init count dict
    # ------------------------------------------------------>>>>>>    
    base_count_dict = {
        "total_base_count": 0,
        "total_mut_base_count": 0,
        "total_mut_bg_pval": 0.0001,
        "query_total_base_count": 0,
        "query_mut_base_count": 0,
        "query_mut_bg_pval": 0.0001
    }

    for from_base in "ATCGN":
        for to_base in "ATCGN":
            if from_base == to_base:
                key = from_base
            else:
                key = from_base + to_base

            base_count_dict[key] = 0

    # ------------------------------------------------------>>>>>>    
    # count mutation and cover info
    # ------------------------------------------------------>>>>>>    
    for run_index, line in enumerate(input_file):
        if run_index % 1000000 == 0:
            logging.debug("Processing %s bmat file... %s" % (ref_name, run_index))

        line_list = line.strip().split("\t")

        chr_name = line_list[0]
        ref_base = line_list[2]

        # check chr name
        if chr_name != ref_name:
            logging.error(".bmat file chr_name does not match <ref_name>!")
            raise IOError(".bmat file chr_name does not match <ref_name>!")

        # init record cover num
        cover_base_count = 0

        for index, to_base in enumerate("AGCTN"):
            base_count = int(line_list[3 + index])
            cover_base_count += base_count

            # record mutation info 
            if ref_base != to_base:
                key = ref_base + to_base
                base_count_dict[key] += base_count

        # accumulation total cover
        base_count_dict[ref_base] += cover_base_count

    # ------------------------------------------------------>>>>>>    
    # count query and total, calculate bg p.val
    # ------------------------------------------------------>>>>>>
    for query_mut_type in query_mut_type_list:
        query_ref_base = query_mut_type[0]

        # accumulation query number
        base_count_dict["query_total_base_count"] += base_count_dict[query_ref_base]
        base_count_dict["query_mut_base_count"] += base_count_dict[query_mut_type]

    for ref_base in "ATCGN":
        for to_base in "ATCGN":
            if ref_base == to_base:
                key = ref_base
                base_count_dict["total_base_count"] += base_count_dict[key]

            else:
                key = ref_base + to_base
                base_count_dict["total_mut_base_count"] += base_count_dict[key]

    # calculate background pval
    base_count_dict["total_mut_bg_pval"] = base_count_dict["total_mut_base_count"] / 1.0 / base_count_dict[
        "total_base_count"]
    base_count_dict["query_mut_bg_pval"] = base_count_dict["query_mut_base_count"] / 1.0 / base_count_dict[
        "query_total_base_count"]

    # ------------------------------------------------------>>>>>>    
    # return and output
    # ------------------------------------------------------>>>>>>
    # output JSON result
    if json_out is not None:
        json_obj = json.dumps(base_count_dict, encoding="utf-8", sort_keys=True, indent=4, separators=(', ', ': '))
        json_out.write(json_obj)

        # out log files
    logging.info("Processing bmat file. DONE!")

    # close files
    if (json_out_filename != "stdout") and (json_out is not None):
        json_out.close()

    input_file.close()

    # return 
    return base_count_dict


#################################################################################
# FUN
#################################################################################
# function split files
def split_bmat_by_chr_name(input_bmat_filename, temp_dir=None, force_temp_dir=True, log_verbose=3, out_gzip=False):
    """
    INPUT
        <input_bmat_filename>
            str, .bmat filename, support .gz OR .gzip 
            
        <temp_dir>
            str, a dir to store temp files, None means the same dir with <input_bmat_filename>
        
        <log_verbose>
            int, log output info level, bigger means more log info
                    
    OUTPUT
        Split files by chr_name
    
    RETURN
        A dict, structure like:
        
        dict = {
            "chr1" : "file.chr1.saFSDfjsj91.bmat",
            "chr2" : "file.chr2.saFSDasjfj2.bmat"
            ... ... 
        }
    
    """
    # --------------------------------------------------->>>>>>
    # log setting 
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # set temp dir 
    # --------------------------------------------------->>>>>>
    if temp_dir is None:
        temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
    else:
        temp_dir = os.path.abspath(temp_dir)

    if temp_dir[-1] != "/":
        temp_dir += "/"

    # temp dir check and create
    if not os.path.exists(temp_dir):
        if force_temp_dir:
            logging.warning("<temp_dir> setting is not exist \t %s " % temp_dir)
            logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % temp_dir)

            try:
                os.makedirs(os.path.abspath(temp_dir))
            except:
                logging.warning("Temp dir creating error: \t %s" % temp_dir)
                logging.warning("set <temp_dir> as the same dir with <input_bmat_filename>")
                temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))

        else:
            temp_dir = os.path.abspath(os.path.dirname(input_bmat_filename))
            logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % temp_dir)
    else:
        temp_dir = os.path.abspath(temp_dir)

    # --------------------------------------------------->>>>>>
    # get input basename 
    # --------------------------------------------------->>>>>>
    input_file_basename = os.path.basename(input_bmat_filename)

    # --------------------------------------------------->>>>>>
    # make record dict
    # --------------------------------------------------->>>>>>
    record_dict = {
        "chr_name_order": [],
    }

    # --------------------------------------------------->>>>>>
    # split file
    # --------------------------------------------------->>>>>>
    logging.info("Try to split bmat file...")
    logging.info("Output dir is %s" % temp_dir)

    # open input bmat file
    if (input_bmat_filename[-3:] == ".gz") or (input_bmat_filename[-5:] == ".gzip"):
        input_file = gzip.open(input_bmat_filename, "r")
    else:
        input_file = open(input_bmat_filename, "r")

    # set init 
    cur_chr_name = None
    cur_out_file = None

    for line in input_file:
        line_list = line.strip().split("\t")
        chr_name = line_list[0]

        if chr_name == "chr_name":
            continue

        if cur_chr_name is None:
            # log
            logging.info("Processing %s .bmat file" % chr_name)

            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(temp_dir, temp_file_basename)

            if out_gzip:
                temp_file_name += ".gz"

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            if not out_gzip:
                cur_out_file = open(temp_file_name, "w")
            else:
                cur_out_file = gzip.open(temp_file_name, "w")

            cur_out_file.write(line)

        elif cur_chr_name == chr_name:
            cur_out_file.write(line)

        elif cur_chr_name != chr_name:
            cur_out_file.close()

            # log
            logging.info("Processing %s .bmat file" % chr_name)

            # set next chr_name
            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(temp_dir, temp_file_basename)

            if out_gzip:
                temp_file_name += ".gz"

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            if not out_gzip:
                cur_out_file = open(temp_file_name, "w")
            else:
                cur_out_file = gzip.open(temp_file_name, "w")

            cur_out_file.write(line)

    # close all files
    try:
        cur_out_file.close()
        input_file.close()
    except:
        logging.error("Error occurs at close file step.")
        raise IOError("Error occurs at close file step.")

    # log 
    logging.info("Try to split bmat file. DONE!")

    # --------------------------------------------------->>>>>>
    # return
    # --------------------------------------------------->>>>>>
    return record_dict


#################################################################################
# FUN
#################################################################################
# function
def multi_calculate_SNP_bg_pval(split_bmat_dict, threads=1, query_mut_type_list=["CT", "GA"], log_verbose=3):
    """
    INPUT:
        <split_bmat_dict>
            dict, contain split temp files, format like:
            
            {
                "chr1": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr1.k1LKEjDWSagudoQI", 
                "chr2": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr2.KUrM1AL2dCoHYRlE", 
                "chr3": "./temp_test.merge_chr1_chr2_chr3.1K.bmat.gz.chr3.0TJKF4pHAPErWDMi", 
                "chr_name_order": [
                    "chr1", 
                    "chr2", 
                    "chr3"
                ]
            }
        
        <threads>
            int, set threads number.

        <query_mut_type_list>
            list, inherit from FUN: count_split_chrom_mut_bg
    
    RETURN
        <mutation_bg_dict>
            dict, contain genome mutation background info, key is ref_name, val is mutation info
            
            with format like:
            {
                "chr1":{
                        'A': 523999,
                        'AC': 136,
                        'AG': 553,
                        'AT': 145,
                        'C': 566244,
                        'CA': 593,
                        'CG': 190,
                        'CT': 3845,
                        'G': 540117,
                        'GA': 3405,
                        'GC': 202,
                        'GT': 551,
                        'T': 494954,
                        'TA': 168,
                        'TC': 777,
                        'TG': 154,
                        'query_mut_base_count': 7250,
                        'query_mut_bg_pval': 0.006553014793543879,
                        'query_total_base_count': 1106361,
                        'total_base_count': 2125314,
                        'total_mut_base_count': 10719,
                        'total_mut_bg_pval': 0.005043490044294632          
                }, ...

            }
            
    """

    # --------------------------------------------------->>>>>>
    # log setting 
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # make record dict
    # --------------------------------------------------->>>>>>
    mutation_bg_dict = {
        "chr_name_order": []
    }

    # --------------------------------------------------->>>>>>
    # run part 
    # --------------------------------------------------->>>>>>
    logging.info("-" * 80)
    logging.info("Start to calculate chromosome mutation background...")

    pool = multiprocessing.Pool(processes=threads)

    count_split_result_list = []

    for chr_ref_name in split_bmat_dict["chr_name_order"]:
        split_bmat_filename = split_bmat_dict[chr_ref_name]

        # pool and results
        count_split_result_list.append(
            pool.apply_async(
                func=count_split_chrom_mut_bg,
                args=(
                    split_bmat_filename,
                    chr_ref_name,
                    query_mut_type_list,
                    None,
                    0,
                )
            )
        )

    pool.close()
    pool.join()

    logging.info("Calculating chromosome mutation background done!")
    logging.info("-" * 80)
    # --------------------------------------------------->>>>>>
    # get result and return part 
    # --------------------------------------------------->>>>>>
    for index, res in enumerate(count_split_result_list):
        chr_count_res = res.get().copy()
        chr_name = split_bmat_dict["chr_name_order"][index]

        # add into dict
        mutation_bg_dict[chr_name] = chr_count_res
        mutation_bg_dict["chr_name_order"].append(chr_name)

    return mutation_bg_dict


#################################################################################
# FUN
#################################################################################
# function split files
def split_mpmat_by_chr_name(input_mpmat_filename, temp_dir=None, force_temp_dir=True, log_verbose=3):
    """
    INPUT
        <input_mpmat_filename>
            str, .mpmat filename, support .gz OR .gzip 
            
        <temp_dir>
            str, a dir to store temp files, None means the same dir with <input_mpmat_filename>
        
        <log_verbose>
            int, log output info level, bigger means more log info
                    
    OUTPUT
        Split files by chr_name
    
    RETURN
        A dict, structure like:
        
        dict = {
            "chr1" : "file.chr1.saFSDfjsj91.mpmat",
            "chr2" : "file.chr2.saFSDasjfj2.mpmat"
            ... ... 
        }
    
    """
    # --------------------------------------------------->>>>>>
    # log setting 
    # --------------------------------------------------->>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # --------------------------------------------------->>>>>>
    # set temp dir 
    # --------------------------------------------------->>>>>>
    if temp_dir is None:
        temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))
    else:
        temp_dir = os.path.abspath(temp_dir)

    # add back slash
    if temp_dir[-1] != "/":
        temp_dir += "/"

    # temp dir check and create
    if not os.path.exists(temp_dir):
        if force_temp_dir:
            logging.warning("<temp_dir> setting is not exist \t %s " % temp_dir)
            logging.warning("<force_temp_dir> set as True, try to create temp dir \t %s" % temp_dir)

            try:
                os.makedirs(os.path.abspath(temp_dir))
            except:
                logging.warning("Temp dir creating error: \t %s" % temp_dir)
                logging.warning("set <temp_dir> as the same dir with <input_mpmat_filename>")
                temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))

        else:
            temp_dir = os.path.abspath(os.path.dirname(input_mpmat_filename))
            logging.warning("<temp_dir> setting is not exist, set <temp_dir> as %s" % temp_dir)
    else:
        temp_dir = os.path.abspath(temp_dir)

    # --------------------------------------------------->>>>>>
    # get input basename 
    # --------------------------------------------------->>>>>>
    input_file_basename = os.path.basename(input_mpmat_filename)

    # --------------------------------------------------->>>>>>
    # make record dict
    # --------------------------------------------------->>>>>>
    record_dict = {
        "chr_name_order": []
    }

    # --------------------------------------------------->>>>>>
    # split file
    # --------------------------------------------------->>>>>>
    logging.info("Try to split mpmat file...")
    logging.info("Output dir is %s" % temp_dir)

    # open input mpmat file
    if (input_mpmat_filename[-3:] == ".gz") or (input_mpmat_filename[-5:] == ".gzip"):
        input_file = gzip.open(input_mpmat_filename, "r")
    else:
        input_file = open(input_mpmat_filename, "r")

    # set init 
    cur_chr_name = None
    cur_out_file = None

    for line in input_file:
        line_list = line.strip().split("\t")
        chr_name = line_list[0]

        if chr_name == "chr_name":
            continue

        if cur_chr_name is None:
            # log
            logging.info("Processing %s .mpmat file" % chr_name)

            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(temp_dir, temp_file_basename)

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            cur_out_file = open(temp_file_name, "w")
            cur_out_file.write(line)

        elif cur_chr_name == chr_name:
            cur_out_file.write(line)

        elif cur_chr_name != chr_name:
            cur_out_file.close()

            # log
            logging.info("Processing %s .mpmat file" % chr_name)

            # set next chr_name
            cur_chr_name = chr_name

            # make temp filename
            temp_file_basename = "temp_" + input_file_basename + "." + cur_chr_name + "." + "".join(
                random.sample(string.ascii_letters + string.digits, 16))
            temp_file_name = os.path.join(temp_dir, temp_file_basename)

            # record info into dict
            record_dict["chr_name_order"].append(cur_chr_name)
            record_dict[cur_chr_name] = temp_file_name

            # write
            cur_out_file = open(temp_file_name, "w")
            cur_out_file.write(line)

    # close all files
    try:
        cur_out_file.close()
        input_file.close()
    except:
        logging.error("Error occurs at close file step.")
        raise IOError("Error occurs at close file step.")

    # log 
    logging.info("Try to split .mpmat file. DONE!")

    # --------------------------------------------------->>>>>>
    # return
    # --------------------------------------------------->>>>>>
    return record_dict


#################################################################################
# FUN
# Info: calculating genome Poisson lambda background
#################################################################################
def _count_mismatch_num(align_mismatch_pairs, region_mut_type_list=["CT", "GA"]):
    """
    INPUT:
        <align_mismatch_pairs>
            list, generate from

    RETURN:
        <mismatch_count_dict>

    """
    total_mut_num = 0
    query_mut_num = 0
    other_mut_num = 0

    for ref_idx, align_idx, ref_base, align_base in align_mismatch_pairs:
        total_mut_num += 1
        mut_type = ref_base + align_base

        if mut_type in region_mut_type_list:
            query_mut_num += 1
        else:
            other_mut_num += 1

    return query_mut_num, other_mut_num, total_mut_num


def count_chr_bam_mut_count(
        input_bam_filename,
        select_chr_name,
        ref_genome_filename,
        query_mut_type_list=["CT", "GA"],
        query_mut_min_cutoff=1,
        query_mut_max_cutoff=16,
        total_mut_max_cutoff=20,
        other_mut_max_cutoff=16,
        log_verbose=3
):
    """
    INPUT:
        <input_bam_filename>
            str, bam filename

        <select_chr_name>
            str, like chr1, chr2, chr3 ...

        <ref_genome_filename>
            str, reference genome FASTA filename

        <query_mut_type_list>
            list, like ["CT","GA"]

        <query_mut_min_cutoff>
            int, if mutation number >= query_mut_min_cutoff in the mpmat region,
                will be counted as 'region_mut_count'.

        <query_mut_max_cutoff>
            int, larger than this will be marked.

        <total_mut_max_cutoff>
            int, larger than this will be marked.

        <other_mut_max_cutoff>
            int, larger than this will be marked.

    RETURN:
        <align_count_dict>
            dict, key and value like:
                {
                    "all_align_count": 0,
                    "region_mut_count": 0,
                    "region_non_mut_count": 0,
                    "total_high_mismatch_count": 0,
                    "other_high_mismatch_count": 0,
                    "query_high_mismatch_count": 0,
                    "all_filter_count": 0
                }

    """
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ---------------------------------------------------------->>>>>>>>>>
    # load genome and open BAM file
    # ---------------------------------------------------------->>>>>>>>>>
    input_bam = pysam.AlignmentFile(input_bam_filename, "rb")

    ref_genome_dict = load_reference_fasta_as_dict(ref_fasta_path=ref_genome_filename,
                                                   ref_name_list=[select_chr_name])
    # ---------------------------------------------------------->>>>>>>>>>
    # init var
    # ---------------------------------------------------------->>>>>>>>>>
    logging.info("Starting to count alignment mutation info on %s ..." % select_chr_name)

    # define count dict
    align_count_dict = {
        "all_align_count": 0,
        "region_mut_count": 0,
        "region_non_mut_count": 0,
        "total_high_mismatch_count": 0,
        "other_high_mismatch_count": 0,
        "query_high_mismatch_count": 0,
        "all_filter_count": 0
    }

    # ---------------------------------------------------------->>>>>>>>>>
    # iter align info
    # ---------------------------------------------------------->>>>>>>>>>
    for align in input_bam.fetch(reference=select_chr_name):
        # count total
        align_count_dict["all_align_count"] += 1

        # make sure MD state
        MD_state = True
        try:
            MD_str_tag = align.get_tag("MD")
        except:
            MD_state = False

        # get mismatch pairs
        if MD_state:
            align_mismatch_pairs = get_align_mismatch_pairs(align)
        else:
            align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)

        # analysis of mismatch pairs
        if align_mismatch_pairs is None:
            align_count_dict["region_non_mut_count"] += 1

        else:
            # load mut num
            align_query_mut_num, align_other_mut_num, align_total_mut_num = _count_mismatch_num(
                align_mismatch_pairs=align_mismatch_pairs,
                region_mut_type_list=query_mut_type_list)

            # filter high mismatch reads
            if align_total_mut_num >= total_mut_max_cutoff:
                align_count_dict["total_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            elif align_other_mut_num >= other_mut_max_cutoff:
                align_count_dict["other_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            elif align_query_mut_num >= query_mut_max_cutoff:
                align_count_dict["query_high_mismatch_count"] += 1
                align_count_dict["all_filter_count"] += 1

            else:
                if align_total_mut_num == 0:
                    align_count_dict["region_non_mut_count"] += 1

                elif align_query_mut_num >= query_mut_min_cutoff:
                    align_count_dict["region_mut_count"] += 1

    # close file
    input_bam.close()

    logging.info("Counting Poisson background on %s ... Done!" % select_chr_name)

    return align_count_dict


def multi_thread_chr_bam_mut_count(
        input_bam_filename,
        select_chr_name_list,
        ref_genome_filename,
        thread_num=1,
        query_mut_type_list=["CT", "GA"],
        query_mut_min_cutoff=1,
        query_mut_max_cutoff=16,
        total_mut_max_cutoff=20,
        other_mut_max_cutoff=16,
        log_verbose=3):
    """
    INPUT:
        <input_bam_filename>
            str, bam filename

        <select_chr_name_list>
            str, like ["chr1", "chr2", "chr3"]

        <ref_genome_filename>
            str, reference genome FASTA filename

        <thread_num>
            int, thread number to use

        <query_mut_type_list>
            list, like ["CT","GA"]

        <query_mut_min_cutoff>
            int, if mutation number >= query_mut_min_cutoff in the mpmat region,
                will be counted as 'region_mut_count'.

        <query_mut_max_cutoff>
            int, larger than this will be marked.

        <total_mut_max_cutoff>
            int, larger than this will be marked.

        <other_mut_max_cutoff>
            int, larger than this will be marked.

    RETURN:
        <meta_count_dict>
            meta_count_dict = {
                "meta_data" : {
                    "create_date" : time,
                    "input_BAM" : in_bam_filename,
                    "reference" : ref_genome_filename,
                    "threads" : thread_num,
                    "query_mutation_type" : val,
                    "query_mutation_min_cutoff" : val,
                    "query_mutation_max_cutoff" : val,
                    "other_mismatch_max_cutoff" : val,
                    "total_mismatch_max_cutoff" : val
                },
                "count_dict" : {
                    "chr1" : {
                        "all_align_count": 0,
                        "region_mut_count": 0,
                        "region_non_mut_count": 0,
                        "total_high_mismatch_count": 0,
                        "other_high_mismatch_count": 0,
                        "query_high_mismatch_count": 0,
                        "all_filter_count": 0
                    },
                    "chr2" : {
                        "all_align_count": 0,
                        "region_mut_count": 0,
                        "region_non_mut_count": 0,
                        "total_high_mismatch_count": 0,
                        "other_high_mismatch_count": 0,
                        "query_high_mismatch_count": 0,
                        "all_filter_count": 0
                    }....
                }
            }
    """

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    # log setting
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    logging.info("Starting to count genome Poisson lambada background...")

    pool = multiprocessing.Pool(processes=thread_num)

    # record info
    input_BAM_count_result = []

    for ref_name in select_chr_name_list:
        input_BAM_count_result.append(
            pool.apply_async(
                func=count_chr_bam_mut_count,
                args=(
                    input_bam_filename,
                    ref_name,
                    ref_genome_filename,
                    query_mut_type_list,
                    query_mut_min_cutoff,
                    query_mut_max_cutoff,
                    total_mut_max_cutoff,
                    other_mut_max_cutoff,
                    log_verbose,
                )
            )
        )

    pool.close()
    pool.join()

    # out pool result
    input_BAM_count_dict = {}
    for index, res in enumerate(input_BAM_count_result):
        run_res = res.get()
        input_BAM_count_dict[select_chr_name_list[index]] = run_res.copy()

    logging.info("Calculation done!")

    # make meta dict
    meta_count_dict = {
        "meta_data": {
            "create_date": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            "input_BAM": os.path.abspath(input_bam_filename),
            "reference": os.path.abspath(ref_genome_filename),
            "threads": thread_num,
            "query_mutation_type": ",".join(query_mut_type_list),
            "query_mutation_min_cutoff": query_mut_min_cutoff,
            "query_mutation_max_cutoff": query_mut_max_cutoff,
            "other_mismatch_max_cutoff": other_mut_max_cutoff,
            "total_mismatch_max_cutoff": total_mut_max_cutoff
        },
        "count_dict": input_BAM_count_dict.copy()
    }

    return meta_count_dict


def count_dict_sum_up(count_dict, key_ref_name_list="All", ignore_key_list=["total"]):
    """
    INPUT
        <count_dict>
            key = chr1, chr2, chr3 ...
            value = dict and key is

                  non_mut_align_count
                  all_align_count
                  total_high_mismatch_count
                  query_high_mismatch_count
                  mut_align_count
                  other_high_mismatch_coun
                  all_filter_count

        <key_ref_name_list>
            use those key to sum up and count

        <ignore_key_list>
            ignore

    RETURN
        Dict like <count_dict> and add 'total' key

    INFO
        Final changes, 2020-11-05
    """

    if key_ref_name_list == "All":
        key_ref_name_list = list(count_dict.keys())

    total_dict = {}
    inner_key_list = list(count_dict[key_ref_name_list[0]].keys())

    for inner_key in inner_key_list:
        total_dict[inner_key] = 0

    for ref_key in key_ref_name_list:
        if ref_key not in ignore_key_list:
            for ref_inner_key in inner_key_list:
                total_dict[ref_inner_key] += count_dict[ref_key][ref_inner_key]

    count_dict["total"] = total_dict.copy()

    return count_dict


def calculate_bg_scale_dict(meta_data_dict, to_large_state=False, scale_reads_count=None):
    """
    INPUT:
        <meta_data_dict>
            A dict generated from raw BAM file, and data structure like
                meta_count_dict = {
                    "meta_data" : {
                        "create_date" : ,
                        "ctrl_BAM" : ,
                        "treat_BAM" : ,
                        "reference" : ,
                        "threads" : ,
                        "query_mutation_type" : ,
                        "query_mutation_min_cutoff" : ,
                        "query_mutation_max_cutoff" : ,
                        "other_mismatch_max_cutoff" : ,
                        "total_mismatch_max_cutoff" : ,
                        "filter_log" :
                    },
                    "ctrl" : ,
                    "treat" : ,
                }

        <to_large_state>
            Default is False
                Means scale the larger sample up to the smaller sample.
                Usually, scaling down will bring down background noise, which is good for enrichment detection.

    RETURN:
        <scale_dict>


    INFO:
        Final changes, 2020-11-05
    """

    # init vars
    ctrl_scale_factor = 1
    treat_scale_factor = 1
    scale_dict = {"ctrl": {}, "treat": {}}

    select_key_list = [
        "all_align_count",
        "region_mut_count",
        "region_non_mut_count"
    ]

    if scale_reads_count == 0:
        raise ValueError("scale_reads_count can't be 0!")

    elif scale_reads_count is not None:
        ctrl_scale_factor = meta_data_dict["ctrl"]["total"]["all_align_count"] / 1.0 / int(scale_reads_count)
        treat_scale_factor = meta_data_dict["treat"]["total"]["all_align_count"] / 1.0 / int(scale_reads_count)

    for ref_name in meta_data_dict["ctrl"].keys():
        # make all count scale factor
        ctrl_all_count = meta_data_dict["ctrl"][ref_name]["all_align_count"]
        treat_all_count = meta_data_dict["treat"][ref_name]["all_align_count"]

        if to_large_state:
            all_scale_count = max(ctrl_all_count, treat_all_count)
        else:
            all_scale_count = min(ctrl_all_count, treat_all_count)

        if all_scale_count == 0:
            raise ValueError("%s all_align_count is 0!!!" % ref_name)

        ctrl_all_scale_factor = ctrl_all_count / 1.0 / all_scale_count
        treat_all_scale_factor = treat_all_count / 1.0 / all_scale_count

        # make dict
        if scale_dict["ctrl"].get(ref_name) is None:
            scale_dict["ctrl"][ref_name] = {}

        if scale_dict["treat"].get(ref_name) is None:
            scale_dict["treat"][ref_name] = {}

        for key in meta_data_dict["ctrl"][ref_name].keys():
            ctrl_count = meta_data_dict["ctrl"][ref_name][key]
            treat_count = meta_data_dict["treat"][ref_name][key]

            if scale_reads_count is None:
                if to_large_state:
                    scale_value = max(ctrl_count, treat_count)
                else:
                    scale_value = min(ctrl_count, treat_count)

                if scale_value == 0:
                    if key in select_key_list:
                        raise ValueError("%s %s scale_value is 0!!!" % (ref_name, key))

                    else:
                        sys.stderr.write("Warning %s %s scale_value is 0!!!\n" % (ref_name, key))
                        sys.stderr.write("Use all_align_count scale info instead." + "\n")
                        ctrl_scale_factor = ctrl_all_scale_factor
                        treat_scale_factor = treat_all_scale_factor

                else:
                    ctrl_scale_factor = ctrl_count / 1.0 / scale_value
                    treat_scale_factor = treat_count / 1.0 / scale_value

            # make new dict key
            if key in select_key_list:
                scale_dict_key = key + ".scale_factor"
                scale_dict["ctrl"][ref_name][scale_dict_key] = ctrl_scale_factor
                scale_dict["treat"][ref_name][scale_dict_key] = treat_scale_factor

    # return part
    return scale_dict


def calculate_effective_genome(ref_genome_filename, log_verbose=3):
    """
    INPUT:
        <ref_genome_filename>
            str, ref FASTA filename

    RETURN:
        <genome_base_count_dict>
            dict, key is chr1, chr2 ... and total
                  value is A,G,C,T,N count
    """
    # log setting
    logging.basicConfig(level=(4 - log_verbose) * 10,
                        format='%(levelname)-5s @ %(asctime)s: %(message)s ',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=sys.stderr,
                        filemode="w")

    logging.info("Starting to calculate genome effective length...")

    # init vars
    genome_base_count_dict = {"total": {"A": 0, "T": 0, "C": 0, "G": 0, "N": 0}}

    # open genome
    genome_fa = SeqIO.parse(handle=ref_genome_filename, format="fasta")

    for ref in genome_fa:
        logging.debug("Run on %s ..." % ref.id)

        chr_seq = ref.seq.upper()

        if genome_base_count_dict.get(ref.id) is None:
            genome_base_count_dict[ref.id] = {}

        # add into dict
        genome_base_count_dict[ref.id]["A"] = chr_seq.count("A")
        genome_base_count_dict[ref.id]["T"] = chr_seq.count("T")
        genome_base_count_dict[ref.id]["C"] = chr_seq.count("C")
        genome_base_count_dict[ref.id]["G"] = chr_seq.count("G")
        genome_base_count_dict[ref.id]["N"] = chr_seq.count("N")

        # add to total
        genome_base_count_dict["total"]["A"] += genome_base_count_dict[ref.id]["A"]
        genome_base_count_dict["total"]["T"] += genome_base_count_dict[ref.id]["T"]
        genome_base_count_dict["total"]["C"] += genome_base_count_dict[ref.id]["C"]
        genome_base_count_dict["total"]["G"] += genome_base_count_dict[ref.id]["G"]
        genome_base_count_dict["total"]["N"] += genome_base_count_dict[ref.id]["N"]

    # log
    logging.info("Calculating genome effective length, done!")

    return genome_base_count_dict


def back_all_normalization_scale_dict(meta_data_dict,
                                      genome_base_count_dict,
                                      seq_reads_length=150,
                                      norm_scale_reads_count=1e6):
    """
    INPUT:
        <meta_dict>
            dict, contain genome wide alignment value info, that is a combination of ctrl and treat result
                  generated from FUN multi_thread_chr_bam_mut_count

        <genome_base_count_dict>
            dict, generate from FUN calculate_effective_genome

        <norm_scale_reads_count>
            int, if set <scale_reads_count> is 1e6, the normalization value equals to CPM.

    RETURN:
        <scale_factor_dict>

        <normalize_scale_factor_dict>

        <genome_bg_dict>

    INFO
        Final changes, 2020-11-06
    """
    # -------------------------------------------------->>>>>>>>>>
    # init vars
    # -------------------------------------------------->>>>>>>>>>
    genome_bg_dict = {"ctrl": {
        "scale_mut_bg": {},
        "norm_scale_mut_bg": {},
        "scale_all_bg": {},
        "norm_scale_all_bg": {}
    }, "treat": {
        "scale_mut_bg": {},
        "norm_scale_mut_bg": {},
        "scale_all_bg": {},
        "norm_scale_all_bg": {}
    }}

    # -------------------------------------------------->>>>>>>>>>
    # make scale factor dict
    # -------------------------------------------------->>>>>>>>>>
    scale_factor_dict = calculate_bg_scale_dict(meta_data_dict, to_large_state=False)

    # -------------------------------------------------->>>>>>>>>>
    # make normalize scale factor
    # -------------------------------------------------->>>>>>>>>>
    normalize_scale_factor_dict = calculate_bg_scale_dict(meta_data_dict,
                                                          to_large_state=False,
                                                          scale_reads_count=int(norm_scale_reads_count))

    # -------------------------------------------------->>>>>>>>>>
    # genome background
    # -------------------------------------------------->>>>>>>>>>
    genome_effective_len = 0
    for dna_base in "AGCT":
        genome_effective_len += genome_base_count_dict["total"][dna_base]

    # genome level background
    ctrl_scale_factor = scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]
    treat_scale_factor = scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]

    ctrl_norm_scale_factor = normalize_scale_factor_dict["ctrl"]["total"]["all_align_count.scale_factor"]
    treat_norm_scale_factor = normalize_scale_factor_dict["treat"]["total"]["all_align_count.scale_factor"]

    # genome level mutation background
    ctrl_mut_bg_genome_lambda = meta_data_dict["ctrl"]["total"]["region_mut_count"] * seq_reads_length / 1.0 / ctrl_scale_factor / genome_effective_len
    treat_mut_bg_genome_lambda = meta_data_dict["treat"]["total"]["region_mut_count"] * seq_reads_length / 1.0 / treat_scale_factor / genome_effective_len

    ctrl_norm_mut_bg_genome_lambda = meta_data_dict["ctrl"]["total"]["region_mut_count"] * seq_reads_length / 1.0 / ctrl_norm_scale_factor / genome_effective_len
    treat_norm_mut_bg_genome_lambda = meta_data_dict["treat"]["total"]["region_mut_count"] * seq_reads_length / 1.0 / treat_norm_scale_factor / genome_effective_len

    genome_bg_dict["ctrl"]["scale_mut_bg"]["genome_bg"] = ctrl_mut_bg_genome_lambda
    genome_bg_dict["treat"]["scale_mut_bg"]["genome_bg"] = treat_mut_bg_genome_lambda

    genome_bg_dict["ctrl"]["norm_scale_mut_bg"]["genome_bg"] = ctrl_norm_mut_bg_genome_lambda
    genome_bg_dict["treat"]["norm_scale_mut_bg"]["genome_bg"] = treat_norm_mut_bg_genome_lambda

    # genome level all align background
    ctrl_all_bg_genome_lambda = meta_data_dict["ctrl"]["total"]["all_align_count"] * seq_reads_length / 1.0 / ctrl_scale_factor / genome_effective_len
    treat_all_bg_genome_lambda = meta_data_dict["treat"]["total"]["all_align_count"] * seq_reads_length / 1.0 / treat_scale_factor / genome_effective_len

    ctrl_norm_all_bg_genome_lambda = meta_data_dict["ctrl"]["total"]["all_align_count"] * seq_reads_length / 1.0 / ctrl_norm_scale_factor / genome_effective_len
    treat_norm_all_bg_genome_lambda = meta_data_dict["treat"]["total"]["all_align_count"] * seq_reads_length / 1.0 / treat_norm_scale_factor / genome_effective_len

    genome_bg_dict["ctrl"]["scale_all_bg"]["genome_bg"] = ctrl_all_bg_genome_lambda
    genome_bg_dict["treat"]["scale_all_bg"]["genome_bg"] = treat_all_bg_genome_lambda

    genome_bg_dict["ctrl"]["norm_scale_all_bg"]["genome_bg"] = ctrl_norm_all_bg_genome_lambda
    genome_bg_dict["treat"]["norm_scale_all_bg"]["genome_bg"] = treat_norm_all_bg_genome_lambda

    # -------------------------------------------------->>>>>>>>>>
    # chromosome background
    # -------------------------------------------------->>>>>>>>>>
    for ref_name in meta_data_dict["ctrl"].keys():
        # chr effective length
        chr_effective_len = 0
        for dna_base in "AGCT":
            chr_effective_len += genome_base_count_dict[ref_name][dna_base]

        # chr scale factor
        chr_ctrl_scale_factor = scale_factor_dict["ctrl"][ref_name]["all_align_count.scale_factor"]
        chr_treat_scale_factor = scale_factor_dict["treat"][ref_name]["all_align_count.scale_factor"]

        chr_ctrl_norm_scale_factor = normalize_scale_factor_dict["ctrl"][ref_name]["all_align_count.scale_factor"]
        chr_treat_norm_scale_factor = normalize_scale_factor_dict["treat"][ref_name]["all_align_count.scale_factor"]

        # chr mutation background
        ctrl_mut_bg_chr_lambda = meta_data_dict["ctrl"][ref_name]["region_mut_count"] * seq_reads_length / 1.0 / chr_ctrl_scale_factor / chr_effective_len
        treat_mut_bg_chr_lambda = meta_data_dict["treat"][ref_name]["region_mut_count"] * seq_reads_length / 1.0 / chr_treat_scale_factor / chr_effective_len

        norm_ctrl_mut_bg_chr_lambda = meta_data_dict["ctrl"][ref_name]["region_mut_count"] * seq_reads_length / 1.0 / chr_ctrl_norm_scale_factor / chr_effective_len
        norm_treat_mut_bg_chr_lambda = meta_data_dict["treat"][ref_name]["region_mut_count"] * seq_reads_length / 1.0 / chr_treat_norm_scale_factor / chr_effective_len

        # store in dict
        genome_bg_dict["ctrl"]["scale_mut_bg"][ref_name] = ctrl_mut_bg_chr_lambda
        genome_bg_dict["treat"]["scale_mut_bg"][ref_name] = treat_mut_bg_chr_lambda

        genome_bg_dict["ctrl"]["norm_scale_mut_bg"][ref_name] = norm_ctrl_mut_bg_chr_lambda
        genome_bg_dict["treat"]["norm_scale_mut_bg"][ref_name] = norm_treat_mut_bg_chr_lambda

        # chr all background
        ctrl_all_bg_chr_lambda = meta_data_dict["ctrl"][ref_name]["all_align_count"] * seq_reads_length / 1.0 / chr_ctrl_scale_factor / chr_effective_len
        treat_all_bg_chr_lambda = meta_data_dict["treat"][ref_name]["all_align_count"] * seq_reads_length / 1.0 / chr_treat_scale_factor / chr_effective_len

        norm_ctrl_all_bg_chr_lambda = meta_data_dict["ctrl"][ref_name]["all_align_count"] * seq_reads_length / 1.0 / chr_ctrl_norm_scale_factor / chr_effective_len
        norm_treat_all_bg_chr_lambda = meta_data_dict["treat"][ref_name]["all_align_count"] * seq_reads_length / 1.0 / chr_treat_norm_scale_factor / chr_effective_len

        # store in dict
        genome_bg_dict["ctrl"]["scale_all_bg"][ref_name] = ctrl_all_bg_chr_lambda
        genome_bg_dict["treat"]["scale_all_bg"][ref_name] = treat_all_bg_chr_lambda

        genome_bg_dict["ctrl"]["norm_scale_all_bg"][ref_name] = norm_ctrl_all_bg_chr_lambda
        genome_bg_dict["treat"]["norm_scale_all_bg"][ref_name] = norm_treat_all_bg_chr_lambda

    return scale_factor_dict, normalize_scale_factor_dict, genome_bg_dict
