# _*_ coding: UTF-8 _*_

from __future__ import division
from Bio import SeqIO

import pysam
import sys
import os
import logging

# Version information START ----------------------------------------------------
VERSION_INFO = \
    """

Author: MENG Howard

Version-01:
  2019-12-27 find significant mpmat signal

E-Mail: meng_howard@126.com

"""
# Version information END ------------------------------------------------------


# Function List START ----------------------------------------------------------
COPY_RIGHT = \
"""
Belong to MENG Haowei @ YiLab in Peking University

"""
# Function List END ------------------------------------------------------------


#################################################################################
# load file and file check
#################################################################################
def load_reference_fasta_as_dict(ref_fasta_path, ref_name_list="All", log_verbose=30):
    """
    INPUT:
        <ref_fasta_path>
            Reference fasta file path
        
        <ref_name_list>
            If set All, load all seq info in reference, else only try to load seq_id in the list
    
    RETURN
        <ref_seq_dict>
            A dict, key is seq_id and value is sequence with  .upper()
        
        None
            If occur error, return None.
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
    # load genome as dict
    # ---------------------------------------------------------->>>>>>>>>>>>>>>>>>>>
    try:
        genome_fa = SeqIO.parse(handle=ref_fasta_path, format="fasta")
    except:
        raise IOError("Load file error! %s" % ref_fasta_path)

    # init var
    ref_seq_dict = {}
    ref_name_set = set(ref_name_list)

    logging.info("Starting to load the reference genome...")

    for ref in genome_fa:
        if ref_name_list == "All":
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

        elif ref.id in ref_name_list:
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

            # remove already loaded seq
            ref_name_set.remove(ref.id)

            # load all info
            if len(ref_name_set) == 0:
                break

    logging.info("Loading genome done!")

    return ref_seq_dict


def check_input_bam_file(input_bam_filename):
    """
    HELP
    
        1. check exist
        
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


def get_BAM_ref_name(bam_filename, ref_count_cutoff=0):
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

    return ref_name_list






