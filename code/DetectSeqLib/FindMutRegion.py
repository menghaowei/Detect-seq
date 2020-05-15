# _*_ coding: UTF-8 _*_

from __future__ import division
from scipy import stats
from statsmodels.stats import proportion
from Bio import SeqIO

import pysam
import numpy as np
import sys
import os

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
VERSION_INFO = \
"""

<get_region_mut_count>
    Input a mpmat region and back mutant reads count, non-mutant reads count ...



"""
# Function List END ------------------------------------------------------------

#################################################################################
# load file and file check
#################################################################################
def load_reference_fasta_as_dict(ref_fasta_path, ref_name_list="All", show_load_id = True):
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
    
    ref_seq_dict = {}
    genome_fa =  SeqIO.parse(handle=ref_fasta_path, format="fasta")

    if ref_name_list != "All":
        ref_name_set = set(ref_name_list)

    for ref in genome_fa:
        if ref_name_list == "All":
            ref_seq_dict[ref.id] = ref.seq.upper()
            if show_load_id:
                sys.stderr.write("Loading...\t" + ref.id + "\n")

        elif ref.id in ref_name_set:
            ref_seq_dict[ref.id] = ref.seq.upper()
            
            if show_load_id:
                sys.stderr.write("Loading...\t" + ref.id + "\n")

            # remove already loaded seq
            ref_name_set.remove(ref.id)

            # load all info
            if len(ref_name_set) == 0:
                break

        else:
            pass
            # if show_load_id:
            #     sys.stderr.write("Skip...\t" + ref.id + "\n")
                        
    return(ref_seq_dict)


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


def get_BAM_ref_name(bam_filename, ref_count_cutoff = 0):
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
    
    return(ref_name_list)

#################################################################################
# Function part
#################################################################################
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
            
            ref_index is the same coordinate with UCSC genome browser
            
            When NM == 0, return None
    """
    # No mismatch 
    try:
        if align.get_tag("NM") == 0:
            return None
    except:
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
                    mismatch_pair_list.append([cur_index + align.reference_start, cur_index-1, ref_base, align.query_sequence[fix_index-1]])
                else:
                    return(None)
                
    return(mismatch_pair_list)


def get_No_MD_align_mismatch_pairs(align, ref_genome_dict):
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
            
            ref_index is the same coordinate with UCSC genome browser
            
            When NM == 0, return None
            
    """
    # No mismatch 
    try:
        if align.get_tag("NM") == 0:
            return None
    except:
        return None
    
    mismatch_pair_list = []
    for align_idx, ref_idx in align.get_aligned_pairs():
        if (align_idx != None) and (ref_idx != None):
            align_base = align.query_sequence[align_idx]
            ref_base = ref_genome_dict[align.reference_name][ref_idx]
            
            if align_base != ref_base:
                mismatch_pair_list.append([ 
                    ref_idx + 1,
                    align_idx,
                    ref_base,
                    align_base
                ])
                
    return(mismatch_pair_list)


def make_mismatch_count_dict(align_mismatch_pairs, start = None, end = None, region_mut_type_list = ["CT","GA"]):
    """
    INPUT:
            <start> and <end> 
                Genome coordinate same with UCSC genome browser system,
                if not set as None, means only mismatch in this region should be considered.
                default value is None.
    
    RETURN:
            <mismatch_count_dict>
                Key like CT means C to T mismatch 
            
    """
    count_dict = {}
    
    count_region_state = False
    if (start != None) and (end != None):
        count_region_state = True
        
    for ref_idx, align_idx, ref_base, align_base in align_mismatch_pairs:
        key = ref_base + align_base
        
        if count_dict.get(key) == None:
            count_dict[key] = 1
        else:
            count_dict[key] += 1
        
        if count_region_state:
            if key in region_mut_type_list:
                if start <= ref_idx <= end:
                    new_key = "region." + key
                    
                    if count_dict.get(new_key) == None:
                        count_dict[new_key] = 1
                    else:
                        count_dict[new_key] += 1 
            
    return(count_dict)


def count_mismatch_dict(count_dict, query_key_list = "All", other_key_list = "All", region_key_list = ["CT","GA"]):
    """
    HELP:
    
        <count_dict> 
            generate from Function <make_mismatch_count_dict>
        
        <query_key>
            default = "All", return all mismatch count 
            or only return sum of  select key count in query_key_list
        
        <other_key_list>
            default = "All"
            means return all mismatch count not in  <query_key_list>
        
        <region_key_list>
            default = ["CT","GA"]
            count CT GA mutation number in selected region
    
    RETURN
        [query_mis_count, other_mis_count, total_mis_count, region_mis_count]
    
    """
    query_mis_count = 0
    other_mis_count = 0
    total_mis_count = 0
    region_mis_count = 0
   
    # all key set
    all_key_set = set(count_dict.keys())

    # first count region info
    for key in region_key_list:
        region_key = "region." + key
        
        # count
        if count_dict.get(region_key) != None:
            region_mis_count += count_dict[region_key] 
        
        # rm region key
        if region_key in all_key_set:
            all_key_set.remove(region_key)
    
    # set All key
    if query_key_list == "All":
        query_key_list = all_key_set
    elif query_key_list == None:
        query_key_list = set()
      
    if other_key_list == "All":
        other_key_list == all_key_set - set(query_key_list)
    elif other_key_list == None:
        other_key_list = set()
    
    for key in all_key_set:
        count = count_dict[key]
        total_mis_count += count
        
        if key in query_key_list:
            query_mis_count += count
            
        elif key in other_key_list:
            other_mis_count += count
    
    return([query_mis_count, other_mis_count, total_mis_count, region_mis_count])


def get_mpmat_region_count(
    in_bam_obj, 
    mpmat_chr_name,
    mpmat_region_start,
    mpmat_region_end,
    ref_genome_dict,
    query_mismatch_type_list = ["CT","GA"], 
    other_mismatch_type_list = "All",
    region_mismatch_type_list = ["CT","GA"],
    region_query_mut_min_cutoff = 1,
    query_mut_max_cutoff = 16,
    total_mut_max_cutoff = 20,
    other_mut_max_cutoff = 16):
    """
    INPUT:
        <in_bam_obj>
            pysam.libcalignmentfile.AlignmentFile object

        <mpmat_chr_name>
        
        <mpmat_region_start>
        
        <mpmat_region_end>
    
        <region_mismatch_type_list>
            like ["CT","GA"], only those type of mutations in mpmat region considered as region mutation

        <region_query_mut_min_cutoff>
            no less than this cutoff consider as a region mutant read

        <query_mismatch_type_list>
            like ["CT","GA"], and set "All" means count every type of mismatch mutations

        <other_mismatch_type_list>
            like ["GT"], and set "All" means count every type of mismatch mutations, which not in <query_mismatch_type_list>
        
        <query_mut_max_cutoff>
            higher than this, considered as "query_high_mismatch"

        <total_mut_max_cutoff>
            higher than this, considered as "total_high_mismatch"
        
        <other_mut_max_cutoff>
            higher than this, considered as "other_high_mismatch"
        
            !!! The priority is  total_mut_max_cutoff > other_mut_max_cutoff > query_mut_max_cutoff

    
    RETURN:

        dict = {
            "all_align_count" : val,
            "region_mut_count" : 0,
            "region_non_mut_count" : 0,
            "all_filter_count" : val,
            "total_high_mismatch_count" : val,
            "other_high_mismatch_count" : val,
            "query_high_mismatch_count" : val
        }
 
    """
    
    # define count dict
    align_count_dict = {
        "all_align_count" : 0,
        "region_mut_count" : 0,
        "region_non_mut_count" : 0,
        "total_high_mismatch_count" : 0,
        "other_high_mismatch_count" : 0,
        "query_high_mismatch_count" : 0,
        "all_filter_count" : 0
    }

    for align in in_bam_obj.fetch(reference = mpmat_chr_name, start = mpmat_region_start - 1, end = mpmat_region_start + 1):    
        # count total
        align_count_dict["all_align_count"] +=1

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
            align_mismatch_pairs = get_No_MD_align_mismatch_pairs(align, ref_genome_dict)

        # count mismatch info
        if align_mismatch_pairs != None:

            count_dict = make_mismatch_count_dict(
                align_mismatch_pairs, 
                start = mpmat_region_start, 
                end = mpmat_region_end, 
                region_mut_type_list=["CT","GA"]
            )

            mut_count_res = count_mismatch_dict(
                count_dict, 
                query_key_list=["CT","GA"], 
                other_key_list="All", 
                region_key_list=["CT","GA"]
            )

            query_mut_count, other_mut_count, total_mut_count, region_mut_count = mut_count_res

            filter_reason = ""

            # filter part
            if total_mut_count >= total_mut_max_cutoff:
                align_count_dict["total_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Total High Mismatch Count"

            elif other_mut_count >= other_mut_max_cutoff:
                align_count_dict["other_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Other High Mismatch Count"

            elif query_mut_count >= query_mut_max_cutoff:
                align_count_dict["query_high_mismatch_count"] +=1
                align_count_dict["all_filter_count"] +=1
                filter_reason = "Query High Mismatch Count"

            else:
                if region_mut_count < region_query_mut_min_cutoff:
                    # non mut 
                    align_count_dict["region_non_mut_count"] +=1
                else:
                    # mut count 
                    align_count_dict["region_mut_count"] +=1

        else:
            # no mismatch align
            align_count_dict["region_non_mut_count"] +=1    


    # return part 
    return(align_count_dict)


#################################################################################
# Statistics Function part
#################################################################################
def _zstat_generic_parse(value, std_diff, alternative):
    '''
    INPUT:
        <value>
            Stats

        <std_diff>
            Std
        
        <alternative>
            Alternative string 
    
    RETURN:
        z-score, pvalue
    
    HELP:
        Copied and fixed from statsmodels.stats.weightstats
    '''

    zstat = value / std_diff

    if alternative in ['two-sided', '2-sided', '2s']:
        pvalue = stats.norm.sf(np.abs(zstat))*2

    elif alternative in ['larger', 'l']:
        pvalue = stats.norm.sf(zstat)

    elif alternative in ['smaller', 's']:
        pvalue = stats.norm.cdf(zstat)

    else:
        raise ValueError('invalid alternative')

    return zstat, pvalue
    

def poisson_test(lambda_1, lambda_2, method='score', alternative='two-sided'):
    '''
    INPUT:
    
        <lambda_1>, <lambda_2>
            Poisson parameter in case1 and case2
        
        <method>
            Support method:
                'wald': method W1A, wald test, variance based on separate estimates
                'score': method W2A, score test, variance based on estimate under Null
                'sqrt': W5A, based on variance stabilizing square root transformation

                'exact-cond': exact conditional test based on binomial distribution
                'cond-midp': midpoint-pvalue of exact conditional test
                    
        <alternative>
            'two-sided' 
                means H0: lambda_1 = lambda_2
                
            'larger' 
                means H0: lambda_1 > lambda_2
            
            'smaller'
                means H0: lambda_1 < lambda_2

    HELP:
    
        Reference paper:
            Gu, Ng, Tang, Schucany 2008: Testing the Ratio of Two Poisson Rates,
            Biometrical Journal 50 (2008) 2, 2008
        
        Raw code is copied from html and then fixed realted to DETECT-Seq project requirement
            https://stackoverflow.com/questions/33944914/implementation-of-e-test-for-poisson-in-python
    '''

    # calculate stat
    if method in ['score']:
        stat = (lambda_1 - lambda_2) / np.sqrt((lambda_1 + lambda_2))
        dist = 'normal'

    elif method in ['wald']:
        stat = (lambda_1 - lambda_2) / np.sqrt(lambda_1 + lambda_2)
        dist = 'normal'

    elif method in ['sqrt']:
        stat = 2 * (np.sqrt(lambda_1 + 3 / 1.0 / 8) - np.sqrt((lambda_2 + 3 / 1.0 / 8)))
        stat /= np.sqrt(2)
        dist = 'normal'

    elif method in ['exact-cond', 'cond-midp']:
        lambda_total = lambda_1 + lambda_2
        stat = None
        pvalue = proportion.binom_test(lambda_1, lambda_total, prop=0.5, alternative=alternative)
        
        if method in ['cond-midp']:
            # not inplace in case we still want binom pvalue
            pvalue = pvalue - 0.5 * stats.binom.pmf(lambda_1, lambda_total, 0.5)

        dist = 'binomial'

    # return part 
    if dist == 'normal':
        return _zstat_generic_parse(stat, 1, alternative)

    else:
        return stat, pvalue



