# _*_ coding: UTF-8 _*_

#####################################################################################################
# Data 2022-07-21
# Author Howard MENG
# E-mail meng_howard@126.com
# Detect-seq BE4max-VEGFA DdCBE-ND6
#####################################################################################################
# --------------------------------------------------------------->>>>>>>
# pipeline
# --------------------------------------------------------------->>>>>>>
# 1. hisat3n alignment -> BAM1 + unmapped.FASTQ
# 2. BAM1 -> BAM1.MAPQ20 + BAM1.MAPQ_lower_20;
# 3. unmapped.FASTQ + BAM1.MAPQ_lower_20.FASTQ -> BAM2
# 4. BAM1.realignment + BAM1.MAPQ_lower_20.FASTQ -> final BAM2

# --------------------------------------------------------------->>>>>>>
# software
# --------------------------------------------------------------->>>>>>>
hisat3n = "/your_conda_path/envs/DetectSeq/bin/hisat3n"

PICARD = "/your_app_path/picard.jar"
JAVA = "/your_conda_path/bin/java"

# --------------------------------------------------------------->>>>>>>
# index and files
# --------------------------------------------------------------->>>>>>>
HISAT_3N_HG38_IDX = "/your_ref_path/reference/hisat3n_hg38_CT/hg38_only_chromosome.fa"

BWA_HG38_IDX = "/your_ref_path/reference/bwa_hg38/hg38_only_chromosome.fa"

HG38_FA = "/your_ref_path/reference/hisat3n_hg38_CT/hg38_only_chromosome.fa"

# rm INDEL info 
# if you don't have vcf file, please omit this
SNV_BED_INFO = "snv_info/293T-snv_info.list1.bed"
SNV_VCF_FILE = "snv_info/293T-snv_info.list2.vcf"


# --------------------------------------------------------------->>>>>>>
# vars
# --------------------------------------------------------------->>>>>>>
SAMPLES = [
    "293T-BE4max-mCherry-PD",
    "293T-BE4max-VEGFA-All-PD",
    "293T-DdCBE-GFP-PD",
    "293T-DdCBE-ND6-All-PD"
]

READ_IDX = ["1","2"]


rule all:
    input:
        expand("pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.pmat.gz",sample=SAMPLES),
        expand("pmat_and_mpmat/{sample}_hg38.MAPQ20_{ref_base_type}.pmat",sample=SAMPLES, ref_base_type=["C","G"]),
        expand("pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.CT.mpmat",sample=SAMPLES),
        expand("pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.GA.mpmat",sample=SAMPLES)
        

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# from bam to pmat  
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule bam_to_pmat:
    input:
        bam = "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam",
        bam_index = "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam.bai"
    output:
        "pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.pmat"
    log:
        "pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam2pmat.log"
    shell:
        """
        python ./bin/bam2pmat-V01.py -i {input.bam} -r {HG38_FA} -o {output} -p 20 \
        --out_format pmat --bed_like_format True --mut_type ALL --block_size 100000 \
        --cover_num_cutoff 0 --mut_num_cutoff 0 --mut_ratio_cutoff 0 --keep_temp_file False \
        --out_header False > {log}  2>&1
        """


rule pigz_pmat_file:
    input:
        "pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.pmat"
    output:
        "pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.pmat.gz"
    shell:
        "pigz -p 20 {input}"




rule select_pmat_C_and_G:
    input:
        "pmat_and_mpmat/{sample}_hg38_merge_sort_rmdup.MAPQ20.pmat.gz"
    output:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20_{ref_base_type}.pmat"
    wildcard_constraints:   
        ref_base_type = "C|G"
    params:
        "{ref_base_type}"
    shell:
        """zcat {input} | awk '$10 == "{params}" {{print $0}}' > {output}"""


rule pmat_merge_CT:
    input:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20_C.pmat"
    output:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.CT.mpmat"
    log:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.CT.mpmat.log"        
    shell:
        "python ./bin/pmat-merge-V05.py -i {input} -f C -t T -r {HG38_FA} -d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP {SNV_BED_INFO},{SNV_VCF_FILE} -o {output} 2> {log}" 



rule pmat_merge_GA:
    input:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20_G.pmat"
    output:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.GA.mpmat"
    log:
        "pmat_and_mpmat/{sample}_hg38.MAPQ20.merge_d50_D100.GA.mpmat.log"
    shell:
        "python ./bin/pmat-merge-V05.py -i {input} -f G -t A -r {HG38_FA} -d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP {SNV_BED_INFO},{SNV_VCF_FILE} -o {output} 2> {log}"         






