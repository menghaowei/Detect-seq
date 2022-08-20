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
        expand("bam.hisat3n/{sample}_hisat3n_hg38.bam",sample=SAMPLES),
        expand("bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.bam",sample=SAMPLES),
        expand("bam.hisat3n/{sample}_R{read_idx}_unmapped_and_LowerMAPQ20.fq.gz",sample=SAMPLES,read_idx=READ_IDX),
        expand("bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam",sample=SAMPLES),
        expand("bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam.bai",sample=SAMPLES),
        expand("bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.mapping_stats",sample=SAMPLES)


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# cutadapter
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule cutadapt:
    input:
        "raw.fastq/Detect-{sample}_R1.fastq.gz",
        "raw.fastq/Detect-{sample}_R2.fastq.gz"
    output:
        "fix.fastq/Detect-{sample}_R1_cutadapt.fq.gz",
        "fix.fastq/Detect-{sample}_R2_cutadapt.fq.gz"
    log:
        "fix.fastq/Detect-{sample}_cutadapt.log"
    shell:
        "cutadapt -j 20 --times 1  -e 0.1  -O 3  --quality-cutoff 25 \
        -m 55 -a AGATCGGAAGAGCACACGT \
        -A  AGATCGGAAGAGCGTCGTG \
        -o {output[0]} -p {output[1]} {input[0]} {input[1]} > {log} 2>&1 "    


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# hisat3n mapping 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule hisat3n_mapping:
    input:
        R1 = "fix.fastq/Detect-{sample}_R1_cutadapt.fq.gz",
        R2 = "fix.fastq/Detect-{sample}_R2_cutadapt.fq.gz"
    output:
        "bam.hisat3n/{sample}_hisat3n_hg38.bam"
    log:
        "bam.hisat3n/{sample}_hisat3n.log"
    params:
        unmapped = "bam.hisat3n/{sample}_hisat3n_hg38_unmapped.fq.gz",
        taginfo = "{sample}"
    shell:
        """
        hisat-3n -x {HISAT_3N_HG38_IDX} -1 {input.R1} -2 {input.R2} -p 20 \
        --sensitive --base-change C,T --unique-only --repeat-limit 1000 --no-spliced-alignment -X 700  \
        --un-conc-gz {params.unmapped} --summary-file {log} \
        --rg-id {params.taginfo} --rg "PL:ILLUMINA" --rg "ID:{params.taginfo}" --rg "SM:{params.taginfo}" | samtools view -hb > {output} 
        """

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# MAPQ select 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule high_MAPQ_hisat3n_bam:
    input:
        "bam.hisat3n/{sample}_hisat3n_hg38.bam"
    output:
        "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.bam"
    shell:
        "samtools view -@ 15 -h -b -q 20 -f 3 -F 256 -o {output} {input}"


rule lower_MAPQ_hisat3n_bam:
    input:
        "bam.hisat3n/{sample}_hisat3n_hg38.bam"
    output:
        "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.LowerMAPQ20.bam"
    shell:
        r"""samtools view -h -@ 4 {input} | awk '$1~"@" || $5 <= 20  {{print $0}}' |  samtools view -@ 4 -hb > {output}"""


rule sort_by_name:
    input:
        "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.LowerMAPQ20.bam"
    output:
        "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.LowerMAPQ20.SortName.bam"
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp -@ 15 -m 2G -n {input}"


rule lower_MAPQ20_make_FASTQ:
    input:
        "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.LowerMAPQ20.SortName.bam"
    output:
        fq1 = "bam.hisat3n/{sample}_hisat3n_hg38.LowerMAPQ20_R1.fq.gz",
        fq2 = "bam.hisat3n/{sample}_hisat3n_hg38.LowerMAPQ20_R2.fq.gz"
    log:
        "bam.hisat3n/{sample}_hisat3n_hg38.LowerMAPQ20.bamTofq.log"
    shell:
        "samtools fastq -@ 15 -0 /dev/null -s /dev/null -n -F 0x900 -1 {output.fq1} -2 {output.fq2} --reference {HG38_FA} {input} > {log} 2>&1"


rule merge_FASTQ:
    input:
        low_fq = "bam.hisat3n/{sample}_hisat3n_hg38.LowerMAPQ20_R{read_idx}.fq.gz",
        unmap_fq = "bam.hisat3n/{sample}_hisat3n_hg38_unmapped.fq.{read_idx}.gz"
    output:
        "bam.hisat3n/{sample}_R{read_idx}_unmapped_and_LowerMAPQ20.fq.gz"
    shell:
        "cat {input.low_fq} {input.unmap_fq} > {output}"


# ------------------------------------------------------------------------------------------>>>>>>>>>>
# lowerMAPQ and unmapped bwa realignment
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule bwa_mapping:
    input:
        "bam.hisat3n/{sample}_R1_unmapped_and_LowerMAPQ20.fq.gz",
        "bam.hisat3n/{sample}_R2_unmapped_and_LowerMAPQ20.fq.gz"
    output:
        "bam.hisat3n/{sample}_bwa_hg38_realign.bam"
    log:
        "bam.hisat3n/{sample}_bwa_hg38_realign.log"
    params:
        "'@RG\\tID:{sample}\\tPL:ILLUMINA\\tSM:{sample}'"
    shell:
        "bwa mem {BWA_HG38_IDX} {input[0]} {input[1]} -t 20 -M -R {params} 2>{log} | samtools view -hb > {output} "


rule high_MAPQ_bwa_bam:
    input:
        "bam.hisat3n/{sample}_bwa_hg38_realign.bam"
    output:
        "bam.hisat3n/{sample}_bwa_hg38_realign.MAPQ20.bam"
    shell:
        "samtools view -@ 15 -h -b -q 20 -f 3 -F 256 -o {output} {input}"

# ------------------------------------------------------------------------------------------>>>>>>>>>>
# merge BAM and make final output 
# ------------------------------------------------------------------------------------------>>>>>>>>>>
rule merge_BAM_file:
    input:
        hisat3n = "bam.hisat3n/{sample}_hisat3n_hg38.MAPQ20.bam",
        bwa = "bam.hisat3n/{sample}_bwa_hg38_realign.MAPQ20.bam"
    output:
        "bam.hisat3n/{sample}_hg38_merge.MAPQ20.bam"
    shell:
        "samtools cat -o {output} {input.hisat3n} {input.bwa}"


rule merge_BAM_sort_by_position:
    input:
        "bam.hisat3n/{sample}_hg38_merge.MAPQ20.bam"
    output:
        "bam.hisat3n/{sample}_hg38_merge_sort.MAPQ20.bam"
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp -@ 20 -m 2G {input}"


rule merge_BAM_rmdup:
    input:
        "bam.hisat3n/{sample}_hg38_merge_sort.MAPQ20.bam"
    output:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.WithClip.bam",
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.matrix"
    log:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.log"
    shell:
        "{JAVA} -Xms50g -Xmx50g -XX:ParallelGCThreads=20 -jar {PICARD} MarkDuplicates I={input} O={output[0]} M={output[1]} ASO=coordinate REMOVE_DUPLICATES=true 2>{log}"


rule merge_BAM_filter_clip:
    input:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.WithClip.bam"
    output:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam"
    shell:
        """samtools view -@ 4 -h {input} | samclip --ref {BWA_HG38_IDX} --max 3 --progress 0 | awk 'function abs(v) {{return v < 0 ? -v : v}} $1~"@" || ($7 == "=" && abs($9) <= 2500 ) {{print $0}}' | samtools view -@ 4 -hb > {output}"""


rule merge_BAM_sort_index:
    input:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam"
    output:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam.bai"
    shell:
        "samtools index -@ 10 {input} {output}"


rule merge_BAM_stats:
    input:
        bam = "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam",
        bam_index = "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.bam.bai"
    output:
        "bam.hisat3n/{sample}_hg38_merge_sort_rmdup.MAPQ20.mapping_stats"
    shell:
        "samtools stats -@ 10 --remove-overlaps --reference {BWA_HG38_IDX} {input.bam} > {output}"    

