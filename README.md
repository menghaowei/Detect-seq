# What is Detect-seq?
Detect-seq is short for ***dU-detection enabled by C to T transition during sequencing***. 

This method performs the genome-wide identification of CBE-induced off-targets in the cellular context. And Detect-seq is based on chemical labeling and enrichment of the direct editing products of CBE, which can trace the in vivo editing events in an unbiased manner.

This technique is developed by YiLab @ School of Life Sciences, Peking University, Beijing, China.

For full informations please click [www.detect-seq.com](http://www.detect-seq.com)


# What is Detect-seq tools?
Detect-seq tools are a collection of Python scripts, which is designed for the analysis of Detect-seq data. The tools can help to perform Detect-seq analysis including but are not limited to off-targets finding sgRNA alignment, and visualization of results.


# Requirement, Download and Usage

## 1. Requirement 

### Python version
- Python = 2.7x

### Necessary Python packages
- Biopython >= 2.2.4
- pysam >= 0.15
- pandas >= 0.24.2
- numpy >= 1.16.2
- matplotlib >= 1.74
- scikit-learn >= 0.20.3

## 2. Download and usage

All Python code in this repertory can be directly downloaded and used like:

```
python + sepecific_cmd 
```

## 3. Future 
We are going to distribute Detect-seq tools on `pip`, which will be coming soon! 


- - - - - - 


# The best practice of Detect-seq analysis

We provide some test data in `test` dir. So you can download the BAM files and find Detect-seq signal with the following steps.

## Contents
- [0. A general analysis pipeline](#0-a-general-analysis-pipeline)
- [1. From .BAM to .mpileup file](#1-from-bam-to-mpileup-file)
- [2. From .mpileup to .pmat file](#2-from-mpileup-to-pmat-file)
- [3. Merge .pmat file into .mpmat file](#3-merge-pmat-file-into-mpmat-file)
- [4. Run enrichment test with .mpmat file](#4-run-enrichment-test-with-mpmat-file)
- [5. Select signicicant regions and run sgRNA alignment](#5-select-signicicant-regions-and-run-sgrna-alignment)
- [6. Plot sgRNA alignment results](#6-plot-sgrna-alignment-results)

## 0. A general analysis pipeline
When you obtain the `BAM` files, you can follow this analysis pipeline to get your final off-target list and make a sgRNA alignment plot.

![](./image/bioinfo_analysis_pipeline.png)


## 1. From .BAM to .mpileup file

### 1.1 Requirement
1. FILE: sorted `BAM`
2. FILE: reference genome `FASTA` (Here we use `hg38.fa` as an example, and this ref file have to match your `BAM`)
3. CMD: `samtools` and version >= 1.9

### 1.2 Run code
You can generate `.mpileup` file from a sorted `.BAM` file with the following command:

```
samtools mpileup -q 20 -Q 20 --reference hg38.fa -o detect_seq.mpileup  detect_seq.sort.bam
```

The output file will be like:

```
chr1	1302588	G	1	^K.	e
chr1	1302589	T	1	.	j
chr1	1302590	G	1	.	o
chr1	1302591	T	1	.	o
chr1	1302592	G	1	.	o
chr1	1302593	T	1	.	s
chr1	1302594	C	1	.	s
chr1	1302595	C	1	.	s
chr1	1302596	A	1	.	s
chr1	1302597	T	1	.	s
```
The `.mpileup` format explain please check the HTML 
[mpileup explain](http://samtools.sourceforge.net/pileup.shtml)


## 2. From .mpileup to .pmat file
### 2.1 Requirement
1. FILE: `.mpileup`
2. CMD: `parse-mpileup-V04.py` (* This command support multiple threads.)
3. CMD: `bmat2pmat-V02.py`

### 2.2 From `.mpileup` to `.bmat` file
You can generate `.bmat` file from a `.mpileup` file with the following command:

```
parse-mpileup-V04.py -i detect_seq.mpileup -o detect_seq.bmat -p 1 -n 0
```

For help info, please run `python parse-mpileup-V04.py -h`:

```
python parse-mpileup-V04.py  -h
usage: parse-mpileup-V04.py [-h] -i INPUT [-o OUTPUT] [-p THREADS] [-n MUTNUM]
                        [--TempDir TEMPDIR]

convert mpileup file to info file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --Input INPUT
                        samtools mpileup format file
  -o OUTPUT, --Output OUTPUT
                        Output parsed file
  -p THREADS, --Threads THREADS
                        Multiple threads number, default=1
  -n MUTNUM, --MutNum MUTNUM
                        Only contain mutation info go to the output, set 0
                        mean output all site, default=0
  --TempDir TEMPDIR     Where to keep temp files, default is the same dir with
                        --Input
```

The `.bmat`format will looks like:

```
chr_name	chr_index	ref_base	A	G	C	T	del_count	insert_count	ambiguous_count	deletion	insertion	ambiguous	mut_num
chr1	1307682	C	0	0	83	0	0	0	0	.	.	.	0
chr1	1307683	G	0	85	0	0	0	0	0	.	.	.	0
chr1	1307684	C	0	0	83	0	0	0	0	.	.	.	0
chr1	1307685	T	0	0	0	84	0	0	0	.	.	.	0
chr1	1307686	G	0	82	0	0	0	0	0	.	.	.	0
chr1	1307687	A	77	0	0	0	1	0	0	C	.	.	0
chr1	1307688	C	0	0	50	31	0	0	0	.	.	.	31
```

### 2.3 select `C-base` and `G-base` info
As we all know, `BAM` files convert all sequence information into refrence strand, which also called foward strand (+). So if an CBE edit occur on the reverse strand (-), `BAM` file record that information as `G-to-A` rather than `C-to-T`. So it is necessary to split results into twe parts, a `C-base` part and the other is `G-base` part, represent foward strand (+) edits and reverse strand (-) edits respectively.

```
# select C-base info (possible foward strand edits)
awk '$3 == "C" {print $0}' detect_seq.bmat > detect_seq.C.bmat 

# select G-base info (possible reverse strand edits)
awk '$3 == "G" {print $0}' detect_seq.bmat > detect_seq.G.bmat 
```

### 2.4 From `.bmat` to `.pmat` file

`.bmat` and `.pmat` files are presented with a little difference. You can generate `.pmat` file from a `.bmat` file with the following command:

```
# For C-base info
bmat2pmat-V02.py -i detect_seq.C.bmat -o detect_seq.C.pmat --InHeader False --InLikeBED False --OutHeader True

# For G-base info
bmat2pmat-V02.py -i detect_seq.G.bmat -o detect_seq.G.pmat --InHeader False --InLikeBED False --OutHeader True
```

For help info, you can run `python bmat2pmat-V02.py -h`

```
usage: bmat2pmat-V02.py [-h] -i INPUT [-o OUTPUT] [-c COVERNUMCUTOFF]
                        [-m MUTNUMCUTOFF] [-r MUTRATIOCUTOFF] [-t MUTTYPE]
                        [--InHeader INHEADER] [--InLikeBED INLIKEBED]
                        [--OutHeader OUTHEADER]

convert bmat file to pmat file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --Input INPUT
                        Input bmat file
  -o OUTPUT, --Output OUTPUT
                        Output BED format file
  -c COVERNUMCUTOFF, --CoverNumCutoff COVERNUMCUTOFF
                        Site coverage number cutoff default=0
  -m MUTNUMCUTOFF, --MutNumCutoff MUTNUMCUTOFF
                        Site mutation number cutoff default=0
  -r MUTRATIOCUTOFF, --MutRatioCutoff MUTRATIOCUTOFF
                        Site mutation ratio cutoff default=0
  -t MUTTYPE, --MutType MUTTYPE
                        Select mutation type, ALL means no selection, can set
                        like CT, default=ALL
  --InHeader INHEADER   If contain header line in input file, default=True
  --InLikeBED INLIKEBED
                        If input bmat file looks like bed file, default=False
  --OutHeader OUTHEADER
                        If contain header line in output file, default=True
```

## 3. Merge .pmat file into .mpmat file
### 3.1 Requirement
1. FILE: reference genome `FASTA` (Here we use `hg38.fa` as an example, and this ref file have to match your `BAM`)
2. FILE: `.pmat` file
3. CMD: `pmat-merge-V04.py`

### 3.2 Explain
The `pmat` file only record sigle site information on each line, so we next try to search the tandem `C-to-T` (as `G-to-A` for reverse strand) pattern in the whole genome by `pmat-merge-V04.py` command. 

### 3.3 Run code

```
# For C-to-T .pmat merge
pmat-merge-V04.py -f C -t T -r hg38.fa \
-i detect_seq.C.pmat -o detect_seq.CT.mpmat \
-d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP SNP_info.vcf & 

# For G-to-A .pmat merge
pmat-merge-V04.py -f G -t A -r hg38.fa \
-i detect_seq.G.pmat -o detect_seq.GA.mpmat \
-d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP SNP_info.vcf & 
```


## 4. Run enrichment test with .mpmat file 
### 1.1 Requirement

### 1.2 Run code

## 5. Select signicicant regions and run sgRNA alignment
### 1.1 Requirement

### 1.2 Run code

## 6. Plot sgRNA alignment results
### 1.1 Requirement

### 1.2 Run code


## Tool list

### From `BAM` file to `.mpmat` file

- `parse-mpileup.py` parse mpileup file into `.bmat` file
 
- `bmat2pmat.py` convet `.bmat` file into `.pmat` file

- `pmat-merge.py` merge `.pmat` file and generate `.mpmat` file

### Find the significant off-target region

- `mpmat-select.py` select `.mpmat` file with the setting criterion, which can help you distinguish the real Detect-seq signals from the background.

- `filter-mpmat-high-mismatch.py` Remove `.mpmat` regions with severe mismatch issue, commonly are highly repetitive regions. 

- `find-significant-mpmat.py` Run a statistical test and report significant `.mpmat` file.

### Plot off-target alignment result table (`.art` file)

- `mpmat-to-art.py` Alignment sgRNA with off-target, generate `.art` file

- `plot-art.py` Plot `.art` file, output image like:

![art-plot](./image/art-plot.png)



