# What is Detect-seq?
Detect-seq is short for ***dU-detection enabled by C to T transition during sequencing***. 

This method performs the genome-wide identification of CBE-induced off-targets in the cellular context. And Detect-seq is based on chemical labeling and enrichment of the direct editing products of CBE, which can trace the in vivo editing events in an unbiased manner.

This technique is developed by YiLab @ School of Life Sciences, Peking University, Beijing, China.

For full informations please click [www.detect-seq.com](http://www.detect-seq.com)


# What is Detect-seq tools?
Detect-seq tools are a collection of Python scripts, which is designed for the analysis of Detect-seq data. The tools can help to perform Detect-seq analysis including but are not limited to off-targets finding sgRNA alignment, and visualization of results.


# Requirement, Download, and Usage

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

We provide some test data in the `test` dir. So you can download the BAM files and find Detect-seq signal with the following steps.

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
chr1    1302588    G    1    ^K.    e
chr1    1302589    T    1    .    j
chr1    1302590    G    1    .    o
chr1    1302591    T    1    .    o
chr1    1302592    G    1    .    o
chr1    1302593    T    1    .    s
chr1    1302594    C    1    .    s
chr1    1302595    C    1    .    s
chr1    1302596    A    1    .    s
chr1    1302597    T    1    .    s
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
chr_name    chr_index    ref_base    A    G    C    T    del_count    insert_count    ambiguous_count    deletion    insertion    ambiguous    mut_num
chr1    1307682    C    0    0    83    0    0    0    0    .    .    .    0
chr1    1307683    G    0    85    0    0    0    0    0    .    .    .    0
chr1    1307684    C    0    0    83    0    0    0    0    .    .    .    0
chr1    1307685    T    0    0    0    84    0    0    0    .    .    .    0
chr1    1307686    G    0    82    0    0    0    0    0    .    .    .    0
chr1    1307687    A    77    0    0    0    1    0    0    C    .    .    0
chr1    1307688    C    0    0    50    31    0    0    0    .    .    .    31
```

### 2.3 select `C-base` and `G-base` info
As we all know, `BAM` files convert all sequence information into the reference strand, which also called forward strand (+). So if a CBE edit occurs on the reverse strand (-), `BAM` file record that information as `G-to-A` rather than `C-to-T`. So it is necessary to split results into two parts, a `C-base` part and the other is `G-base` part, represent forward strand (+) edits and reverse strand (-) edits respectively.

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
3. FILE: SNP annotation as `.vcf` format (not necessary)
4. CMD: `pmat-merge-V04.py`

### 3.2 Explain
The `pmat` file only records single site information on each line, so we next try to search the tandem `C-to-T` (as `G-to-A` for reverse strand) pattern in the whole genome by `pmat-merge-V04.py` command. 

And you can use the `--SNP` option to set SNP information to ignore SNP or SNV signals during your Detect-seq analysis. 

### 3.3 Run code

```
# For C-to-T .pmat merge
pmat-merge-V04.py -f C -t T -r hg38.fa --OutHeader False \
-i detect_seq.C.pmat -o detect_seq.CT.mpmat \
-d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP SNP_info.vcf & 

# For G-to-A .pmat merge
pmat-merge-V04.py -f G -t A -r hg38.fa --OutHeader False \
-i detect_seq.G.pmat -o detect_seq.GA.mpmat \
-d 50 -D 100 --NoMutNumCutoff 2 --OmitTandemNumCutoff 2 --SNP SNP_info.vcf & 
```

For help info, you can run `python pmat-merge-V04.py -h`

```
usage: pmat-merge-V04.py [-h] -i INPUT [-o OUTPUT] [-f FROMBASE] [-t TOBASE]
                         -r REFERENCE [-d MAXSITEDISTANCE]
                         [-D MAXREGIONDISTANCE]
                         [--NoMutNumCutoff NOMUTNUMCUTOFF]
                         [--OmitTandemNumCutoff OMITTANDEMNUMCUTOFF]
                         [--SNP SNP] [--OutHeader OUTHEADER]
                         [--InHeader INHEADER]

merge pmat file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --Input INPUT
                        Input bmat file
  -o OUTPUT, --Output OUTPUT
                        Output BED format file
  -f FROMBASE, --FromBase FROMBASE
                        Ref base, accept A,G,C,T default=C
  -t TOBASE, --ToBase TOBASE
                        Mut base, accept A,G,C,T default=T
  -r REFERENCE, --reference REFERENCE
                        Reference genome fasta file
  -d MAXSITEDISTANCE, --MaxSiteDistance MAXSITEDISTANCE
                        Max distance between two sites in one region,
                        default=50
  -D MAXREGIONDISTANCE, --MaxRegionDistance MAXREGIONDISTANCE
                        Max length of a mutation region, default=100
  --NoMutNumCutoff NOMUTNUMCUTOFF
                        The number of site without mutation --ToBase signal in
                        a mutation region, default=2
  --OmitTandemNumCutoff OMITTANDEMNUMCUTOFF
                        The omit tande site cutoff, default=2
  --SNP SNP             SNP file with vcf or bed format, if use multiple file,
                        use ',' to separate, default=None
  --OutHeader OUTHEADER
                        If contain header line in output file, default=True
  --InHeader INHEADER   If contain header line in input file, default=True
```


## 4. Run enrichment test with .mpmat file 
### 4.1 Requirement
1. FILE: `.mpmat`
2. FILE: Genome background file in `.json` format
3. CMD: `calculate-mut-stats-V02.py` (* This command support multiple threads.)
4. CMD: `find-significant-mpmat-V02.py`

### 4.2 Merge `C-to-T` and `G-to-A` `.mpmat` files
In this step, we should merge two strands `.mpmat` file together.

```
# merge file 
cat detect_seq.CT.mpmat detect_seq.GA.mpmat > detect_seq.merge.mpmat

# sort file 
bedtools sort -g hg38.fa.fai -i detect_seq.merge.mpmat  > detect_seq.merge.sort.mpmat
```

### 4.3 Count genome background and generate `.json` files
Next to calculate genome background before running a statistical test. 

```
# ctrl sample
calculate-mut-stats-V02.py -i Full.Ctrl.bam -r hg38.fa -p 1 -o background_ctrl.json

# treat (PD) sample
calculate-mut-stats-V02.py -i Full.Treat.bam -r hg38.fa -p 1 -o background_treat.json
```

This step has to provide a full-size bam with all genome mapping reads into a script, so here we just show a pretend demo. And you can find `background_ctrl.json` and `background_treat.json` files in the `test` dir.

### 4.4 Run Poisson test
Run Poisson test with `find-significant-mpmat-V02.py` and this idea refers to [MACS2](https://github.com/taoliu/MACS)

```
find-significant-mpmat-V02.py \
-i detect_seq.merge.sort.mpmat \
-c detect_seq.ctrl.sort.bam \
-t detect_seq.sort.bam \
-r hg38.fa \
-m background_ctrl.json \
-n background_treat.json \
-g hg38.json \
-o detect_seq.StatsTest.table \
--region_mutation_min_cutoff 2 \
--query_mutation_type CT,GA \
--query_mutation_min_cutoff 2 \
--query_mutation_max_cutoff 18 \
--other_mutation_max_cutoff 12 \
--total_mutation_max_cutoff 26
```

Then you will have a test result like:

```
chr_name    region_start    region_end    ctrl_count    treat_count    ctrl_mut_count    treat_mut_countctrl_count.norm    treat_count.norm    ctrl_mut_count.norm    treat_mut_count.norm    log2_FC    log2_FC_mut    p_value    FDR
chr1    1302672    1302672    4    5    0    0    0.0    0.0    2.88415875157    2.54780515605    NA    -0.178895625038    0.9433256859064737    0.9433256859064737
chr1    1307688    1307704    0    108    0    105    0.0    53.5039082771    0.0    55.0325913707    5.74157237418    5.78221435867    1.4956224254049805e-16    6.89481938111696e-14
chr1    93448684    93448686    7    85    0    1    0.0    0.50956103121    5.04727781524    43.3126876529    -0.972673143489    3.10121229415    0.8426930101168323    0.9433256859064737
chr1    93448687    93448705    6    85    1    75    0.721039687892    38.2170773408    4.32623812735    43.3126876529    5.72799497057    3.32360471549    6.577261620968718e-12    2.2740882054499343e-09
chr1    109629399    109629422    6    157    0    144    0.0    73.3767884943    4.32623812735    80.0010819    6.19725185795    4.20883452825    1.4200100472603131e-22    8.416602408690056e-20
chr1    109629400    109629406    7    157    0    3    0.0    1.52868309363    5.04727781524    80.0010819    0.612289357232    3.98644210691    0.5724525302864157    0.9433256859064737
chr2    396529    396546    5    36    0    31    0.0    15.7963919675    3.60519843946    18.3441971236    3.9815231669    2.34717318663    2.2714678759049305e-05    0.004712160108564779
chr6    68640125    68640148    6    171    0    136    0.0    69.3003002446    4.32623812735    87.1349363369    6.11478969776    4.33206629424    1.4101914795343944e-21    7.313605560735253e-19
chr11    69704140    69704158    0    47    0    42    0.0    21.4015633108    0.0    23.9493684669    4.41964427929    4.58191570819    4.333639994373281e-07    9.463301229818284e-05
chr12    105659914    105659927    8    304    0    263    0.0    134.014551208    5.76831750314    154.906553488    7.0662458458    4.74710379352    3.1570106740925235e-41    5.849243543839085e-38
chr17    66303049    66303069    0    53    0    44    0.0    22.4206853732    0.0    27.0067346541    4.48675847515    4.75524731107    2.6013830159213885e-07    5.9961878516988e-05
```

The output table column explain is:

- `ctrl_count` reads count in control `BAM` file
- `treat_count` reads count in treat `BAM` file
- `ctrl_mut_count` reads count in control `BAM` file with tandem mutation signals
- `treat_mut_count` reads count in treat `BAM` file with tandem mutation signals
- `ctrl_count.norm` nomalized reads count in control `BAM` file
- `treat_count.norm` nomalized reads count in treat `BAM` file 
- `ctrl_mut_count.norm` nomalized mutation reads count in control `BAM` file
- `treat_mut_count.norm` nomalized mutation reads count in treat `BAM` file
- `log2_FC` calculate as `log2_FC = log2(treat_count.norm / ctrl_count.norm)`
- `log2_FC_mut` calculate as `log2_FC_mut = log2(treat_mut_count.norm / ctrl_mut_count.norm)`
- `p_value` Poisson test pvalue
- `FDR` Adjust pvalue with Benjamini-Hochberg method.

## 5. Select signicicant regions and run sgRNA alignment
### 5.1 Requirement
1. FILE: `detect_seq.StatsTest.table`
2. CMD: `mpmat-to-art-V04.py`

### 5.2 Select significant region
Related to our experience, the regions with the following criterion can be considered as enriched off-targets region:

1. `FDR` < 0.05;
2. `log2_FC_mut` > 1;
3. `ctrl_mut_count` < 3 OR `ctrl_mut_count.norm` < 3;

After this select step, we can obtain `detect_seq.StatsTest.Sign.table`.

### 5.2 Select `.mpmat` file 
```
bedtools intersect -a detect_seq.merge.sort.mpmat -b detect_seq.StatsTest.Sign.table -wa > detect_seq.merge.sort.Sign.mpmat
```

### 5.3 Run sgRNA alignment
You can run sgRNA alignment as the following command. Only one thing you should take care of, the `--sgRNA` sequence has to include `PAM` info.

```
mpmat-to-art-V04.py \
-i detect_seq.merge.sort.Sign.mpmat -r hg38.fa --sgRNA GACCCCCTCCACCCCGCCTCCGG > detect_seq.merge.sort.Sign.art
```
For more alignment options you can run `python mpmat-to-art-V04.py -h`

```
usage: mpmat-to-art-V04.py [-h] -i MPMAT_TABLE --sgRNA SGRNA
                           [-o OUT_ALIGNMENT_RESULT] -r REFERENCE
                           [-s SAMPLE_NAME] [--extend_method EXTEND_METHOD]
                           [--align_dist_to_signal ALIGN_DIST_TO_SIGNAL]
                           [--bedtools_path BEDTOOLS_PATH]
                           [--align_method ALIGN_METHOD]
                           [--align_settings ALIGN_SETTINGS]
                           [--align_min_score ALIGN_MIN_SCORE]
                           [--input_header INPUT_HEADER]
                           [--input_sep INPUT_SEP]
                           [--more_colname MORE_COLNAME]

From .mpmat to alignment result table (.art) file

optional arguments:
  -h, --help            show this help message and exit
  -i MPMAT_TABLE, --mpmat_table MPMAT_TABLE
                        .mpmat table file, generated from <pmat-merge.py> or
                        <mpmat-select.py>
  --sgRNA SGRNA         sgRNA sequence with PAM (NGG/NAG) sequence
  -o OUT_ALIGNMENT_RESULT, --out_alignment_result OUT_ALIGNMENT_RESULT
                        Output alignment result table (.art) filename,
                        default=Stdout
  -r REFERENCE, --reference REFERENCE
                        Reference genome fasta file
  -s SAMPLE_NAME, --sample_name SAMPLE_NAME
                        Sample name of this .mpmat, default=run_mpmat
  --extend_method EXTEND_METHOD
                        Select and extend region of mpmat file to get FASTA
                        file, <region> <upstream_site> <downstream_site>
                        <highest_site>, default=highest_site
  --align_dist_to_signal ALIGN_DIST_TO_SIGNAL
                        If distance too far between mutation signal and
                        alignment, consider as there is no appropriate
                        alignment, default=20
  --bedtools_path BEDTOOLS_PATH
                        Software <bedtools> PATH, default considered
                        <bedtools> is already included in current PATH
                        environment.
  --align_method ALIGN_METHOD
                        Provide two methods, <PAM_first> <fair_align>,
                        'PAM_first' treat PAM region with a higher weight,
                        default=fair_align
  --align_settings ALIGN_SETTINGS
                        Set <align_match_score> <align_mismatch_score>
                        <align_gap_open_score> <align_gap_extension_score>,
                        default=5,-4,-24,-8
  --align_min_score ALIGN_MIN_SCORE
                        If alignment score lower than this, consider as no
                        appropriate alignment, default=15
  --input_header INPUT_HEADER
                        If .mpmat file contain header, default=False
  --input_sep INPUT_SEP
                        default=\t
  --more_colname MORE_COLNAME
                        More info you want to include in output table, so you
                        can write the column names like col1,col2...
                        default=None
```

## 6. Plot sgRNA alignment results
### 6.1 Requirement
1. FILE: `detect_seq.merge.sort.Sign.art`
2. CMD: `plot-art-V01.py`

### 6.2 Run code
```
plot-art-V01.py -i detect_seq.merge.sort.Sign.art --sgRNA GACCCCCTCCACCCCGCCTCCGG --out_figure_format png -o out.art.png
```

After a few seconds, you can obtain an image with a match and mismatch information between the sgRNA sequence and off-target region.

![art-plot](./test/out.art.png)