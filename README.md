# Detect-seq tools

Final edition, 2020-05-15 By MENG Haowei


## Download and usage

All Python code in this repertory can be directly downloaded and used like:

```
python + sepecific_cmd 

```
## General analysis pipeline
![](./image/bioinfo_analysis_pipeline.png)

## Tool list

### From `BAM` file to `.mpmat` file

- `parse-mpileup.py` parse mpileup file into `.bmat` file
 
- `bmat2pmat.py` convet `.bmat` file into `.pmat` file

- `pmat-merge.py` merge `.pmat` file and generate `.mpmat` file

### Find significant off-target region

- `mpmat-select.py` select `.mpmat` file with setting cretirion, which can help you distinguish the real Detect-seq signals from background.

- `filter-mpmat-high-mismatch.py` Remove `.mpmat` regions with severe mismatch issue, commonly are highly repetative regions. 

- `find-significant-mpmat.py` Run a statistical test and report significant `.mpmat` file.

### Plot off-target alignment result table (`.art` file)

- `mpmat-to-art.py` Alignment sgRNA with off-target, generate `.art` file

- `plot-art.py` Plot `.art` file, output image like:

![art-plot](./image/art-plot.png)


## Requirement 

### Python version
- Python = 2.7x

### Necessary Python packages
- Biopython >= 2.2.4
- pysam >= 0.15
- pandas >= 0.24.2
- numpy >= 1.16.2
- matplotlib >= 1.74
- scikit-learn >= 0.20.3
