# CANDO+P reactor - 16s rRNA analysis
This document summarizes the workflow and scripts for analyzing 16s rRNA sequence data from the lab-scale CANDO+P reactor using the V4/V5 primer set. The QIIME workflow was performed on the Northwestern High Performance Computing Cluster, and commands used in a Quest interactive job are shown below. Data analysis was performed using R.


## QIIME2 workflow:
1) create manifest file
- created in Excel and saved to manifest_v4.csv

2) import paired end reads - these will be output as demultiplexed since they are imported with the manifest file

```
qiime tools import --input-path /projects/b1052/mckenna/cando_16s/v4/qiime/manifest_v4.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/b1052/mckenna/cando_16s/v4/qiime/v4.qza \
--type SampleData[PairedEndSequencesWithQuality]
```

3) trim primers

```
qiime cutadapt trim-paired --i-demultiplexed-sequences /projects/b1052/mckenna/cando_16s/v4/qiime/v4.qza  \
--p-front-f GTGYCAGCMGCCGCGGTAA \
--p-front-r CCGYCAATTYMTTTRAGTTT \
--p-error-rate 0.1 \
--p-overlap 3 \
--o-trimmed-sequences /projects/b1052/mckenna/cando_16s/v4/qiime/v4_trimmed.qza
```
error rate and overlap are default parameters - written out in the command for informational purposes

4) visualize read quality

```
qiime demux summarize --i-data /projects/b1052/mckenna/cando_16s/v4/qiime/v4.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/v4/qiime/v4_readquality.qzv

qiime demux summarize --i-data /projects/b1052/mckenna/cando_16s/v4/qiime/v4_trimmed.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/v4/qiime/v4_trimmed_readquality.qza
```
