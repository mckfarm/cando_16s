# CANDO+P reactor - 16s rRNA analysis
This document summarizes the workflow and scripts for analyzing 16s rRNA sequence data from the lab-scale CANDO+P reactor using the V1/V3 primer set. The QIIME workflow was performed on the Northwestern High Performance Computing Cluster, and commands used in a Quest interactive job are shown below. Data analysis was performed using R.


## QIIME2 workflow:
1) create manifest file
- created in Excel and saved to manifest_v1v3.csv

2) import paired end reads - these will be output as demultiplexed since they are imported with the manifest file

```
qiime tools import --input-path /projects/b1052/mckenna/cando_16s/qiime/manifest_v1v3.csv \
--input-format PairedEndFastqManifestPhred33 \
--output-path /projects/b1052/mckenna/cando_16s/qiime/v1v3.qza \
--type SampleData[PairedEndSequencesWithQuality]
```

3) trim primers

```
qiime cutadapt trim-paired --i-demultiplexed-sequences /projects/b1052/mckenna/cando_16s/qiime/v1v3.qza  \
--p-front-f AGAGTTTGATCCTGGCTCAG \
--p-front-r ATTACCGCGGCTGCTGG \
--p-error-rate 0.1 \
--p-overlap 3 \
--o-trimmed-sequences /projects/b1052/mckenna/cando_16s/qiime/v1v3_trimmed.qza
```
error rate and overlap are default parameters - written out in the command for informational purposes

4) visualize read quality

```
qiime demux summarize --i-data /projects/b1052/mckenna/cando_16s/qiime/v1v3.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_readquality.qzv

qiime demux summarize --i-data /projects/b1052/mckenna/cando_16s/qiime/v1v3_trimmed.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_trimmed_readquality.qzv
```

5) [dada2.sh](https://github.com/mckfarm/cando_16s/blob/main/v1v3/scripts/dada2.sh)
- denoise and trim - trimming is based on read quality statistics from the previous step - keep sequences with average read quality of >20
- only using forward reads -> using the reverse reads ends up having dada2 filter out basically every sequence :/
- this command creates three files:
  - dada2 quality filtering table (stats)
  - data table of read info that can be coupled to metadata (table)
  - list of amplicon sequence variants that will be used for blast or other commands (rep_seqs)

6) create phylogenetic tree with mafft

```
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qza \
--o-alignment /projects/b1052/mckenna/cando_16s/qiime/v1v3_aligned_rep_seqs_dada2.qza \
--o-masked-alignment /projects/b1052/mckenna/cando_16s/qiime/v1v3_masked_aligned_rep_seqs_dada2.qza \
--o-tree /projects/b1052/mckenna/cando_16s/qiime/v1v3_unrooted_tree.qza \
--o-rooted-tree /projects/b1052/mckenna/cando_16s/qiime/v1v3_rooted_tree.qza
```

7) assign taxonomy with [taxa.sh](https://github.com/mckfarm/wssc/blob/main/scripts/taxa.sh)
- assigns taxa from Midas and Silva classifers
- also produces output qzv file for viewing

8) make alpha rarefaction curves
- measure of how diversity changes with sequencing depth

```
qiime diversity alpha-rarefaction \
--i-table /projects/b1052/mckenna/cando_16s/qiime/v1v3_table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/cando_16s/qiime/v1v3_rooted_tree.qz \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_rarefaction.qzv \
--p-max-depth 10000
```

9) rarefy samples
- picking a depth of 5000 based on rough estimate of plateau in faith_pd rarefaction curve
- rarefaction is important for equally comparing sequence data from different sampling dates and DNA extractions

```
qiime diversity core-metrics-phylogenetic \
--i-table /projects/b1052/mckenna/cando_16s/qiime/v1v3_table_dada2.qza \
--i-phylogeny /projects/b1052/mckenna/cando_16s/qiime/v1v3_rooted_tree.qza \
--p-sampling-depth 5000 \
--m-metadata-file /projects/b1052/mckenna/cando_16s/qiime/v1v3_metadata.txt \
--output-dir /projects/b1052/mckenna/cando_16s/qiime/core-metrics-results-5000
```

# Data analysis
[analysis.R](https://github.com/mckfarm/s2ebpr_16s/blob/main/analysis.R)
- Data analysis performed in R with phyloseq and other R functions
