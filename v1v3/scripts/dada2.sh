#!/bin/bash
#SBATCH --job-name="dada2"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 03:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=dada2.out
#SBATCH --error=dada2.err

module purge all
module load qiime2/2021.11

# dada2 command - only using the forward read :/
qiime dada2 denoise-single --verbose \
--i-demultiplexed-seqs /projects/b1052/mckenna/cando_16s/qiime/v1v3_trimmed.qza \
--p-trunc-len 250 \
--o-representative-sequences /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qza \
--o-table /projects/b1052/mckenna/cando_16s/qiime/v1v3_table_dada2.qza \
--o-denoising-stats /projects/b1052/mckenna/cando_16s/qiime/v1v3_stats_dada2.qza


# output visuals
qiime feature-table tabulate-seqs \
--i-data /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qzv

qiime metadata tabulate --m-input-file /projects/b1052/mckenna/cando_16s/qiime/v1v3_stats_dada2.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_stats_dada2.qzv
