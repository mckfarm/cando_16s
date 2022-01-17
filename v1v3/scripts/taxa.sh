#!/bin/bash
#SBATCH --job-name="taxa"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 05:00:00
#SBATCH -N 1
#SBATCH --mem=10G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=taxa.out
#SBATCH --error=taxa.err

module purge all
module load qiime2/2021.11

# midas classifier
qiime feature-classifier classify-sklearn --i-classifier /projects/b1052/shared/qiime/midas_4.8.1_classifier.qza \
--i-reads /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qza --o-classification /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_midas.qza

qiime metadata tabulate --m-input-file /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_midas.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_midas.qzv
