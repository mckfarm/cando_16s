#!/bin/bash
#SBATCH --job-name="taxa_silva"
#SBATCH -A b1042
#SBATCH -p genomicsguest
#SBATCH -t 01:00:00
#SBATCH -N 1
#SBATCH --mem=30G
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=mckennafarmer2023@u.northwestern.edu
#SBATCH --output=taxa_silva.out
#SBATCH --error=taxa_silva.err


module purge all
module load qiime2/2021.11
# silva classifier
qiime feature-classifier classify-sklearn --i-classifier /projects/b1052/shared/qiime/silva-138-99-nb-classifier.qza \
--i-reads /projects/b1052/mckenna/cando_16s/qiime/v1v3_rep_seqs_dada2.qza --o-classification /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_silva.qza

qiime metadata tabulate --m-input-file /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_silva.qza \
--o-visualization /projects/b1052/mckenna/cando_16s/qiime/v1v3_taxonomy_silva.qzv
