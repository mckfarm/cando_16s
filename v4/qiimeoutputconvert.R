### packages and working directory ---------
setwd("~/Github/cando_16s/v4")

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)

physeq <- qza_to_phyloseq(
  features="rarefied_table.qza",
  tree="rooted_tree.qza",
  taxonomy="taxonomy_midas.qza",
  metadata = "mar01_metadata.txt"
)

physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
write.csv(rel@tax_table,"all_tax.csv")
write.csv(rel@otu_table, "abundances.csv")


