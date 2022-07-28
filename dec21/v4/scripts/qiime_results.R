# packages and working directory ---------
setwd("~/Github/cando_16s/dec21/v4")

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lubridate)
library(MetBrewer)
library(ggcorrplot)
library(microbiome)

# data import --------
physeq <- qza_to_phyloseq(
  features="./qiime_outputs/table_dada2.qza",
  tree="./qiime_outputs/rooted_tree.qza",
  taxonomy="./qiime_outputs/taxonomy.qza",
  metadata = "./qiime_outputs/metadata.txt"
)


phases <- data.frame(x1=ymd("2021-05-26"),x2=ymd("2021-08-06"),
                     x3=ymd("2021-12-3"),x4=ymd("2022-04-04"),
                     x5=ymd("2022-05-24"))

# cleaning -------
# removing eukaryotes
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# rarefying
## this is an optional step where reads are randomly sampled from each sample
## the number of reads is equal to minimum number of reads 
physeq2 <- rarefy_even_depth(physeq,sample.size=min(sample_sums(physeq)),
                             rngseed=1)
rarefy_level <- min(sample_sums(physeq2))

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
write.csv(rel@tax_table,"./results/all_tax.csv")

# specific functional groups -------
## nitrifiers, PAOs/GAOs
rel_nit <- subset_taxa(rel,Genus=="Nitrotoga" | 
                         Genus=="Nitrospira" | 
                         Genus=="Nitrobacter" | 
                         Genus=="Nitrosomonas" | 
                         Genus=="Nitrosospira")
rel_phos <- subset_taxa(rel,Genus=="Ca_Accumulibacter" | 
                          Genus=="Tetrasphaera" |
                          Genus=="Dechloromonas" |
                          Genus=="Ca_Competibacter" |
                          Genus=="Micropruina")

rel_nit <- tax_glom(rel_nit,"Genus")
rel_phos <- tax_glom(rel_phos,"Genus")

# prep for plotting with ggplot - convert to dataframe format
rel_nit_df <- psmelt(rel_nit) 
rel_phos_df <- psmelt(rel_phos) 
rel_phos_df$Genus <- factor(rel_phos_df$Genus,
                            levels=c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas",
                                     "Ca_Competibacter","Micropruina"))

# date formatting
rel_nit_df$date <- mdy(rel_nit_df$date)
rel_phos_df$date <- mdy(rel_phos_df$date)


# plotting nitrifiers -----------
ggplot(data=rel_nit_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  geom_bar(stat="identity") +   
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Hiroshige", 3)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") +
  labs(fill="Genus")
ggsave("./results/relative_ab_nit.tiff",width=5,height=4,units="in")



# plotting PAOs/GAOs ------------
ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  geom_bar(stat="identity") + 
  theme_classic() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 6)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %y") +
  labs(y="Relative Abundance [%]",x="Date") + 
  geom_vline(xintercept=phases$x2,color="chocolate4") +
  ylim(0,40) +
  theme(legend.position="none")
ggsave("./results/relative_ab_phos.tiff",width=5,height=4,units="in")


# correlation matrix PAO/GAO ------
# reshape to wide format for correlation function
rel_phos_df_subset <- spread(rel_phos_df, key = Genus, value = Abundance) 
write.csv(rel_phos_df_subset,"./results/phos_taxa.csv")

# save correlation matrix
spear_rel_phos <- cor(select(rel_phos_df_subset,-c(date)),method="spearman")
pear_rel_phos <- cor(select(rel_phos_df_subset,-c(date)),method="pearson")

# save significance matrix
spear_rel_phos_pmat <- cor_pmat(select(rel_phos_df_subset,-c(date)),method="spearman",exact=F)
pear_rel_phos_pmat <- cor_pmat(select(rel_phos_df_subset,-c(date)),method="pearson")

# correlation plots
ggcorrplot(spear_rel_phos, p.mat=spear_rel_phos_pmat, 
           type="lower", lab="true", sig.level=0.05,
           insig="pch", colors=c("orangered3","white","lightblue"))
ggsave("./results/phos_corr.tiff",width=4,height=3,unit="in",scale=1.5)


# core microbiome -----------
# install microbiome package from BioConductor package manager
# https://www.bioconductor.org/
# https://microbiome.github.io/tutorials/

# genus level 
rel_glom <- tax_glom(rel,taxrank="Genus")

core_test <- core(subset_samples(rel_glom,location=="test"),
                  detection = 0.5, prevalence = 80/100)
core_con <- core(subset_samples(rel_glom,location=="control"),
                 detection = 0.5, prevalence = 80/100)
core_ras <- core(subset_samples(rel_glom,location=="ras"),
                 detection = 0.5, prevalence = 80/100)
plot_bar(core_test,"Sample","Abundance","Genus")


