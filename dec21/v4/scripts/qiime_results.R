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

# data import --------
physeq <- qza_to_phyloseq(
  features="./qiime/rarefied_table.qza",
  tree="./qiime/rooted_tree.qza",
  taxonomy="./qiime/taxonomy_midas.qza",
  metadata = "./qiime/mar01_metadata.txt"
)

phases <- data.frame(x1=ymd("2021-05-26"),x2=ymd("2021-08-06"),
                     x3=ymd("2021-12-3"),x4=ymd("2022-04-04"),
                     x5=ymd("2022-05-24"))

# cleaning -------
# removing eukaryotes
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
# write.csv(rel@tax_table,"./results/all_tax.csv")

topN <- 20
rel_10_names <- sort(taxa_sums(rel), decreasing=TRUE)[1:topN]
rel_10 <- prune_taxa(names(rel_10_names), rel)

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
rel_phos_df$Genus <- factor(rel_phos_df$Genus,levels=c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas",
                                           "Ca_Competibacter","Micropruina"))
rel_10_df <- psmelt(rel_10)

# date formatting
rel_nit_df$date <- mdy(rel_nit_df$date)
rel_phos_df$date <- mdy(rel_phos_df$date)
rel_10_df$date <- mdy(rel_10_df$date)

# subset columns
rel_phos_df <- rel_phos_df %>% select("date","Abundance","Genus")
rel_10_df <- rel_10_df %>% select("date","Abundance","Genus")
# write.csv(rel_10_df,"./results/10_taxa.csv")


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
ggsave("./results/relative_ab_nit.tiff",width=2500,height=1500,unit="px")

# plotting PAOs/GAOs
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
ggsave("./results/relative_ab_phos.tiff",width=7,height=4,units="in")

# top 10 abundant genus
ggplot(data=rel_10_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Redon", 13)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %y") +
  ylab("Relative Abundance (%)") + 
  xlab("") + 
  labs(fill="Genus")

ggsave("./results/top10.tiff",width=3500,height=1500,unit="px",dpi=400)


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




