### ---------
# WSSC Initial 16s rRNA sequence analysis - MiDAS taxonomy assignments
# McKenna Farmer

### packages and working directory ---------
setwd("~/Github/cando_16s/v4v5/")

library(qiime2R)
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(lubridate)
library(MetBrewer)
library(ggcorrplot)

### data import -----------
physeq <- qza_to_phyloseq(
  features="qiime/rarefied_table.qza",
  tree="qiime/rooted_tree.qza",
  taxonomy="qiime/taxonomy_midas.qza",
  metadata = "qiime/mar01_metadata.txt"
)

### summary stats ----------
# cleaning - remove mitochondria and chloroplasts
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
write.csv(rel@tax_table,"all_tax.csv")

topN <- 20
rel_10_names <- sort(taxa_sums(rel), decreasing=TRUE)[1:topN]
rel_10 <- prune_taxa(names(rel_10_names), rel)


# specific genus of interest
## nitrifiers, PAOs/GAOs
rel_nit <- subset_taxa(rel,Genus=="Nitrotoga" | 
                           Genus=="Nitrospira" | 
                           Genus=="Nitrobacter" | 
                           Genus=="Nitrosomonas" | 
                           Genus=="Nitrosospira")
rel_phos <- subset_taxa(rel,Genus=="Ca_Accumulibacter" | 
                            Genus=="Tetrasphaera" |
                            Genus=="Ca_Obscuribacter" |
                            Genus=="Dechloromonas" |
                            Genus=="Halomonas" |
                            Genus=="Microlunatus" |
                            Genus=="Ca_Competibacter" |
                            Genus=="Ca_Contendobacter" |
                            Genus=="Defluviicoccus"|
                            Genus=="Micropruina"|
                            Genus=="Propionvibrio")

rel_nit <- tax_glom(rel_nit,"Genus")
rel_phos <- tax_glom(rel_phos,"Genus")

# prep for plotting with ggplot - convert to dataframe format
rel_nit_df <- psmelt(rel_nit) 
rel_phos_df <- psmelt(rel_phos)
rel_10_df <- psmelt(rel_10)

# date formatting
rel_nit_df$date <- mdy(rel_nit_df$date)
rel_phos_df$date <- mdy(rel_phos_df$date)
rel_10_df$date <- mdy(rel_10_df$date)

### correlation matrix PAO/GAO ------

# subset columns
rel_phos_df <- rel_phos_df %>% select("Abundance","date","Genus")

# reshape to wide format for correlation function
rel_phos_df_subset <- spread(rel_phos_df, key = Genus, value = Abundance) 
write.csv(rel_phos_df_subset,"phos_taxa.csv")

# save correlation matrix
cor_rel_phos <- cor(select(rel_phos_df_subset,-c(date)))

# save significance matrix
cor_rel_phos_pmat <- cor_pmat(select(rel_phos_df_subset,-c(date)))

### plotting ---------------
# nitrifiers
ggplot(data=rel_nit_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  geom_bar(stat="identity") +   
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Hiroshige", 3)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") +
  labs(fill="Genus")

ggsave("./results_midas/relative_ab_nit.tiff",width=2500,height=1500,unit="px")

# PAOs/GAOs
ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,
  fill=factor(Genus,levels=c("Ca_Accumulibacter","Tetrasphaera","Ca_Competibacter")))) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_manual(values=met.brewer("Hiroshige", 3)) +
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") +
  xlab("Date (Month-Year)") +
  labs(fill="Genus")

ggsave("./results_midas/relative_ab_phos.tiff",width=3000,height=1500,unit="px")

ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,
  fill=factor(Genus,levels=c("Ca_Accumulibacter","Tetrasphaera","Ca_Obscuribacter","Dechloromonas","Halomonas","Microlunatus",
                             "Ca_Competibacter","Ca_Contendobacter","Defluviicoccus","Micropruina","Propionvibrio")))) + 
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 11)) + 
  theme(axis.text.x = element_text(angle = 0),legend.text = element_text(face = "italic")) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date") + 
  labs(fill="Genus")

ggsave("relative_ab_phos.tiff",
       width=3000,height=2000,unit="px",dpi=400)

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

ggsave("top10.tiff",width=3500,height=1500,unit="px",dpi=400)

# shannon diversity
shannon <- plot_richness(physeq, x="test", measures=c("Shannon"))
shannon_df <- as.data.frame(shannon$data)
shannon_df$date <- dmy(shannon_df$date)

# used to format date - 
# https://stackoverflow.com/questions/5024798/how-can-a-color-gradient-based-on-date-be-applied-to-a-ggplot2-scatter-plot
# https://stackoverflow.com/questions/21311489/scatter-plot-with-ggplot2-colored-by-dates
as.Date_origin <- function(x){
  as.Date(x, origin = "1970-01-01")
}

ggplot(data=shannon_df) + 
  geom_point(aes(x=test, y=value, color=as.integer(date))) +
  geom_boxplot(aes(x=test, y=value), alpha=0.1) +
  theme_classic() + 
  ylab("Shannon Diversity Index") +
  xlab("Test battery") + 
  scale_colour_gradientn(colors=met.brewer("Lakota", 3), 
                        limits=as.integer(as.Date(c("2021-06-25","2021-12-27"))),
                        labels=as.Date_origin) + 
  labs(color="Date") +
  ylim(4.3,5.3)

ggsave("./results_midas/shannon.tiff",width=1200,height=1000,unit="px")


# correlation plots

ggcorrplot(cor_rel_phos, p.mat=cor_rel_phos_pmat, 
           type="lower", lab="true", insig="blank",colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at1.tiff",width=1200,height=1000,unit="px",scale=1.5)

ggcorrplot(cor_rel_phos_at23, p.mat=cor_rel_phos_at23_pmat, 
           type="lower", lab="true",insig="blank",colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at23.tiff",width=1200,height=1000,unit="px",scale=1.5)
