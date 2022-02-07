### ---------
# WSSC Initial 16s rRNA sequence analysis - MiDAS taxonomy assignments
# McKenna Farmer  
# January 2022

### packages and working directory ---------
setwd("~/Github/cando_16s/v1v3")

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
  tree="qiime/v1v3_rooted_tree.qza",
  taxonomy="qiime/v1v3_taxonomy_midas.qza",
  metadata = "qiime/v1v3_metadata.txt"
)

### summary stats ----------
# cleaning - remove mitochondria and chloroplasts
physeq <- subset_taxa(physeq, !Genus=="Mitochondria" & !Genus=="Chloroplast")

# relative abundance
rel <- transform_sample_counts(physeq, function(x) x*100/sum(x))
write.csv(rel@tax_table,"all_tax.csv")

topN <- 10
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
rel_nit_df$date <- dmy(rel_nit_df$date)
rel_phos_df$date <- dmy(rel_phos_df$date)
rel_10_df$date <- dmy(rel_10_df$date)

### correlation matrix PAO/GAO ------

# subset columns
rel_phos_df_subset_at1 <- rel_phos_df %>% filter(test=="AT1") %>% select("Abundance","date","Genus","test") %>% filter(Genus=="Ca_Accumulibacter" | 
                                                                                                                       Genus=="Tetrasphaera" |
                                                                                                                         Genus=="Ca_Competibacter")
rel_phos_df_subset_at23 <- rel_phos_df %>% filter(test=="AT2/AT3") %>% select("Abundance","date","Genus","test") %>% filter(Genus=="Ca_Accumulibacter" | 
                                                                                                                              Genus=="Tetrasphaera" |
                                                                                                                              Genus=="Ca_Competibacter")

# reshape to wide format for correlation function
rel_phos_df_subset_at1 <- spread(rel_phos_df_subset_at1, key = Genus, value = Abundance) 
rel_phos_df_subset_at23 <- spread(rel_phos_df_subset_at23, key = Genus, value = Abundance) 
write.csv(rbind(rel_phos_df_subset_at1,rel_phos_df_subset_at23),"./results_midas/phos_taxa.csv")


# save correlation matrix
cor_rel_phos_at1 <- cor(select(rel_phos_df_subset_at1,-c(date:test)))
cor_rel_phos_at23 <- cor(select(rel_phos_df_subset_at23,-c(date:test)))

# save significance matrix
cor_rel_phos_at1_pmat <- cor_pmat(select(rel_phos_df_subset_at1,-c(date:test)))
cor_rel_phos_at23_pmat <- cor_pmat(select(rel_phos_df_subset_at23,-c(date:test))) 

### plotting ---------------
# nitrifiers
ggplot(data=rel_nit_df,mapping=aes(x=date,y=Abundance,fill=Genus)) + 
  geom_bar(stat="identity") +   
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Isfahan2", 3)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") +
  labs(fill="Genus")

ggsave("./results_midas/relative_ab_nit.tiff",width=2500,height=1500,unit="px")

# PAOs/GAOs
ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") +
  theme_bw() +
  scale_fill_manual(values=met.brewer("Lakota", 5)) +
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") +
  xlab("Date (Month-Year)") +
  labs(fill="Genus")
ggsave("./results_midas/relative_ab_phos.tiff",width=3000,height=1500,unit="px")

ggplot(data=rel_phos_df,mapping=aes(x=date,y=Abundance,
  fill=factor(Genus,levels=c("Ca_Accumulibacter","Tetrasphaera","Ca_Obscuribacter","Dechloromonas","Halomonas","Microlunatus",
                             "Ca_Competibacter","Ca_Contendobacter","Defluviicoccus","Micropruina","Propionvibrio"))))
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 11)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") + 
  labs(fill="Genus")
ggsave("./results_midas/relative_ab_phos_nonputative.tiff",
       width=3000,height=1500,unit="px",scale=0.8)

# top 10 abundant genus
top10_plt <- ggplot(data=rel_10_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Hiroshige", 10)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("") + 
  labs(fill="Genus")

top10_at23 <- ggplot(data=rel_10_at23_df,mapping=aes(x=date,y=Abundance,fill=Genus)) +
  geom_bar(stat="identity") + 
  theme_bw() + 
  scale_fill_manual(values=met.brewer("Tiepolo", 10)) + 
  theme(axis.text.x = element_text(angle = 0)) +
  scale_x_date(date_labels="%m-%y") +
  ylab("Relative Abundance (%)") + 
  xlab("Date (Month-Year)") + 
  labs(fill="Genus") +
  ggtitle("Control - AT2/AT3")

plot_grid(top10_at1,top10_at23,nrow=2)
ggsave("./results_midas/top10.tiff",width=3000,height=1700,unit="px")

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

ggcorrplot(cor_rel_phos_at1, p.mat=cor_rel_phos_at1_pmat, 
           type="lower", lab="true", colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at1.tiff",width=1200,height=1000,unit="px")

ggcorrplot(cor_rel_phos_at23, p.mat=cor_rel_phos_at23_pmat, 
           type="lower", lab="true",colors=c("red","white","lightblue"))
ggsave("./results_midas/phos_at23.tiff",width=1200,height=1000,unit="px")
