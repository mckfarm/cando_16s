# CANDO+P reactor - 16s rRNA analysis
This repo summarizes the workflow and scripts for analyzing 16s rRNA sequence data from the lab-scale CANDO+P reactor. DNA was extracted from reactor biomass and sequenced at the Rush Genomics Core and Microbiome Core Facility.

As of March 25 2022, GitHub is organized by sequencing date:

## December 21 summary
16s rRNA PCR was performed with two different primer sets to compare the taxonomic results. The V4/V5 workflow uses the primers 515F/926R as recommended in (Parada et al 2016)[https://sfamjournals.onlinelibrary.wiley.com/doi/full/10.1111/1462-2920.13023]. This is what I will use going forward for community analysis.

V4-V5 515F	GTGYCAGCMGCCGCGGTAA
V4-V5 926R	CCGYCAATTYMTTTRAGTTT

The V1/V3 workflow uses the primers 27F/534R primers as recommended by (MiDAS version 4)[https://www.biorxiv.org/content/10.1101/2021.07.06.451231v1]. These results weren't great because the length of the final amplified product wasn't long enough to overlap the forward and reverse reads reliably.

V1-V3 27F	AGAGTTTGATCCTGGCTCAG
V1-V3 534R	ATTACCGCGGCTGCTGG
