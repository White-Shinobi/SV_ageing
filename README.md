# SV_ageing
This repository contains the code for the manuscript "Microbial Adaptation in Healthy Ageing: Insights from Age-associated Structural Variation in the Human Gut Microbiome".
# Paper

- Title: Microbial Adaptation in Healthy Ageing: Insights from Age-associated Structural Variation in the Human Gut Microbiome
- Journal: ^_^ [in revision]
- Year: 2025

# Summary

This repository contains the code for the manuscript "Microbial Adaptation in Healthy Ageing: Insights from Age-associated Structural Variation in the Human Gut Microbiome".

# Contents
The code is organised in three folders:

- SV_analysis contains scripts for the gut microbial SV and host age analysis and downstream microbiome-related analyses
- Gene_and_gene_family_analysis contains scripts for the gut microbial genes/gene families and host age/health-related phenotypes analysis
- Functions contains scripts for all the functions used in paper

## Description of each folder
Here we describe the general workflow, for details please see the comments and descriptions provided in each file.

### SV_analysis folder

- **Figure_1_step1.R**: run associations between SV and age in DMP cohort

- **Figure_1_step2.R**: run associations between SV and age in other cohorts, check replication and do plots; to get the 105 SVs

- **Generation_annotation_file.R**: generate SV gene+function annotation (two softwares)

- **Shortbred_marker_generation.R**: generate protein markers for shortbred

### Gene_and_gene_family_analysis folder

- **1.Shortbred_gene_phenotype_DMP.R**: run shortbred associations in DMP cohort (run on calculation cluster)

- **2.Figure_2.R**: draw Figure 2 plots

- **3.Strain_MAG_short_version.R**: run MAG analysis for Oscillibacter sp. ER4.

- **4.Figure_4.R**: sraw Figure 4 plots+calculate numbers

- **5.Shortbred_gene_phenotype_LLD.R**: run shortbred associations in LLD

- **6.Humann_gene_phenotype_DMP.R**: associations for Humann and age and phenotypes in DMP

- **7.Humann_1.R**: summarize humann gene ~ phenotypes results+volcano plot+Humann gene annotation file+Enrichment analysis for gene families

- **8.Figure4_humann.R**: forest plot for humann data

- **9.Shortbred_Humann_gene_metabolites_LLD.R**: analysis of Shortbred+Humann genes ~ metabolites in LLD

- **10.Gene_metabolites_1.R**: summarize and plot for Figure5

- **11.Gene_metabolites_2.R**: summarize and plot for humann data

### Functions folder

- **Part1_functions.R**: association analysis

- **Part2_functions.R**: replication test

- **Part3_functions.R**: strain analysis

NOTE: Most scripts point to paths on our cluster, so if you want to replicate any of the scripts and need input file format descriptions, please contact us:

- Dr. Yue Zhang, yuezhang.umcg@gmail.com

