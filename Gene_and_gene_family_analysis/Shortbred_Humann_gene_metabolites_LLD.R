source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
######################################################
######### Shortbred gene ~ metabolites in LLD #########
######################################################
#### 0. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
# abundance data
load("./abundance_rawdata/s_abun.RData")
all_sv_info_anno <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/all_sv_info_anno.rds")
# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
# para <- readRDS("~/Downloads/2023_09_02/running_list_higher_30.rds")
# para <- readRDS("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
# physical score
#physical_score <- read.delim("~/Downloads/2023_09_02/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
#Frailty_index <- read.delim("~/Downloads/2023_09_02/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
#trans_ID <- read.delim("~/Downloads/2023_09_02/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index=merge(Frailty_index,trans_ID[,1:3],by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$DAG3_sampleID),c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU","DAG3_sampleID")]%>% `rownames<-`(.[,'DAG3_sampleID']) %>% dplyr::select(-'DAG3_sampleID')
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
full_phen=merge(full_phen,Frailty_index[,c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# pheno_cov
# read counts data
load("./readcounts_rawdata/full_read.RData")
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
sgv_full=merge(dsgv_full, vsgv_full, by.x = "row.names", by.y = "row.names", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
Shortbred <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/Shortbred.rds")
LLD_ID <- read.delim("/scratch/p303998/SV_MWAS/Shortbred/key_LLD_GoNL_1659samples_v2.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
LLD_ID$LLD_bam_id=gsub("fece_","",LLD_ID$LLD_bam_id)
Shortbred_LLD=merge(Shortbred,LLD_ID[,c("LLD_GoNL_all_id","LLD_bam_id")],by.x = "row.names",by.y = "LLD_bam_id",all.x = T)
Shortbred_LLD$Row.names=NULL
Shortbred_LLD=Shortbred_LLD[!is.na(Shortbred_LLD$LLD_GoNL_all_id),]
row.names(Shortbred_LLD)=Shortbred_LLD$LLD_GoNL_all_id
Shortbred_LLD$LLD_GoNL_all_id=NULL
saveRDS(Shortbred_LLD,paste0("/scratch/p303998/SV_MWAS/NM_revision/","Shortbred_LLD",".rds"))
Shortbred_LLD <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/Shortbred_LLD.rds")
#### 1.1 Shortbred cut-phenotypes running ####
cohort = c("lld1")
metabolites_LLD <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/data_1442samples_LLD_baseline_1183plasma_metabolites.txt",sep="\t",row.names = 1,header = T,check.names = F,fill = F)
pheno_cohort = metabolites_LLD

#### 1.2 Template for the script (Habrak) ####
script_template <- '
# Habrok
setwd("/scratch/p303998/SV_MWAS/Rdata_1216/")
# functions
source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#### 1. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")
# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("lld1")]=c("lld1")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# physical score
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID_LLD <- read.delim("/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/key_lld_to_ll2_ids.txt",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
colnames(trans_ID_LLD)[3]=c("PROJECT_PSEUDO_ID")
Frailty_index=merge(Frailty_index,trans_ID_LLD,by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$LLDEEP_ID),c("FI41_B","LLDEEP_ID")]%>% `rownames<-`(.[,"LLDEEP_ID"]) %>% dplyr::select(-"LLDEEP_ID")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
full_phen=merge(full_phen,Frailty_index,by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
# read counts data
load("./readcounts_rawdata/full_read.RData")
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
# pheno_cov
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration", "Age", "Sex","BristolType")], full_read[,c("V1", "log10_counts"), drop=FALSE], by.x = "row.names", by.y = "V1", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# vSV and dSV together
sgv_full=merge(dsgv_full, vsgv_full, by.x = "row.names", by.y = "row.names", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# shortbred
Shortbred_LLD <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/Shortbred_LLD.rds")
metabolites_LLD <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/data_1442samples_LLD_baseline_1183plasma_metabolites.txt",sep="\t",row.names = 1,header = T,check.names = F,fill = F)
pheno_cohort = metabolites_LLD

#### lld1 ####
cohort = c("lld1")
pheno_cov_abundance = pheno_cov
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Shortbred_LLD), row.names(pheno_cov_abundance)))
covar = c("Sex","Age") #,"log10_counts"
result = shortbred_correct_noSV(Shortbred_LLD[sample_name,{NUM_RANGE}], pheno_cohort[sample_name,,drop=F],sgv_full[sample_name,],
                                pheno_cov_abundance[sample_name,], covar,colnames(Shortbred_LLD)[{NUM_RANGE}])
for (i in c("Beta_IDage","SE_IDage","p_IDage", "Beta_phenotype","SE_phenotype","p_phenotype","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`result`[,i]=as.numeric(`result`[,i])}
saveRDS(result,paste0("/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/metabolites_analysis/","_range_{NUM_RANGE}.rds"))
'
#### 1.3 generate scripts ####
generate_scripts <- function(start, end, step) {
  for (i in seq(start, end, step)) {
    num_range <- paste(i, min(i + step - 1, end), sep = ":")
    script_content <- gsub("\\{NUM_RANGE\\}", num_range, script_template)
    script_name <- paste0("analysis_script_", num_range, ".R")
    writeLines(script_content, con = script_name)
  }
}
setwd("/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/metabolites_analysis/")
generate_scripts(1, ncol(Shortbred_LLD), 10)
#### 1.4 linux submit ####
# cd /scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/all_phenos/
# vim submit_jobs.sh
#!/bin/bash
# SCRIPT_DIR="/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/all_phenos/" #####!!!!! Don't forget to change the folder name.
# module load RPlus
# for script in $SCRIPT_DIR/*.R
# do
# ls $script
# sbatch --job-name=myRjob --output=output_%j.txt --time=01:00:00 --cpus-per-task=4 --mem=10G --wrap="module load RPlus; Rscript $script"
# done
# chmod +x submit_jobs.sh

#### 1.5 rbind all phenotypes results ####
setpath <- "/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/metabolites_analysis/"
list = list.files(setpath,pattern = '.rds',recursive = F,full.names = TRUE)%>%unique()
dataframes_list <- lapply(list, readRDS)
combined_dataframe <- do.call(rbind, dataframes_list)
length(unique(combined_dataframe$X))
length(unique(combined_dataframe$Y))
length(unique(combined_dataframe$X))*length(unique(combined_dataframe$Y))==nrow(combined_dataframe)
metabolites_info <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/key_lld_1183meta_annotation.tsv",sep="\t",row.names = 1,header = T,check.names = F,fill = F)
combined_dataframe=merge(combined_dataframe,metabolites_info,by.x = "X", by.y = "row.names", all.x = T)
combined_dataframe$Var.20=NULL;combined_dataframe$Var.21=NULL;combined_dataframe$Var.22=NULL;combined_dataframe$Var.23=NULL
combined_dataframe$cs_description=NULL
saveRDS(combined_dataframe,paste0("/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/","Shortbred_gene_metabolites_LLD",".rds"))

######################################################
######### Humann gene ~ metabolites in LLD #########
######################################################
source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#### 0. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
# abundance data
load("./abundance_rawdata/s_abun.RData")
all_sv_info_anno <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/all_sv_info_anno.rds")
# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
# para <- readRDS("~/Downloads/2023_09_02/running_list_higher_30.rds")
# para <- readRDS("/groups/umcg-fu/tmp01/users/umcg-yzh/R/2023_09_13/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
# physical score
#physical_score <- read.delim("~/Downloads/2023_09_02/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
#Frailty_index <- read.delim("~/Downloads/2023_09_02/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
#trans_ID <- read.delim("~/Downloads/2023_09_02/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Frailty_index=merge(Frailty_index,trans_ID[,1:3],by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$DAG3_sampleID),c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU","DAG3_sampleID")]%>% `rownames<-`(.[,'DAG3_sampleID']) %>% dplyr::select(-'DAG3_sampleID')
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
full_phen=merge(full_phen,Frailty_index[,c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# pheno_cov
# read counts data
load("./readcounts_rawdata/full_read.RData")
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
Humann4.0_LLD <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/Humann4.0_LLD.rds")
#### 1.1 Humann cut-phenotypes running ####
cohort = c("lld1")
metabolites_LLD <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/data_1442samples_LLD_baseline_1183plasma_metabolites.txt",sep="\t",row.names = 1,header = T,check.names = F,fill = F)
pheno_cohort = metabolites_LLD
#### 1.2 Template for the script (Habrak) ####
script_template <- '
# Habrok
setwd("/scratch/p303998/SV_MWAS/Rdata_1216/")
# functions
source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part1_functions.R"))
#### 1. import data ####
setwd("/scratch/p303998/SV_MWAS/Rdata_0828/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")
# sample number data
para <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("lld1")]=c("lld1")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
# physical score
physical_score <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
# Frailty index
Frailty_index <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID_LLD <- read.delim("/scratch/p303998/SV_MWAS/NM_revision/Shortbred/LLD/key_lld_to_ll2_ids.txt",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
colnames(trans_ID_LLD)[3]=c("PROJECT_PSEUDO_ID")
Frailty_index=merge(Frailty_index,trans_ID_LLD,by.x = "PROJECT_PSEUDO_ID",by.y = "PROJECT_PSEUDO_ID",all.x=T)
Frailty_index=Frailty_index[!is.na(Frailty_index$LLDEEP_ID),c("FI41_B","LLDEEP_ID")]%>% `rownames<-`(.[,"LLDEEP_ID"]) %>% dplyr::select(-"LLDEEP_ID")
# Phenotype data
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
full_phen=merge(full_phen,Frailty_index,by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,"Row.names"]) %>% dplyr::select(-"Row.names")
# read counts data
load("./readcounts_rawdata/full_read.RData")
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
# pheno_cov
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration", "Age", "Sex","BristolType")], full_read[,c("V1", "log10_counts"), drop=FALSE], by.x = "row.names", by.y = "V1", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# vSV and dSV together
sgv_full=merge(dsgv_full, vsgv_full, by.x = "row.names", by.y = "row.names", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
# Humann
Humann4.0_LLD <- readRDS("/scratch/p303998/SV_MWAS/NM_revision/Humann4.0_LLD.rds")
metabolites_LLD <- read.delim("/scratch/p303998/SV_MWAS/Rdata_1216/data_1442samples_LLD_baseline_1183plasma_metabolites.txt",sep="\t",row.names = 1,header = T,check.names = F,fill = F)
pheno_cohort = metabolites_LLD

#### lld1 ####
cohort = c("lld1")
pheno_cov_abundance = pheno_cov
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_LLD), row.names(pheno_cov_abundance)))
covar = c("Age","Sex","log10_counts")
result = genes_correct_noSV(Humann4.0_LLD[sample_name,{NUM_RANGE}], pheno_cohort[sample_name,,drop=F],sgv_full[sample_name,],
                                pheno_cov_abundance[sample_name,], covar,colnames(Humann4.0_LLD)[{NUM_RANGE}])
for (i in c("Beta_IDage","SE_IDage","p_IDage", "Beta_phenotype","SE_phenotype","p_phenotype","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`result`[,i]=as.numeric(`result`[,i])}
saveRDS(result,paste0("/scratch/p303998/SV_MWAS/NM_revision/Humann/LLD/all_phenos/metabolites/","_range_{NUM_RANGE}.rds"))
'
#### 1.3 generate scripts ####
generate_scripts <- function(start, end, step) {
  for (i in seq(start, end, step)) {
    num_range <- paste(i, min(i + step - 1, end), sep = ":")
    script_content <- gsub("\\{NUM_RANGE\\}", num_range, script_template)
    script_name <- paste0("analysis_script_", num_range, ".R")
    writeLines(script_content, con = script_name)
  }
}
setwd("/scratch/p303998/SV_MWAS/NM_revision/Humann/LLD/all_phenos/metabolites/")
generate_scripts(1, ncol(Humann4.0_LLD), 20)

#### 1.4 linux submit ####
# cd /scratch/p303998/SV_MWAS/NM_revision/scripts_all_phenos/
# vim submit_jobs.sh
#!/bin/bash
# SCRIPT_DIR="/scratch/p303998/SV_MWAS/NM_revision/scripts_all_phenos/" #####!!!!! Don't forget to change the folder name.
# module load R/4.2.1-foss-2022a
# for script in $SCRIPT_DIR/*.R
# do
# ls $script
# sbatch --job-name=myRjob --output=output_%j.txt --time=01:00:00 --cpus-per-task=4 --mem=50G --wrap="module load R/4.2.1-foss-2022a; Rscript $script"
# done
# chmod +x submit_jobs.sh

#### 1.5 rbind all phenotypes results ####
setpath <- "/scratch/p303998/SV_MWAS/NM_revision/Humann/LLD/all_phenos/metabolites/"
list = list.files(setpath,pattern = '.rds',recursive = F,full.names = TRUE)%>%unique()
dataframes_list <- lapply(list, readRDS)
combined_dataframe <- do.call(rbind, dataframes_list)
length(unique(combined_dataframe$X))
length(unique(combined_dataframe$Y))
length(unique(combined_dataframe$X))*length(unique(combined_dataframe$Y))==nrow(combined_dataframe)
saveRDS(combined_dataframe,paste0("/scratch/p303998/SV_MWAS/NM_revision/","Humann_gene_metabolites_LLD",".rds"))
