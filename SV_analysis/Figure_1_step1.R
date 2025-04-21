# Association analysis in 5 cohorts (model:SV~age+sex+DNA.concentrate+readcounts+species abundance)
# functions
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part1_functions.R"))
#### 1.0 Import and prepare data ####
setwd("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/")
# SV data
load("./SV_rawdata/dsgv_full.RData")
load("./SV_rawdata/vsgv_full.RData")
# abundance data
load("./abundance_rawdata/s_abun.RData")
para <- readRDS("./running_list_higher_30.rds")
sample_number <- read.delim("./Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
physical_score <- read.delim("./202312_physical_activity_DAG3_LLD_To-Yue.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
physical_score$ID=physical_score$DAG3_sampleID
physical_score$ID[is.na(physical_score$ID)]=physical_score$LLDEEP_ID[is.na(physical_score$ID)]
row.names(physical_score)=physical_score$ID
physical_score=physical_score[row.names(physical_score)%in%c(sample_number$New_SV_ID),]
Frailty_index <- read.delim("./Frailty_index.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
trans_ID <- read.delim("./key_DAG3_pseudo_id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
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
# read counts data
load("./readcounts_rawdata/full_read.RData")
# SV_info
dsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("./SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
all_sgv_info_anno_ld = rbind(dsgv_info_anno_ld,vsgv_info_anno_ld)
load("./SV_info/info.RData")
all_sv_info_anno <- readRDS("./SV_annotation/all_sv_info_anno.rds")
#load("/scratch/p303998/SV_MWAS/Rdata_0828/SV_annotation/svAnnoDb.RData")
load("./SV_annotation/svAnnoDb.RData")
# pheno_cov
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
sgv_full=merge(dsgv_full,vsgv_full,by.x = "row.names",by.y = "row.names",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')


#### 2. relative abundance ####
#### DMP ####
sp_re_dmp <- s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("DMP")],]
sp_re_dmp_clr=do_clr_externalWeighting(t(sp_re_dmp),t(sp_re_dmp)) %>% data.frame(.)  %>% t(.); sp_re_dmp_clr=as.data.frame(sp_re_dmp_clr)
#### LLD-baseline ####
sp_re_lld_baseline <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("LLD1")],])
sp_re_lld1_clr=do_clr_externalWeighting(t(sp_re_lld_baseline),t(sp_re_lld_baseline)) %>% data.frame(.)  %>% t(.); sp_re_lld1_clr=as.data.frame(sp_re_lld1_clr)
#### 500FG-FSK ####
sp_re_500FG_FSK <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("500FG_FSK")],])
sp_re_500fg_fsk_clr=do_clr_externalWeighting(t(sp_re_500FG_FSK),t(sp_re_500FG_FSK)) %>% data.frame(.)  %>% t(.); sp_re_500fg_fsk_clr=as.data.frame(sp_re_500fg_fsk_clr)
#### 300OB ####
sp_re_300OB <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("300OB")],])
sp_re_300ob_clr=do_clr_externalWeighting(t(sp_re_300OB),t(sp_re_300OB)) %>% data.frame(.)  %>% t(.); sp_re_300ob_clr=as.data.frame(sp_re_300ob_clr)
#### 300TZFG ####
sp_re_300TZFG <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("300TZFG")],])
sp_re_300tzfg_clr=do_clr_externalWeighting(t(sp_re_300TZFG),t(sp_re_300TZFG)) %>% data.frame(.)  %>% t(.); sp_re_300tzfg_clr=as.data.frame(sp_re_300tzfg_clr)

#### 3.1 Age association analysis in 5 cohorts : covar = c("DNA.Concentration", "Sex", "log10_counts","species abundance") ####
# model:SV~age+sex+DNA.concentrate+readcounts+species abundance
no_cores <- 6  # Leave one core free for system stability
registerDoParallel(cores = no_cores)
foreach(cohort = c("dmp","300ob","300tzfg","500fg_fsk","lld1"), .combine = 'rbind', .packages = c('dplyr'),.verbose = TRUE) %dopar% {
  running_info_cohort = para[para$V1 == cohort,]
  pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],]
  pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
    `rownames<-`(.[, 'Row.names']) %>%
    dplyr::select(-'Row.names')
  sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
  
  # Determine covariates based on cohort
  covar = switch(cohort,
                 "dmp" = c("DNA.Concentration", "Sex", "log10_counts"),
                 "300ob" = c("Sex", "log10_counts"),
                 "300tzfg" = c("Sex", "log10_counts"),
                 "500fg_fsk" = c("Sex", "log10_counts"),
                 "lld1" = c("Sex", "log10_counts"))
  
  result = SV_lm_glm(pheno_cohort[sample_name, "Age",drop=F], sgv_full[sample_name,],
                     pheno_cov_abundance[sample_name,], covar, running_info_cohort,cohort)
  saveRDS(result, paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Results/", "pheno_is_Age_",cohort, ".rds"))
}

#### 3.2 Age (covar = c("Sex") ) ####
no_cores <- 6  # Leave one core free for system stability
registerDoParallel(cores = no_cores)
results <- foreach(cohort = c("dmp","300ob","300tzfg","500fg_fsk","lld1"), 
                   .combine = 'rbind', 
                   .packages = c('dplyr')) %do% {
                     
                     print(cohort)  # Debugging: See if cohort is working
                     
                     running_info_cohort <- para[para$V1 == cohort,]
                     pheno_cohort <- full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],]
                     pheno_cov_abundance <- pheno_cov
                     sample_name <- Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
                     
                     covar <- c("Sex")
                     
                     result <- SV_correct_nospecies(sgv_full[sample_name,], 
                                                    pheno_cohort[sample_name, "Age", drop=F],
                                                    sgv_full[sample_name,],
                                                    pheno_cov_abundance[sample_name,], 
                                                    covar, colnames(sgv_full))
                     
                     result  # Return result
                   }
for (i in c("Beta_IDage","SE_IDage","p_IDage", "Beta_phenotype","SE_phenotype","p_phenotype","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`results`[,i]=as.numeric(`results`[,i])}
saveRDS(results, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/SV_age_all_cohorts_cov_only_sex.rds")

#### 4.checking association ####
#### 4.1 300ob ####
cohort=c("300ob")
running_info_cohort=para[para$V1==cohort,]
pheno_cohort=full_phen[row.names(full_phen)%in%sample_number$New_SV_ID[sample_number$Cohort_2==cohort],]
pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
  `rownames<-`(.[, 'Row.names']) %>%
  dplyr::select(-'Row.names')
sample_name=Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Butyrivibrio crossotus DSM 2876:2006_2007")], X_phenotype = full_phen[sample_name,c("Age")],pheno_cov_abundance[sample_name,c("Sex","log10_counts","Butyrivibrio crossotus DSM 2876")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Ruminococcus bromii L2-63:66_70;80_97")], X_phenotype = full_phen[sample_name,c("BMI")],pheno_cov_abundance[sample_name,c("Sex","Age","log10_counts","Ruminococcus bromii L2-63")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,grep("\\[Eubacterium\\] eligens ATCC 27750",colnames(sgv_full))[92]], X_phenotype = full_phen[sample_name,c("TCA")],pheno_cov_abundance[sample_name,c("Sex","Age","log10_counts","[Eubacterium] eligens ATCC 27750")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
print(lm_input$X_phenotype)

if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}

if (length(unique(lm_input$Y_SV)) > 2){for (b in c("Y_SV")){lm_input[,b]=qtrans(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}

if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
for (a in c("log10_counts")){lm_input[,a]=qtrans(lm_input[,a])}
lm_res <- summary(lm(Y_SV~.,data = lm_input))
lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial"))
rownames(lm_res$coefficients)[2]
lm_input$Y_SV=as.factor(lm_input$Y_SV)
ggplot(lm_input, aes(y = X_phenotype, x = Y_SV)) + 
  geom_violin()+geom_point()
# Butyrivibrio crossotus DSM 2876:2006_2007:
# > model=glm(Y_SV~.,data = lm_input,family = "binomial")
# Warning messages:
#1: glm.fit: algorithm did not converge 
#2: glm.fit: fitted probabilities numerically 0 or 1 occurred 
# Then use library(detectseparation) and glmet to double check
library(detectseparation)
model=glm(Y_SV~.,data = lm_input,family = "binomial")
response <- lm_input$Y_SV
predictors <- model.matrix(~ ., data = lm_input)[,-1]
separation <- detect_separation(y = response, x = predictors)
library(glmnet)
x <- model.matrix(Y_SV ~ ., data = lm_input)[,-1]
y <- lm_input$Y_SV
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)
print(lasso_model);coef(lasso_model, s = "lambda.min")

#### 4.1 lld1 ####
cohort=c("lld1")
running_info_cohort=para[para$V1==cohort,]
pheno_cohort=full_phen[row.names(full_phen)%in%sample_number$New_SV_ID[sample_number$Cohort_2==cohort],]
pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
  `rownames<-`(.[, 'Row.names']) %>%
  dplyr::select(-'Row.names')
sample_name=Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Bacteroides massiliensis B84634 = Timone 84634 = DSM 17679 = JCM 13223:2543_2555;4404_4414")], X_phenotype = full_phen[sample_name,c("Age")],pheno_cov_abundance[sample_name,c("Sex","log10_counts","Bacteroides massiliensis B84634 = Timone 84634 = DSM 17679 = JCM 13223")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Ruminococcus bromii L2-63:66_70;80_97")], X_phenotype = full_phen[sample_name,c("BMI")],pheno_cov_abundance[sample_name,c("Sex","Age","log10_counts","Ruminococcus bromii L2-63")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
print(lm_input$X_phenotype)

if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}

if (length(unique(lm_input$Y_SV)) > 2){for (b in c("Y_SV")){lm_input[,b]=qtrans(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}

if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_SV~.,data = lm_input))
lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial"))
rownames(lm_res$coefficients)[2]
#### 4.3 dmp ####
cohort=c("dmp")
pheno_cohort=full_phen[row.names(full_phen)%in%sample_number$New_SV_ID[sample_number$Cohort_2==cohort],]
pheno_cov_abundance = merge(pheno_cov, get(paste0("sp_re_", cohort, "_clr")), by.x = "row.names", by.y = "row.names", all = FALSE) %>%
  `rownames<-`(.[, 'Row.names']) %>%
  dplyr::select(-'Row.names')
sample_name=Reduce(intersect, list(row.names(pheno_cohort), row.names(sgv_full), row.names(pheno_cov_abundance)))
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Bacteroides massiliensis B84634 = Timone 84634 = DSM 17679 = JCM 13223:2543_2555;4404_4414")], X_phenotype = full_phen[sample_name,c("Age")],pheno_cov_abundance[sample_name,c("DNA.Concentration","Sex","log10_counts","Bacteroides massiliensis B84634 = Timone 84634 = DSM 17679 = JCM 13223")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Ruminococcus bromii L2-63:66_70;80_97")], X_phenotype = full_phen[sample_name,c("BMI")],pheno_cov_abundance[sample_name,c("DNA.Concentration","Age","log10_counts","Sex","Ruminococcus bromii L2-63")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Oscillibacter sp. ER4:2001_2003")], X_phenotype = full_phen[sample_name,c("Age")],pheno_cov_abundance[sample_name,c("DNA.Concentration","log10_counts","Sex","Oscillibacter sp. ER4")],Raw_age = full_phen[sample_name,c("Age")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Bacteroides fragilis NCTC 9343:1392_1397")], X_phenotype = full_phen[sample_name,c("Antibiotics")],pheno_cov_abundance[sample_name,c("Age","DNA.Concentration","log10_counts","Sex","Bacteroides fragilis NCTC 9343")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Subdoligranulum sp. 4_3_54A2FAA:4198_4199")], X_phenotype = full_phen[sample_name,c("PPI")],pheno_cov_abundance[sample_name,c("Age","DNA.Concentration","log10_counts","Sex","Subdoligranulum sp. 4_3_54A2FAA")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Oscillibacter sp. ER4:459_461 and 2 segments")], X_phenotype = full_phen[sample_name,c("OsmoticLaxatives")],pheno_cov_abundance[sample_name,c("Age","DNA.Concentration","log10_counts","Sex","Oscillibacter sp. ER4")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
lm_input<-data.frame(Y_SV = sgv_full[sample_name,c("Parabacteroides distasonis ATCC 8503:852_854;863_867")], X_phenotype = full_phen[sample_name,c("Biguanides","Type2Diabetes")],pheno_cov_abundance[sample_name,c("Age","DNA.Concentration","log10_counts","Sex","Parabacteroides distasonis ATCC 8503")])  %>% sapply(as.numeric) %>% na.omit %>% as.data.frame

print(lm_input$X_phenotype)

if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}

if (length(unique(lm_input$Y_SV)) > 2){for (b in c("Y_SV")){lm_input[,b]=qtrans(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}

if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

length(unique(lm_input$Y_SV))
lm_res <- summary(lm(Y_SV~.,data = lm_input[,1:6]))
lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial"))
lm_input$Y_SV_red=resid(lm(Y_SV~Sex+log10_counts+DNA.Concentration+Oscillibacter.sp..ER4,data=lm_input))

lm_input$Y_SV=as.factor(lm_input$Y_SV)
ggplot(lm_input, aes(x = X_phenotype, fill=Y_SV)) +
  geom_bar(position = "fill") +
  labs(x = "Antibiotics usage", y = "Proportion of \nBacteroides fragilis NCTC 9343:1392_1397")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14),)
lm_input$X_phenotype=as.factor(lm_input$X_phenotype)
ggplot(data = lm_input, mapping = aes(x =`X_phenotype`, y =`Y_SV`))+
  geom_violin(color="brown")+
  geom_jitter(height = 0, width = 0.1,alpha=0.6)+
  labs(x = "Biguanides usage", y = "Parabacteroides distasonis ATCC 8503:\n852_854;863_867")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=18),
        axis.title.x = element_text(face="bold",size=18),)
#### 4.4 Figure 1B ####
ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_SV`))+
  geom_point(color="brown")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_SV`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Oscillibacter sp. ER4:2001_2003")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14),)
ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_SV_red`))+
  geom_point(color="brown")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_SV_red`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of Oscillibacter sp. ER4:2001_2003")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14),)


#### 5. load all files ####
# association results
setpath <- "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Results/"
list = list.files(setpath,pattern = '.rds',recursive = TRUE,full.names = TRUE)%>%unique()
for (i in list) {
  var_name <- gsub(".rds", "", sapply(strsplit(as.character(i), "/"), "[", 9))
  var_name <- gsub("pheno_is_Age_", "", var_name)
  print(var_name)
  assign(var_name, readRDS(i))
}

#### 6. filter the results according to the zero rate ####
# dmp
dmp$p=as.numeric(dmp$p)
dmp$p[which(dmp$p==0)]=NA
for (i in c("X","Y")){dmp[,i]=as.character(dmp[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){dmp[,i]=as.numeric(dmp[,i])}
dmp$species=sapply(strsplit(as.character(dmp$Y), "\\:"), "[", 1)
dmp$species[grep("Phascolarctobacterium sp. CAG:207",dmp$Y)]=c("Phascolarctobacterium sp. CAG:207")
dmp$species[grep("Phascolarctobacterium sp. CAG:266",dmp$Y)]=c("Phascolarctobacterium sp. CAG:266")
# one level ---> 0
dmp$Beta[which(dmp$Beta==c("one level"))]=NA
dmp$SE[which(dmp$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
dmp$p[dmp$N<sum(sample_number$Cohort_2==c("dmp"),na.rm = T)*0.01]=NA
dmp$p[dmp$y_uniq_N==2&dmp$y_non_zero_rate>0.9]=NA
dmp$p[dmp$y_uniq_N==2&dmp$y_non_zero_rate<0.1]=NA
dmp$p[dmp$x_uniq_N==2&dmp$x_non_zero_rate>0.9]=NA
dmp$p[dmp$x_uniq_N==2&dmp$x_non_zero_rate<0.1]=NA
dmp$fdr.p=NULL
dmp$bonferroni.p=NULL

# lld1
lld1$p=as.numeric(lld1$p)
lld1$p[which(lld1$p==0)]=NA
lld1$species=sapply(strsplit(as.character(lld1$Y), "\\:"), "[", 1)
lld1$species[grep("Phascolarctobacterium sp. CAG:207",lld1$Y)]=c("Phascolarctobacterium sp. CAG:207")
lld1$species[grep("Phascolarctobacterium sp. CAG:266",lld1$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){lld1[,i]=as.character(lld1[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){lld1[,i]=as.numeric(lld1[,i])}
# one level ---> 0
lld1$Beta[which(lld1$Beta==c("one level"))]=NA
lld1$SE[which(lld1$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
lld1$p[lld1$N<sum(sample_number$Cohort_2==c("lld1"),na.rm = T)*0.1]=NA
lld1$p[lld1$y_uniq_N==2&lld1$y_non_zero_rate>0.9]=NA
lld1$p[lld1$y_uniq_N==2&lld1$y_non_zero_rate<0.1]=NA
lld1$p[lld1$x_uniq_N==2&lld1$x_non_zero_rate>0.9]=NA
lld1$p[lld1$x_uniq_N==2&lld1$x_non_zero_rate<0.1]=NA
lld1$fdr.p=NULL
lld1$bonferroni.p=NULL

# 300ob
`300ob`$p=as.numeric(`300ob`$p)
`300ob`$p[which(`300ob`$p==0)]=NA
`300ob`$species=sapply(strsplit(as.character(`300ob`$Y), "\\:"), "[", 1)
`300ob`$species[grep("Phascolarctobacterium sp. CAG:207",`300ob`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`300ob`$species[grep("Phascolarctobacterium sp. CAG:266",`300ob`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`300ob`[,i]=as.character(`300ob`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`300ob`[,i]=as.numeric(`300ob`[,i])}
# one level ---> 0
`300ob`$Beta[which(`300ob`$Beta==c("one level"))]=NA
`300ob`$SE[which(`300ob`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`300ob`$p[`300ob`$N<sum(sample_number$Cohort_2==c("300ob"),na.rm = T)*0.1]=NA
`300ob`$p[`300ob`$y_uniq_N==2&`300ob`$y_non_zero_rate>0.9]=NA
`300ob`$p[`300ob`$y_uniq_N==2&`300ob`$y_non_zero_rate<0.1]=NA
`300ob`$p[`300ob`$x_uniq_N==2&`300ob`$x_non_zero_rate>0.9]=NA
`300ob`$p[`300ob`$x_uniq_N==2&`300ob`$x_non_zero_rate<0.1]=NA
`300ob`$fdr.p=NULL
`300ob`$bonferroni.p=NULL

# 300tzfg
`300tzfg`$p=as.numeric(`300tzfg`$p)
`300tzfg`$p[which(`300tzfg`$p==0)]=NA
`300tzfg`$species=sapply(strsplit(as.character(`300tzfg`$Y), "\\:"), "[", 1)
`300tzfg`$species[grep("Phascolarctobacterium sp. CAG:207",`300tzfg`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`300tzfg`$species[grep("Phascolarctobacterium sp. CAG:266",`300tzfg`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`300tzfg`[,i]=as.character(`300tzfg`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`300tzfg`[,i]=as.numeric(`300tzfg`[,i])}
# one level ---> 0
`300tzfg`$Beta[which(`300tzfg`$Beta==c("one level"))]=NA
`300tzfg`$SE[which(`300tzfg`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`300tzfg`$p[`300tzfg`$N<sum(sample_number$Cohort_2==c("300tzfg"),na.rm = T)*0.1]=NA
`300tzfg`$p[`300tzfg`$y_uniq_N==2&`300tzfg`$y_non_zero_rate>0.9]=NA
`300tzfg`$p[`300tzfg`$y_uniq_N==2&`300tzfg`$y_non_zero_rate<0.1]=NA
`300tzfg`$p[`300tzfg`$x_uniq_N==2&`300tzfg`$x_non_zero_rate>0.9]=NA
`300tzfg`$p[`300tzfg`$x_uniq_N==2&`300tzfg`$x_non_zero_rate<0.1]=NA
`300tzfg`$fdr.p=NULL
`300tzfg`$bonferroni.p=NULL

# 500fg_fsk
`500fg_fsk`$p=as.numeric(`500fg_fsk`$p)
`500fg_fsk`$p[which(`500fg_fsk`$p==0)]=NA
`500fg_fsk`$species=sapply(strsplit(as.character(`500fg_fsk`$Y), "\\:"), "[", 1)
`500fg_fsk`$species[grep("Phascolarctobacterium sp. CAG:207",`500fg_fsk`$Y)]=c("Phascolarctobacterium sp. CAG:207")
`500fg_fsk`$species[grep("Phascolarctobacterium sp. CAG:266",`500fg_fsk`$Y)]=c("Phascolarctobacterium sp. CAG:266")
for (i in c("X","Y")){`500fg_fsk`[,i]=as.character(`500fg_fsk`[,i])}
for (i in c("Beta","SE","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){`500fg_fsk`[,i]=as.numeric(`500fg_fsk`[,i])}
# one level ---> 0
`500fg_fsk`$Beta[which(`500fg_fsk`$Beta==c("one level"))]=NA
`500fg_fsk`$SE[which(`500fg_fsk`$SE==c("one level"))]=NA
# dSV level ---> 0.1<zero number<0.9
`500fg_fsk`$p[`500fg_fsk`$N<sum(sample_number$Cohort_2==c("500fg_fsk"),na.rm = T)*0.1]=NA
`500fg_fsk`$p[`500fg_fsk`$y_uniq_N==2&`500fg_fsk`$y_non_zero_rate>0.9]=NA
`500fg_fsk`$p[`500fg_fsk`$y_uniq_N==2&`500fg_fsk`$y_non_zero_rate<0.1]=NA
`500fg_fsk`$p[`500fg_fsk`$x_uniq_N==2&`500fg_fsk`$x_non_zero_rate>0.9]=NA
`500fg_fsk`$p[`500fg_fsk`$x_uniq_N==2&`500fg_fsk`$x_non_zero_rate<0.1]=NA
`500fg_fsk`$fdr.p=NULL
`500fg_fsk`$bonferroni.p=NULL

#### 7. save the final files ####
for( i in c("dmp","300ob","300tzfg","500fg_fsk","lld1")){
  print(i)
  saveRDS(get(i),paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/","filter_zero_",i,".rds"))
}


#### 8 Paragraph1: Gut microbial SV profiling overview across five cohorts ####
#### 8.1 108 species ####
# sample size per cohort
samples=sample_number$New_SV_ID[sample_number$Cohort%in%c("300OB","300TZFG","500FG_FSK","DMP","LLD1")]
vsgv_full[which(row.names(vsgv_full)%in%samples),]%>%removeRowsAllNa(.)%>%dim(.)
vsgv_full[which(row.names(vsgv_full)%in%samples),]%>%removeColsAllNa(.)%>%dim(.)
dsgv_full[which(row.names(dsgv_full)%in%samples),]%>%removeRowsAllNa(.)%>%dim(.)
dsgv_full[which(row.names(dsgv_full)%in%samples),]%>%removeColsAllNa(.)%>%dim(.)
test=pheno_cov[which(row.names(pheno_cov)%in%samples),]
test=test[!is.na(test$Age),]
# final number written in the manuscript
table(sample_number$Cohort)
7954+1134+518+290+319
# reference genome number
nrow(info)
# SV number
sum(info$SVs_number)
sum(info$Deletion_SVs_number)
sum(info$Variable_SVs_number)
# relative abundance mean level
sp_re_5cohorts <- na.omit(s_abun[sample_number$New_SV_ID[sample_number$Cohort==c("DMP")|sample_number$Cohort==c("LLD1")|
                                                           sample_number$Cohort==c("500FG_FSK")|sample_number$Cohort==c("300OB")|
                                                           sample_number$Cohort==c("300TZFG")],])
max(rowSums(sp_re_5cohorts[,-ncol(sp_re_5cohorts)]))
min(rowSums(sp_re_5cohorts[,-ncol(sp_re_5cohorts)]))
mean(rowSums(sp_re_5cohorts[,-ncol(sp_re_5cohorts)]))
# Figure S1
sp_re_5cohorts_sums <- data.frame(sums=rowSums(sp_re_5cohorts[, -ncol(sp_re_5cohorts)]))
ggplot(data = sp_re_5cohorts_sums, aes(x = sums)) +
  geom_density(fill = "lightblue", color = "black", alpha = 0.7) +
  geom_rug(sides = "b", color = "dark blue",alpha = 0.5) +  # Add a rug plot
  geom_vline(xintercept = mean(rowSums(sp_re_5cohorts[,-ncol(sp_re_5cohorts)])), linetype = "dashed", color = "dark red", alpha = 0.5) +  # Add "grass" lines+
  labs(title = "Histogram of Row Sums",
       x = "Abundance",
       y = "Density") +
  theme_minimal()+
  theme(axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black"))
# TableS1
write.table(info_profile, "/Users/helloduck/Desktop/info.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
# Figure 1B
# Graphlan galaxy tutorial: https://github.com/biobakery/biobakery/wiki/graphlan#1-graphlan--galaxy-module
# Graphlan input1: info_profile (taxanomical info)
info$Family <- factor(info$Family);info$Genus <- factor(info$Genus)
# numbers per species
info_profile=info
info_profile=merge(info_profile,as.data.frame(colMeans(sp_re_5cohorts,na.rm = T)),by.x = "organism",by.y = "row.names",all=F)
info_profile$dSV_num_perM=info_profile$Deletion_SVs_number/(info_profile$Length/10^6)
info_profile$vSV_num_perM=info_profile$Variable_SVs_number/(info_profile$Length/10^6)
info_profile$SVs_number_perM=info_profile$SVs_number/(info_profile$Length/10^6)
for (i in (colnames(info_profile)[2:9])){info_profile[,i]=factor(info_profile[,i])}
test=as.data.frame(table(info_profile$Family))
info_profile$Family=factor(info_profile$Family,levels=unique(test$Var1[order(test$Freq,decreasing = T)]))
info_profile <- info_profile %>% arrange(Family,SVs_number_perM)
info_profile$organism=factor(info_profile$organism,levels=info_profile$organism)
info_profile$species_num=seq(1:108)
info_profile$IBD=NULL
info_profile$Total_samples_number=info_profile$DMP+info_profile$LLD1+info_profile$X500FG_FSK+info_profile$X300OB+info_profile$X300TZFG
info_profile$all=paste(info_profile$Superkingdom,info_profile$Phylum,info_profile$Class,
                       info_profile$Order,info_profile$Family,info_profile$Genus,
                       info_profile$species_num,sep=".")
info_profile$all=gsub("\\[","",info_profile$all)
info_profile$all=gsub("\\]","",info_profile$all)
info_profile$all=gsub(" ","_",info_profile$all)
# 63 species association analysis in discovery cohort
info_profile$after_filtering=NA
non_NA_final <- readRDS("~/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
info_profile$after_filtering[info_profile$organism%in%unique(non_NA_final$species[non_NA_final$included_downstream==c("Yes")])]=c("Passed_filtering")

#### 8.2 Table S1. The information for 108 reference genomes. ####
write.table(info_profile$all, "/scratch/p303998/SV_MWAS/Rdata_0828/info_profile.txt",sep = "\t", col.names = F, row.names = F, quote = F)
write.table(info_profile, "/Users/helloduck/Desktop/STable1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
# Graphlan input2: annotation_graphlan (ring annotation file)
SV_samsize_1=data.frame(all=info_profile[,c("all")],
                        feature=c("ring_alpha"),
                        ring_num=c("1"),
                        content=rescale(info_profile[,c("Total_samples_number")]))
SV_samsize_2=data.frame(all=info_profile[,c("all")],
                        feature=c("ring_color"),
                        ring_num=c("1"),
                        content=c("#ff595e"))

SV_re_1=data.frame(all=info_profile[,c("all")],
                   feature=c("ring_alpha"),
                   ring_num=c("2"),
                   content=rescale(info_profile[,c("colMeans(s_abun, na.rm = T)")]))
SV_re_2=data.frame(all=info_profile[,c("all")],
                   feature=c("ring_color"),
                   ring_num=c("2"),
                   content=c("#1982c4"))

SV_num=data.frame(all=info_profile[,c("all")],
                  feature=c("ring_height"),
                  ring_num=c("3"),
                  content=rescale(info_profile[,c("SVs_number_perM")])*5)

basic_info_1=data.frame(all=c("ring_internal_separator_thickness"),
                        feature=c("1","2","3"),
                        ring_num=c("1.0"),
                        content="")
basic_info_2=data.frame(all=c("ring_width"),
                        feature=c("3"),
                        ring_num=c("0.5"),
                        content="")
annotation_graphlan=rbind(basic_info_1,basic_info_2,SV_samsize_1,SV_samsize_2,SV_re_1,SV_re_2,SV_num)
write.table(annotation_graphlan, "/scratch/p303998/SV_MWAS/Rdata_0828/annotation_graphlan.tsv",sep = "\t", col.names = F, row.names = F, quote = F)
#numbers per family
table(info$Phylum)
test=as.data.frame(table(info$Genus))
# Figure S1:numbers per genus
info_genus=aggregate(SVs_number ~ Genus, data = info, sum)
test=as.data.frame(table(info$Genus))
info_genus=merge(info_genus,test,by.x = "Genus",by.y = "Var1",all=T)
info_genus$SV_number_pergenus_control=info_genus$SVs_number/info_genus$Freq
ggplot(info_genus, aes(x = reorder(Genus, -SV_number_pergenus_control), y = SV_number_pergenus_control)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Bar Plot of Genus Counts",
       x = "Genus",
       y = "Count")
#Figure S1: sample size
ggplot(info_profile, aes(x = reorder(Short_name, -Total_samples_number), y = Total_samples_number)) +
  geom_bar(stat = "identity",fill=) +
  coord_flip() +
  labs(title = "Bar Plot of Genus Counts",
       x = "Genus",
       y = "Count")

#### 8.3 Table S2. The information for 14,249 SVs. ####
non_NA_final=as.data.frame(colnames(sgv_full));row.names(non_NA_final)=non_NA_final$`colnames(sgv_full)`
for (cohort in c("dmp","300ob", "300tzfg", "500fg_fsk", "lld1")){
  #cohort=c("lld1")
  test=sgv_full[row.names(sgv_full) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort[1]],]
  #test=test[,"Bacteroides uniformis ATCC 8492:130_131 and 15 segments",drop=F]
  non_NA=as.data.frame(colSums(!is.na(test)))
  colnames(non_NA)[1]=paste0("sample_size","_",cohort[1])#Non_NA_sample_size
  test1=colSums(!is.na(test) & test == 0)%>%as.data.frame(.)
  test2=colSums(!is.na(test) & test == 1)%>%as.data.frame(.)
  detection_0_1=merge(test1,test2,by.x = "row.names",by.y = "row.names",all=T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
  colnames(detection_0_1)[1]=paste0("num_of_0","_",cohort[1])
  colnames(detection_0_1)[2]=paste0("num_of_1","_",cohort[1])
  non_NA=merge(non_NA,detection_0_1,by.x = "row.names",by.y = "row.names",all=T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
  non_NA$rate_of_1=non_NA[,3]/non_NA[,1]
  colnames(non_NA)[4]=paste0("rate_of_1","_",cohort[1])
  print(nrow(non_NA_final))
  non_NA_final=merge(non_NA,non_NA_final,by.x = "row.names",by.y = "row.names",all=T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
}
non_NA_final$`colnames(sgv_full)`=NULL
non_NA_final$Detected_in_how_many_other_cohorts=rowSums(!non_NA_final[,grep("sample_size",colnames(non_NA_final))[1:4]]==0)
non_NA_final=merge(all_sv_info_anno[,c("SV_Name","SV_type")],non_NA_final,by.x = "SV_Name",by.y = "row.names",all=T)
non_NA_final$After_filtering=NA
non_NA_final$After_filtering[which(non_NA_final$SV_type==c("dSV")&non_NA_final$sample_size_dmp>80&non_NA_final$rate_of_1_dmp>=0.1&non_NA_final$rate_of_1_dmp<0.9)]=c("passed")
non_NA_final$After_filtering[non_NA_final$SV_type==c("vSV")&non_NA_final$sample_size_dmp>80]=c("passed")
non_NA_final$Pass_cutoff_in_how_many_other_cohorts=0
for (i in c("300ob", "300tzfg", "500fg_fsk", "lld1")){
  cut_value=sum(sample_number$Cohort_2==i,na.rm = T)*0.1
  for (j in 1:nrow(non_NA_final)){
    if (non_NA_final[j,grep(paste0("sample_size_",i),colnames(non_NA_final))]>=cut_value){
      non_NA_final$Pass_cutoff_in_how_many_other_cohorts[j]=non_NA_final$Pass_cutoff_in_how_many_other_cohorts[j]+1
    }
  }
}
non_NA_final$species=sapply(strsplit(as.character(non_NA_final$SV_Name), "\\:"), "[", 1)
non_NA_final$species[grep("Phascolarctobacterium sp. CAG:207",non_NA_final$SV_Name)]=c("Phascolarctobacterium sp. CAG:207")
non_NA_final$species[grep("Phascolarctobacterium sp. CAG:266",non_NA_final$SV_Name)]=c("Phascolarctobacterium sp. CAG:266")
length(unique(non_NA_final$species[which(non_NA_final$After_filtering==c("passed"))]))#63
length(unique(non_NA_final$species[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]))#50
non_NA_final$included_downstream=NA
non_NA_final$included_downstream[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]=c("Yes")
non_NA_final$included_downstream[is.na(non_NA_final$included_downstream)]=c("No")
saveRDS(non_NA_final, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
double_check=non_NA_final[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts==0),grep("size",colnames(non_NA_final))]
colnames(non_NA_final)=gsub("lld1","LLD",colnames(non_NA_final))
colnames(non_NA_final)=gsub("500fg_fsk","500FG",colnames(non_NA_final))
colnames(non_NA_final)=gsub("300tzfg","300TZFG",colnames(non_NA_final))
colnames(non_NA_final)=gsub("300ob","300OB",colnames(non_NA_final))
colnames(non_NA_final)=gsub("dmp","DMP",colnames(non_NA_final))
write.table(non_NA_final, paste0("/Users/helloduck/Desktop/","table1",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)

#### 9.1 63 species ####
filter_zero_dmp <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/Results/All/filter_zero/filter_zero_dmp.rds")
filter_zero_dmp <- readRDS("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/Results/All/filter_zero/filter_zero_dmp.rds")
filter_zero_dmp=filter_zero_dmp[!is.na(filter_zero_dmp$p),]
filter_zero_dmp <- merge(filter_zero_dmp, all_sv_info_anno, by.x = "Y", by.y="SV_Name",all.x = TRUE)
test=as.data.frame(table(filter_zero_dmp[filter_zero_dmp$SV_type==c("dSV"),c("species","Y")]))
test=test[test$Freq>0,]
test=as.data.frame(table(test$species))
colnames(test)[2]=c("Deletion_SVs_number")
Deletion_SVs_number=test
test=as.data.frame(table(filter_zero_dmp[filter_zero_dmp$SV_type==c("vSV"),c("species","Y")]))
test=test[test$Freq>0,]
test=as.data.frame(table(test$species))
colnames(test)[2]=c("Variable_SVs_number")
Variable_SVs_number=test
info_63=merge(Deletion_SVs_number,Variable_SVs_number,by.x = "Var1",by.y = "Var1",all=T)
info_63$SVs_number=info_63$Deletion_SVs_number+info_63$Variable_SVs_number
info_63=merge(info[,setdiff(colnames(info),c("Deletion_SVs_number","Variable_SVs_number","SVs_number"))],info_63,by.x = "organism",by.y = "Var1",all.y  =F)
info_63$species_num=seq(1:63)
info_63$IBD=NULL
info_63$Total_samples_number=info_63$DMP+info_63$LLD1+info_63$X500FG_FSK+info_63$X300OB+info_63$X300TZFG
write.table(info_63, "/scratch/p303998/SV_MWAS/Rdata_1216/info_63.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
# Graphlan input1: info_63_profile (taxanomical info_63)
# numbers per species
info_63_profile=info_63
info_63_profile=merge(info_63_profile,as.data.frame(colMeans(sp_re_dmp,na.rm = T)),by.x = "organism",by.y = "row.names",all=F)
#info_63_profile$dSV_num_perM=info_63_profile$Deletion_SVs_number/(info_63_profile$Length/10^6)
#info_63_profile$vSV_num_perM=info_63_profile$Variable_SVs_number/(info_63_profile$Length/10^6)
#info_63_profile$SVs_number_perM=info_63_profile$SVs_number/(info_63_profile$Length/10^6)
info_63_profile$all=paste(info_63_profile$Superkingdom,info_63_profile$Phylum,info_63_profile$Class,
                       info_63_profile$Order,info_63_profile$Family,info_63_profile$Genus,
                       info_63_profile$species_num,sep=".")
info_63_profile$all=gsub("\\[","",info_63_profile$all)
info_63_profile$all=gsub("\\]","",info_63_profile$all)
info_63_profile$all=gsub(" ","_",info_63_profile$all)
write.table(info_63_profile$all, "/scratch/p303998/SV_MWAS/Rdata_1216/info_63_profile.txt",sep = "\t", col.names = F, row.names = F, quote = F)
# Graphlan input2: annotation_graphlan (ring annotation file)
SV_samsize_1=data.frame(all=info_63_profile[,c("all")],
                        feature=c("ring_alpha"),
                        ring_num=c("1"),
                        content=rescale(info_63_profile[,c("Total_samples_number")]))
SV_samsize_2=data.frame(all=info_63_profile[,c("all")],
                        feature=c("ring_color"),
                        ring_num=c("1"),
                        content=c("#ff595e"))

SV_re_1=data.frame(all=info_63_profile[,c("all")],
                   feature=c("ring_alpha"),
                   ring_num=c("2"),
                   content=rescale(info_63_profile[,c("colMeans(sp_re_dmp, na.rm = T)")]))
SV_re_2=data.frame(all=info_63_profile[,c("all")],
                   feature=c("ring_color"),
                   ring_num=c("2"),
                   content=c("#1982c4"))

SV_num=data.frame(all=info_63_profile[,c("all")],
                  feature=c("ring_height"),
                  ring_num=c("3"),
                  content=rescale(info_63_profile[,c("SVs_number")])*5)

basic_info_63_1=data.frame(all=c("ring_internal_separator_thickness"),
                        feature=c("1","2","3"),
                        ring_num=c("1.0"),
                        content="")
basic_info_63_2=data.frame(all=c("ring_width"),
                        feature=c("3"),
                        ring_num=c("0.5"),
                        content="")
annotation_graphlan=rbind(basic_info_63_1,basic_info_63_2,SV_samsize_1,SV_samsize_2,SV_re_1,SV_re_2,SV_num)
write.table(annotation_graphlan, "/scratch/p303998/SV_MWAS/Rdata_1216/annotation_graphlan.tsv",sep = "\t", col.names = F, row.names = F, quote = F)

