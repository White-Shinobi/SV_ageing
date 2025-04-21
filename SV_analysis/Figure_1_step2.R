# Pick out age-associated SVs in DMP and replicate in other cohorts
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part1_functions.R"))
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part2_functions.R"))

#### 1. read files ####
object_names <- c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")
directory <- "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/"
data_list <- list()
for (name in object_names) {
  file_path <- paste0(directory, "filter_zero_", name, ".rds")
  data_list[[name]] <- readRDS(file_path)
}

#### 2. Age-SVs in DMP ####
#### 2.1 all SVs: Age-SVs FDR threshold ####
filter_zero_dmp <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/filter_zero_dmp.rds")
non_NA_final <- readRDS("~/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
number_SV=length(non_NA_final$SV_Name[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)])
fdr_threshold_age=0.05/number_SV
saveRDS(fdr_threshold_age, paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/","fdr_threshold_age",".rds"))

#### 3. Replication in other 4 cohorts ####
#### 3.1 merge all files ####
# Reorder the list to make sure 'dmp' is the first element
dmp_age=na.omit(filter_zero_dmp)#7744
dmp_age=dmp_age[which(dmp_age$Y%in%non_NA_final$SV_Name[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]),]#6447
dmp_age_sig=dmp_age[dmp_age$p<fdr_threshold_age,]
data_list[["dmp"]]=dmp_age_sig
data_list <- data_list[c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")]
for(i in seq_along(data_list)) {
  # Rename columns except for the first two (assuming they are "X" and "Y")
  colnames(data_list[[i]]) <- ifelse(1:ncol(data_list[[i]]) <= 2, 
                                     colnames(data_list[[i]]), 
                                     paste(colnames(data_list[[i]]), names(data_list)[i], sep="_"))
}
merged_data <- Reduce(function(x, y) merge(x, y, by = c("X", "Y"), all.x = TRUE), data_list)
merged_data$Beta_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$SE_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$SE_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_lld1[which(is.na(merged_data$p_lld1))]=NA
merged_data$SE_lld1[which(is.na(merged_data$p_lld1))]=NA

#### 3.2 replication ####
replication_result=data.frame(SV=NA,p_other_cohorts=NA,"300ob"=NA,"300tzfg"=NA,"500fg_fsk"=NA,"lld1"=NA)
for (i in 1:nrow(merged_data)){
  print(i)
  replication_result[i,"SV"]=merged_data$Y[i]
  for (j in c("300ob","300tzfg","500fg_fsk","lld1")){
    print(j)
    # j=c("300ob")
    Beta_name=paste0("Beta_",j)
    P_name=paste0("p_",j)
    if (!j==c("lld1")){j=paste0("X",j)}
    if (is.na(merged_data[i,P_name])){merged_data[i,P_name]=1}
    if (merged_data[i,Beta_name]*merged_data[i,"Beta_dmp"]>0&merged_data[i,P_name]<0.05){replication_result[i,j]=1}else{replication_result[i,j]=0}
  }
}
replication_result$p_other_cohorts=replication_result$X300ob+replication_result$X300tzfg+replication_result$X500fg_fsk+replication_result$lld1
dmp_age_sig=dmp_age_sig[dmp_age_sig$Y%in%replication_result$SV[replication_result$p_other_cohorts>0],]
length(unique(dmp_age_sig$species))#16
dsgv_info_anno_ld <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
dmp_age_sig$SV_type=NA
dmp_age_sig$SV_type[dmp_age_sig$Y%in%dsgv_info_anno_ld$SV_Name]=c("dSV")
dmp_age_sig$SV_type[dmp_age_sig$Y%in%vsgv_info_anno_ld$SV_Name]=c("vSV")
saveRDS(dmp_age_sig, paste0("~/Documents/SV_MWAS/R/Raw_data/","meta_age_significant",".rds"))
#### 3.3 Table S4. Replication of DMP age-SV associations in the other four cohorts. ####
saveRDS(replication_result, paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/Results/","replication_result",".rds"))
write.table(replication_result, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
# Old data: ~/Documents/SV_MWAS/R/Raw_data/meta_age_significant_old_threshold.rds: the result using old threshold: SV cluster ()
# Old data: /Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/Results/replication_result_SV_cluster_threshold.rds: the result using old threshold: SV cluster ()

#### 4. Table S3. Association result between host age and SVs in DMP cohort (replication is done in the other four cohorts: 300OB, 300TZFG, 500FG and LLD). ####
object_names <- c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")
directory <- "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/"
data_list <- list()
for (name in object_names) {
  file_path <- paste0(directory, "filter_zero_", name, ".rds")
  data_list[[name]] <- readRDS(file_path)
}
data_list[["dmp"]]=data_list[["dmp"]][data_list[["dmp"]]$X==c("Age"),]
data_list <- data_list[c("dmp", "300ob", "300tzfg", "500fg_fsk", "lld1")]
for(i in seq_along(data_list)) {
  # Rename columns except for the first two (assuming they are "X" and "Y")
  colnames(data_list[[i]]) <- ifelse(1:ncol(data_list[[i]]) <= 2, 
                                     colnames(data_list[[i]]), 
                                     paste(colnames(data_list[[i]]), names(data_list)[i], sep="_"))
}
merged_data <- Reduce(function(x, y) merge(x, y, by = c("X", "Y"), all.x = TRUE), data_list)
merged_data=merged_data[!is.na(merged_data$p_dmp),]
merged_data$Replicable=NA
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
merged_data$Replicable[merged_data$Y%in%meta_age_significant$Y]=c("Replicable in at least one extra cohort")
merged_data <- merged_data[, c("Replicable", names(merged_data)[-which(names(merged_data) == "Replicable")])]  
merged_data$Beta_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$SE_300ob[which(is.na(merged_data$p_300ob))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$SE_500fg_fsk[which(is.na(merged_data$p_500fg_fsk))]=NA
merged_data$Beta_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$SE_300tzfg[which(is.na(merged_data$p_300tzfg))]=NA
merged_data$Beta_lld1[which(is.na(merged_data$p_lld1))]=NA
merged_data$SE_lld1[which(is.na(merged_data$p_lld1))]=NA
# add the dSV/vSV information
dsgv_info_anno_ld <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_dsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
vsgv_info_anno_ld <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_vsgv_info_anno_ld.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) 
merged_data$SV_type=NA
merged_data$SV_type[merged_data$Y%in%dsgv_info_anno_ld$SV_Name]=c("dSV")
merged_data$SV_type[merged_data$Y%in%vsgv_info_anno_ld$SV_Name]=c("vSV")
non_NA_final <- readRDS("~/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
merged_data=merged_data[which(merged_data$Y%in%non_NA_final$SV_Name[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]),]
dim(merged_data)#6447
for (i in c("dmp", "300ob", "300tzfg", "500fg_fsk")){
  #i=c("dmp")
  print(i)
  cohort_name=paste0("Cohort","_",i)
  print(cohort_name)
  merged_data[,grep(cohort_name,colnames(merged_data))]=toupper(i)
}
for (i in c("lld1")){
  #i=c("lld1")
  print(i)
  cohort_name=paste0("Cohort","_",i)
  print(cohort_name)
  merged_data[,grep(cohort_name,colnames(merged_data))]=c("LLD")
}
for (i in c("dmp", "300ob", "300tzfg", "500fg_fsk")){
  #i=c("dmp")
  print(i)
  colnames(merged_data)=gsub(i,toupper(i),colnames(merged_data))
}
for (i in c("lld1")){
  #i=c("dmp")
  print(i)
  colnames(merged_data)=gsub(i,"LLD",colnames(merged_data))
}
colnames(merged_data)=gsub("\\_FSK","",colnames(merged_data))
merged_data$Cohort_500FG=c("500FG")
write.table(merged_data, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 5.1 Figure 1A ####
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
tmp <- as.data.frame(table(meta_age_significant[,c("species"),drop=F]));
tmp=aggregate(Freq ~ species, data = tmp, sum);tmp=tmp[order(tmp$Freq,decreasing = F),]
tmp=tmp[tmp$Freq>0,]
tmp$species=factor(tmp$species,levels = tmp$species)
ggplot(tmp,aes(x=species,y=Freq))+
  geom_bar(stat='identity',position="stack",fill="#219ebc")+
  #geom_text(aes(label=Freq, vjust=0.5,hjust=-0.5)) +
  coord_flip()+
  theme(axis.text.y = element_text(face="plain",size=12,colour = "black"),
        axis.text.x = element_text(face="plain",size=16,colour = "black"),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16))+
  theme_classic()
  
#### 5.2 Correlation of Beta from different cohorts ####
file=merged_data[which(merged_data$Replicable==c("Replicable in at least one extra cohort")),c(grep("Beta",colnames(merged_data)),3)]
cor_result=data.frame(pheno=NA,rho=NA,p=NA,n=NA)
for (i in c("300ob", "300tzfg", "500fg_fsk", "lld1")){
  print(i)
  SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
  test=file[file$Y%in%SV_name,c(1,grep(i,colnames(file)))]%>%na.omit(.)
  cor_test <- cor.test(test[,1], test[,2], method = "spearman")
  cor_result[i,"pheno"]=i
  cor_result[i,"rho"]=cor_test$estimate
  cor_result[i,"p"]=cor_test$p.value
  cor_result[i,"n"]=nrow(test)
}
write.table(cor_result, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 5.3 Correlation of SE from different cohorts ####
replication_result=readRDS("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/Results/replication_result.rds")
file=merged_data[which(merged_data$Replicable==c("Replicable in at least one extra cohort")),c(grep("SE_",colnames(merged_data)),grep("Beta_",colnames(merged_data)),3)]
cor_result=data.frame(pheno=NA,rho=NA,p=NA,n=NA)
for (i in c("300ob", "300tzfg", "500fg_fsk", "lld1")){
  #i=c("300ob")
  print(i)
  SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
  test=file[file$Y%in%SV_name,c(1,grep(i,colnames(file)))]%>%na.omit(.)
  cor_test <- cor.test(test[,1], test[,2], method = "spearman")
  cor_result[i,"pheno"]=i
  cor_result[i,"rho"]=cor_test$estimate
  cor_result[i,"p"]=cor_test$p.value
  cor_result[i,"n"]=nrow(test)
}
write.table(cor_result, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 5.4 Figure 1C ####
i=c("500fg_fsk")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=file[file$Y%in%SV_name,c(grep("dmp",colnames(file)),grep(i,colnames(file)))]%>%na.omit(.)
test$Beta_dmp_min=test$Beta_dmp-1.96*test$SE_dmp
test$Beta_dmp_max=test$Beta_dmp+1.96*test$SE_dmp
test$Beta_500fg_fsk_min=test$Beta_500fg_fsk-1.96*test$SE_500fg_fsk
test$Beta_500fg_fsk_max=test$Beta_500fg_fsk+1.96*test$SE_500fg_fsk
cor_test <- cor.test(test[,"Beta_dmp"], test[,"Beta_500fg_fsk"], method = "spearman")
p=ggplot(data = test, mapping = aes(x =`Beta_dmp`, y =`Beta_500fg_fsk`))+
  labs(x = "Effect size in DMP cohort", y = "Effect size in 500FG cohort")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14))+
  #coord_fixed(( max(test[,1])-min(test[,1]))/(max(test[,2])-min(test[,2])))+
  geom_vline(xintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) + 
  geom_errorbar(aes(ymin = Beta_500fg_fsk_min,ymax = Beta_500fg_fsk_max),color="#DFD345") + 
  geom_errorbarh(aes(xmin = Beta_dmp_min,xmax = Beta_dmp_max),color="#DFD345")+
  geom_point(color="#62A9C4",size=2)+
  geom_text(x = 8, y = -5, label = paste0("R=",round(cor_test$estimate,2),", P"," = ", round(cor_test$p.value,3)), color = "black", size = 5)+   # Add text annotation for x = 0
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)
p_500fg=p+scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)+xlim(-1,1)+ylim(-1,1)
#### 5.5 Figure 1D ####
i=c("lld1")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=file[file$Y%in%SV_name,c(grep("dmp",colnames(file)),grep(i,colnames(file)))]%>%na.omit(.)
test$Beta_dmp_min=test$Beta_dmp-1.96*test$SE_dmp
test$Beta_dmp_max=test$Beta_dmp+1.96*test$SE_dmp
test$Beta_lld1_min=test$Beta_lld1-1.96*test$SE_lld1
test$Beta_lld1_max=test$Beta_lld1+1.96*test$SE_lld1
cor_test <- cor.test(test[,"Beta_dmp"], test[,"Beta_lld1"], method = "spearman")
p=ggplot(data = test, mapping = aes(x =`Beta_dmp`, y =`Beta_lld1`))+
  labs(x = "Effect size in DMP cohort", y = "Effect size in LLD cohort")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14))+
  #coord_fixed(( max(test[,1])-min(test[,1]))/(max(test[,2])-min(test[,2])))+
  geom_vline(xintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) + 
  geom_errorbar(aes(ymin = Beta_lld1_min,ymax = Beta_lld1_max),color="#DFD345") + 
  geom_errorbarh(aes(xmin = Beta_dmp_min,xmax = Beta_dmp_max),color="#DFD345")+
  geom_point(color="#62A9C4",size=2)+
  geom_text(x = 8, y = -5, label = paste0("R=",round(cor_test$estimate,2),", P"," = ", round(cor_test$p.value,3)), color = "black", size = 5)  # Add text annotation for x = 0
p_lld1=p+scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)+xlim(-1,1)+ylim(-1,1)

#### 5.4 Figure 1E ####
i=c("300ob")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=file[file$Y%in%SV_name,c(grep("dmp",colnames(file)),grep(i,colnames(file)))]%>%na.omit(.)
test$Beta_dmp_min=test$Beta_dmp-1.96*test$SE_dmp
test$Beta_dmp_max=test$Beta_dmp+1.96*test$SE_dmp
test$Beta_300ob_min=test$Beta_300ob-1.96*test$SE_300ob
test$Beta_300ob_max=test$Beta_300ob+1.96*test$SE_300ob
cor_test <- cor.test(test[,"Beta_dmp"], test[,"Beta_300ob"], method = "spearman")
p=ggplot(data = test, mapping = aes(x =`Beta_dmp`, y =`Beta_300ob`))+
  labs(x = "Effect size in DMP cohort", y = "Effect size in 300OB cohort")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14))+
  #coord_fixed(( max(test[,1])-min(test[,1]))/(max(test[,2])-min(test[,2])))+
  geom_vline(xintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) + 
  geom_errorbar(aes(ymin = Beta_300ob_min,ymax = Beta_300ob_max),color="#DFD345") + 
  geom_errorbarh(aes(xmin = Beta_dmp_min,xmax = Beta_dmp_max),color="#DFD345")+
  geom_point(color="#62A9C4",size=2)+
  geom_text(x = 8, y = -5, label = paste0("R=",round(cor_test$estimate,2),", P"," = ", round(cor_test$p.value,3)), color = "black", size = 5)+   # Add text annotation for x = 0
  scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)
p_300OB=p+scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)+xlim(-1,1)+ylim(-1,1)
#### 5.5 Figure 1F ####
i=c("300tzfg")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=file[file$Y%in%SV_name,c(grep("dmp",colnames(file)),grep(i,colnames(file)))]%>%na.omit(.)
test$Beta_dmp_min=test$Beta_dmp-1.96*test$SE_dmp
test$Beta_dmp_max=test$Beta_dmp+1.96*test$SE_dmp
test$Beta_300tzfg_min=test$Beta_300tzfg-1.96*test$SE_300tzfg
test$Beta_300tzfg_max=test$Beta_300tzfg+1.96*test$SE_300tzfg
cor_test <- cor.test(test[,"Beta_dmp"], test[,"Beta_300tzfg"], method = "spearman")
p=ggplot(data = test, mapping = aes(x =`Beta_dmp`, y =`Beta_300tzfg`))+
  labs(x = "Effect size in DMP cohort", y = "Effect size in 300TZFG cohort")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14))+
  #coord_fixed(( max(test[,1])-min(test[,1]))/(max(test[,2])-min(test[,2])))+
  geom_vline(xintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "#A52B28", size = 1) + 
  geom_errorbar(aes(ymin = Beta_300tzfg_min,ymax = Beta_300tzfg_max),color="#DFD345") + 
  geom_errorbarh(aes(xmin = Beta_dmp_min,xmax = Beta_dmp_max),color="#DFD345")+
  geom_point(color="#62A9C4",size=2)+
  geom_text(x = 8, y = -5, label = paste0("R=",round(cor_test$estimate,2),", P"," = ", round(cor_test$p.value,3)), color = "black", size = 5)  # Add text annotation for x = 0
p_300TZFG=p+scale_x_continuous(n.breaks = 5) +
  scale_y_continuous(n.breaks = 5)+xlim(-1,1)+ylim(-1,1)
#### 5.6 Figure 1C-F ####
library(patchwork)
p_500fg+p_lld1|(p_300OB+p_300TZFG)
#### 6. Meta-analysis for 300tzfg and 3000OB ####
i=c("300tzfg")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
file=merged_data[which(merged_data$Replicable==c("Replicable in at least one extra cohort")),c(grep("Beta",colnames(merged_data)),3)]
file=file[file$Y%in%SV_name,]
#test=merged_data[merged_data$Y%in%SV_name,c("Beta_dmp","SE_dmp","Beta_300tzfg","SE_300tzfg","Y","X")]%>%na.omit(.)
test=merged_data[merged_data$Y%in%file$Y,c("Beta_dmp","SE_dmp","Beta_300tzfg","SE_300tzfg","Y","X")]%>%na.omit(.)
meta_300TZFG=my_batch_meta_lm(test,c("DMP", "300TZFG"),c(1,3),c(2,4),row_var_col = 5,col_var_col = 6)
meta_300TZFG=meta_300TZFG[["table"]]
write.table(meta_300TZFG, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

i=c("300ob")
SV_name=replication_result$SV[replication_result[,grep(i,colnames(replication_result))]==1]
test=merged_data[merged_data$Y%in%SV_name,c("Beta_dmp","SE_dmp","Beta_300ob","SE_300ob","Y","X")]%>%na.omit(.)
meta_300OB=my_batch_meta_lm(test,c("DMP", "300OB"),c(1,3),c(2,4),row_var_col = 5,col_var_col = 6)
meta_300OB=meta_300OB[["table"]]



