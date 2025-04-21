source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part1_functions.R"))
#### 1.0 import data ####
setwd("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/")
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
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
#### 2.Numbers ####
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[which(!humann_results$significance==c("non-significant")),]
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results.rds")
UniRef90_not_list <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/UniRef90_not_list.rds")
Final_results=Final_results[which(Final_results$Y%in%UniRef90_not_list),]
unique(Final_results$Y)%>%length(.)#40280 genes
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)#39791 genes
length(unique(Final_results$Y))#39791
length(unique(Final_results$X))#112
Final_results=Final_results[which(Final_results$Y%in%humann_results_sig$Y),]
Final_results=Final_results[which(!Final_results$X==c("Age")),]
length(unique(Final_results$Y))#10486
length(unique(Final_results$X))#111
fdr_threshold=0.05/(10486*111)# 5.823535e-08
Final_results$Significance=NA
Final_results$Significance[which(Final_results$p_phenotype<fdr_threshold&Final_results$Beta_phenotype>0)]=c("significant-positive")
Final_results$Significance[which(Final_results$p_phenotype<fdr_threshold&Final_results$Beta_phenotype<0)]=c("significant-negative")
Final_results$Significance[which(Final_results$p_phenotype>fdr_threshold)]=c("non-significant")
Final_results$Significance=as.factor(Final_results$Significance)
Final_results$Significance=factor(Final_results$Significance,levels = c("significant-positive","significant-negative","non-significant"))
saveRDS(Final_results,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableS21",".rds"))
#Final_results too big, so just save the significant rows
write.table(Final_results[which(!Final_results$Significance==c("non-significant")),], "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

TableS21 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS21.rds")
TableS21_sig=TableS21[which(!TableS21$Significance==c("non-significant")),]
length(unique(TableS21_sig$Y))#4029
length(unique(TableS21_sig$X))#71 phenotypes
length(unique(TableS21_sig$Y))/10486

#### Figure S6 ####
# barplot for phenotype association number
tmp <- as.data.frame(table(TableS21_sig[,c("X"),drop=F]));
tmp=aggregate(Freq ~ X, data = tmp, sum);tmp=tmp[order(tmp$Freq,decreasing = F),]
tmp=tmp[tmp$Freq>0,]
tmp$X=factor(tmp$X,levels = tmp$X)
ggplot(tmp,aes(x=X,y=Freq))+
  geom_bar(stat='identity',position="stack",fill="#219ebc")+
  geom_text(aes(label=Freq, vjust=0.5,hjust=-0.5)) +
  coord_flip()+
  theme(axis.text.y = element_text(face="plain",size=8,colour = "black"),
        axis.text.x = element_text(face="plain",size=16,colour = "black"),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16))
#### Figure 4A ####
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")

tmp1=humann_results[which(!humann_results$significance==c("non-significant")),]
tmp2=tmp1$Y[tmp1$Beta_phenotype<0]
tmp3=tmp1$Y[tmp1$Beta_phenotype>0]
length(tmp2)#6242
length(tmp3)#4244

TableS21 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS21.rds")
for (i in c("Y","X")){`TableS21`[,i]=as.character(`TableS21`[,i])}
TableS21_sig=TableS21[which(!TableS21$Significance==c("non-significant")),]
TableS21_sig$X[which(TableS21_sig$X==c("FI41_FU"))]=c("FI41")

tmp4=TableS21_sig[TableS21_sig$X==c("FI41"),]
tmp5=tmp4$Y[tmp4$Beta_phenotype<0]
tmp6=tmp4$Y[tmp4$Beta_phenotype>0]
length(tmp5)#461
length(tmp6)#373

#Group1:716
#Group1a:344
intersect(tmp3,tmp5)%>%unique(.)%>%length(.)
one_a=intersect(tmp3,tmp5)%>%unique(.)
saveRDS(one_a, paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","one_a_humann",".rds"))
#Group1b:372
intersect(tmp2,tmp6)%>%unique(.)%>%length(.)
one_b=intersect(tmp2,tmp6)%>%unique(.)

#Group2:118
#Group2a:1
intersect(tmp3,tmp6)%>%unique(.)%>%length(.)
two_a=intersect(tmp3,tmp6)%>%unique(.)

#Group2b:117
intersect(tmp2,tmp5)%>%unique(.)%>%length(.)
two_b=intersect(tmp2,tmp5)%>%unique(.)


#### Table S25. The genes classified as the geroprotective pattern and the geropathogenic pattern (genes are extracted from non-age-associated SVs).  ####
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
TableS21 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS21.rds")

table1=humann_results[humann_results$Y%in%c(one_a,one_b),c(1:5)]
table1$type=c("geroprotective")
table2=humann_results[humann_results$Y%in%c(two_a,two_b),c(1:5)]
table2$type=c("geropathogenic")

table3=TableS21[TableS21$Y%in%c(one_a,one_b),c(1:5)]
table3$type=c("geroprotective")
table4=TableS21[TableS21$Y%in%c(two_a,two_b),c(1:5)]
table4$type=c("geropathogenic")

all_table=rbind(table1,table2,table3,table4)
all_table=all_table[all_table$X%in%c("Age","FI41_FU"),]

Uniref90_name <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/Uniref90_name.rds")
all_table=merge(all_table,Uniref90_name,by.x = "Y",by.y = "V1",all.x = T)
colnames(all_table)[7]=c("description")
saveRDS(all_table,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableS_gene_type_humann",".rds"))
write.table(all_table, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 4.3 GO enrichment analysis using clusterProfiler ####
library(clusterProfiler)
#### 4.3.1 geroprotective ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#24516

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[humann_results$Y%in%c(one_a,one_b),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#539
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
GO_reasult_geroprotective <- data.frame(go_rich@result)
saveRDS(GO_reasult_geroprotective, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_geroprotective.rds")
write.table(GO_reasult_geroprotective, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(go_rich)  #富集柱形图
dotplot(go_rich,  showCategory = 5)  #富集气泡图
cnetplot(go_rich) #网络图展示富集功能和基因的包含关系
emapplot(go_rich) #网络图展示各富集功能之间共有基因关系
heatplot(go_rich) #热图展示富集功能和基因的包含关系
#### 4.3.2 geropathogenic ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#24516

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[humann_results$Y%in%c(two_a,two_b),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#71
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.2, 
                    maxGSSize = 500)
GO_reasult_geropathogenic <- data.frame(go_rich@result)
barplot(go_rich)  #富集柱形图
dotplot(go_rich,  showCategory = 16)  #富集气泡图
cnetplot(go_rich) #网络图展示富集功能和基因的包含关系
emapplot(go_rich) #网络图展示各富集功能之间共有基因关系
heatplot(go_rich) #热图展示富集功能和基因的包含关系
saveRDS(GO_reasult_geropathogenic, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_positive.rds")
write.table(GO_reasult_geropathogenic, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 4.3.3 one_a+two_b ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#11693
11639/20058 #58.0%

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[humann_results$Y%in%c(one_a,two_b),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#247
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'BH', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
GO_reasult_geroprotective <- data.frame(go_rich@result)
saveRDS(GO_reasult_geroprotective, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_geroprotective.rds")
write.table(GO_reasult_negative, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(go_rich)  #富集柱形图
dotplot(go_rich,  showCategory = 3)  #富集气泡图
cnetplot(go_rich) #网络图展示富集功能和基因的包含关系
emapplot(go_rich) #网络图展示各富集功能之间共有基因关系
heatplot(go_rich) #热图展示富集功能和基因的包含关系

#### 4.3.4 208 LLD-replicated FI-associated genes ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#24516

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
rep_LLD <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_Humann_LLD.rds")
rep_LLD=rep_LLD[which(rep_LLD$replicated_in_LLD==c("Yes")),]
humann_results_sig=humann_results[humann_results$Y%in%rep_LLD$Y[rep_LLD$X==c("FI41")],]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#155
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
GO_reasult_208_Fi <- data.frame(go_rich@result)
saveRDS(GO_reasult_208_Fi, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_208_Fi.rds")
write.table(GO_reasult_208_Fi, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(go_rich)  #富集柱形图
dotplot(go_rich,  showCategory = 3)  #富集气泡图
cnetplot(go_rich) #网络图展示富集功能和基因的包含关系
emapplot(go_rich) #网络图展示各富集功能之间共有基因关系
heatplot(go_rich) #热图展示富集功能和基因的包含关系


#### 4.4 KO enrichment analysis using clusterProfiler ####
library(clusterProfiler)
#### 4.4.1 208 LLD-replicated FI-associated genes ####
#KO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Humann_Uniref90_to_KO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#19254
19254/39791 #48.0%

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[humann_results$Y%in%rep_LLD$Y[rep_LLD$X==c("FI41")],]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#116
# Enrichment analysis
KO_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("KO_number","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("KO_number","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
KO_reasult_208_FI <- data.frame(KO_rich@result)
saveRDS(KO_reasult_negative, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/KO_reasult_negative.rds")
write.table(KO_reasult_negative, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(KO_rich)  #富集柱形图






#### 4.4.2 geroprotective ####
#KO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Humann_Uniref90_to_KO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))#19254
19254/39791 #48.0%

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[humann_results$Y%in%c(one_a,one_b),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)#452
# Enrichment analysis
KO_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("KO_number","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("KO_number","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
KO_reasult_208_FI <- data.frame(KO_rich@result)
saveRDS(KO_reasult_negative, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/KO_reasult_negative.rds")
write.table(KO_reasult_negative, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(KO_rich)  #富集柱形图

#### Figure 4B : LLD replication ####
#DMP
TableS21 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS21.rds")
for (i in c("Y","X")){`TableS21`[,i]=as.character(`TableS21`[,i])}
TableS21_sig=TableS21[which(!TableS21$Significance==c("non-significant")),]
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
for (i in c("Y","X")){`humann_results`[,i]=as.character(`humann_results`[,i])}
humann_results_sig=humann_results[which(!humann_results$significance==c("non-significant")),]
TableS21_sig=TableS21_sig[,1:16]
humann_results_sig=humann_results_sig[,1:16]
colnames(humann_results_sig)
colnames(TableS21_sig)
colnames(humann_results_sig)[16]=c("Significance")
humann_results_sig=humann_results_sig[,colnames(TableS21_sig)]
humann_results=rbind(TableS21_sig,humann_results_sig)
length(unique(humann_results$Y))#10486
length(unique(humann_results$X))#72
humann_results$X[which(humann_results$X==c("FI41_FU"))]=c("FI41")

#LLD
Final_results_humann_LLD <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_humann_LLD.rds")
Final_results_humann_LLD$Beta_IDage=NULL;Final_results_humann_LLD$SE_IDage=NULL;Final_results_humann_LLD$p_IDage=NULL
for (i in c("Y","X")){`Final_results_humann_LLD`[,i]=as.character(`Final_results_humann_LLD`[,i])}
Final_results_humann_LLD$X[grep("FI41_",Final_results_humann_LLD$X)]=c("FI41")
length(unique(Final_results_humann_LLD$Y))#19793
length(unique(Final_results_humann_LLD$X))#253
replication_result=merge(humann_results, Final_results_humann_LLD, by = c("X", "Y"), all.x = TRUE)
colnames(replication_result)=gsub("\\.x","_DMP",colnames(replication_result))
colnames(replication_result)=gsub("\\.y","_LLD",colnames(replication_result))
replication_result$replicated_in_LLD=NA
replication_result$replicated_in_LLD[replication_result$p_phenotype_LLD<0.05&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0]=c("Yes")
replication_result$replicated_in_LLD[replication_result$p_phenotype_LLD>0.05]=c("No")
replication_result$replicated_in_LLD[replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD<0]=c("No")
length(which(replication_result$X==c("FI41")&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0))
477/834
length(which(replication_result$X==c("FI41")&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0&replication_result$p_phenotype_LLD<0.05))#
208/477
cor_test <- cor.test(replication_result[replication_result$X==c("FI41"),"Beta_phenotype_DMP"], replication_result[replication_result$X==c("FI41"),"Beta_phenotype_LLD"], method = "spearman")
p=ggplot(data = replication_result[replication_result$X==c("FI41"),], mapping = aes(x =`Beta_phenotype_DMP`, y =`Beta_phenotype_LLD`))+
  geom_point(aes(color=replicated_in_LLD),size=3)+
  labs(x = "Effect size in DMP", y = "Effect size in LLD")+
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=14),
        axis.text.x = element_text(face="plain",size=14),
        axis.title.y = element_text(face="bold",size=14),
        axis.title.x = element_text(face="bold",size=14))+
  #coord_fixed(( max(test[,1])-min(test[,1]))/(max(test[,2])-min(test[,2])))+
  geom_vline(xintercept = 0, linetype = "dotted", color = "#415a77", size = 1) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "#415a77", size = 1) +  
  scale_color_manual(values = c("Yes" = "#B78045", "No" = "#7D8C86"))+
  #geom_text(x = 0.1, y = -0.4, label = paste0("R=",round(cor_test$estimate,2),", P = "," 0.001"), color = "black", size = 5)+   # Add text annotation for x = 0
  theme(axis.text.y = element_text(face="plain",size=16,colour = "black"),
        axis.text.x = element_text(face="plain",size=16,colour = "black"),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16),)
p+xlim(-0.5,0.3)+ylim(-0.5,0.3)
#### Table S17. LLD replication results for the genes both associated with age and health-realted metrices discovered in DMP cohort. ####
saveRDS(replication_result,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","replication_result_Humann_LLD",".rds"))
write.table(replication_result, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 8.0 Forest plots (Figure 4 and E) ####
rep_LLD <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_Humann_LLD.rds")
rep_LLD=rep_LLD[which(rep_LLD$replicated_in_LLD==c("Yes")|rep_LLD$X==c("Age")),]
rep_LLD$group <- NA
rep_LLD$group[rep_LLD$Y%in%c(one_a,one_b)]=c("Geroprotective genes")
rep_LLD$group[rep_LLD$Y%in%c(two_a,two_b)]=c("Geropathogenic genes")
rep_LLD$group=as.factor(rep_LLD$group)
rep_LLD$conf.low=rep_LLD$Beta_phenotype_DMP-1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$conf.high=rep_LLD$Beta_phenotype_DMP+1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$model=paste0(rep_LLD$Y,"~",rep_LLD$X)
rep_LLD$Y=as.factor(rep_LLD$Y)
Uniref90_name <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/Uniref90_name.rds")
rep_LLD=merge(rep_LLD,Uniref90_name,by.x = "Y",by.y = "V1",all.x = T)

#geroprotective  
#plot=rep_LLD[rep_LLD$X%in%c("Age","FI41")&rep_LLD$Y%in%c("UniRef90_R5P1L6","UniRef90_D4JXP6","UniRef90_A0A2A7AAA8"),] #geropathogenic 
# UniRef90_R5P1L6: Putative signal-transducing histidine kinase
# UniRef90_D4JXP6: DNA or RNA helicases of superfamily II
# UniRef90_A0A2A7AAA8: DNA-binding protein
plot=rep_LLD[rep_LLD$Y%in%c("UniRef90_A0A174AJ31","UniRef90_D4S234"),]

plot$model <- as.factor(as.character(plot$model))
plot$Y <- as.factor(as.character(plot$Y))
plot$row_id <- as.numeric(plot$model)
library(wesanderson)
ggplot(plot,aes(x=fct_rev(model),y=Beta_phenotype_DMP,ymin=conf.low,ymax=conf.high)) +
  geom_rect(aes(xmin = row_id - 0.5, xmax = row_id + 0.5, ymin = -Inf, ymax = Inf), 
            fill = rep(c("gray", "white"), length.out = nrow(plot)), 
            alpha = 1, inherit.aes = FALSE)+
  geom_pointrange(aes(color = Y), size = .5) +
  scale_color_manual(values = wes_palette(2, name = "AsteroidCity2", type = "discrete"))+
  geom_hline(yintercept = 0, color = "steelblue") +  
  xlab(" ") +
  ylab("Coefficient (95% Confidence Interval)") +
  labs(title ="Linear Regression Models Estimating the Effects of Vehicle Weight \n on Fuel Efficiency") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))+
  coord_flip()
#facet_grid(Y ~ ., scales = "free", switch = "both")
#geropathogenic 
rep_LLD <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_Humann_LLD.rds")
rep_LLD$group <- NA
rep_LLD$group[rep_LLD$Y%in%c(one_a,one_b)]=c("Geroprotective genes")
rep_LLD$group[rep_LLD$Y%in%c(two_a,two_b)]=c("Geropathogenic genes")
rep_LLD$group=as.factor(rep_LLD$group)
rep_LLD$conf.low=rep_LLD$Beta_phenotype_DMP-1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$conf.high=rep_LLD$Beta_phenotype_DMP+1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$model=paste0(rep_LLD$Y,"~",rep_LLD$X)
rep_LLD$Y=as.factor(rep_LLD$Y)
plot=rep_LLD[rep_LLD$X%in%c("Age","FI41")&rep_LLD$Y%in%c("UniRef90_D4K3J1","UniRef90_A0A174GWM5","UniRef90_A0A2J4JK13"),]
#UniRef90_D4K3J1: DNA-directed RNA polymerase specialized sigma subunit, sigma24 homolog
#UniRef90_A0A174GWM5: D-gamma-glutamyl-meso-diaminopimelic acid endopeptidase CwlS
#UniRef90_A0A2J4JK13 :Histidine phosphatase family protein

plot$model <- as.factor(as.character(plot$model))
plot$Y <- as.factor(as.character(plot$Y))
plot$row_id <- as.numeric(plot$model)
library(wesanderson)
ggplot(plot,aes(x=fct_rev(model),y=Beta_phenotype_DMP,ymin=conf.low,ymax=conf.high,fill=Y)) +
  geom_rect(aes(xmin = row_id - 0.5, xmax = row_id + 0.5, ymin = -Inf, ymax = Inf), 
            fill = rep(c("gray", "white"), length.out = nrow(plot)), 
            alpha = 1, inherit.aes = FALSE)+
  geom_pointrange(aes(color = Y),  size = .5) +
  scale_color_manual(values = wes_palette(3, name = "Rushmore1", type = "discrete"))+
  geom_hline(yintercept = 0, color = "steelblue") +  
  xlab(" ") +
  ylab("Coefficient (95% Confidence Interval)") +
  labs(title ="Linear Regression Models Estimating the Effects of Vehicle Weight \n on Fuel Efficiency") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))+
  coord_flip()
#### 9.0 Enrichment analysis:Geroprotective genes ####
#### KO enrichment ####
Humann_gene_info <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Humann_gene_info.rds")
Humann_gene_info=na.omit(Humann_gene_info[,c("UniRef90","KO_number")])
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
Humann_gene_info$UniRef90%>%unique(.)%>%length(.)#1793

Geroprotective_genes_results=Humann_gene_info[Humann_gene_info$UniRef90%in%c(one_a,one_b),]
res.enrich=data.frame()
for (KO in unique(na.omit(Humann_gene_info$KO_number))){
  print(KO)
  #KO=unique(na.omit(Humann_gene_info$KO_number))[1]
  #KO=c("K06223")
  tmp=Geroprotective_genes_results[which(Geroprotective_genes_results$KO_number==KO),]
  # the table
  inputDf <- data.frame(Have_this_KO = c(NA, NA), Not_have_this_KO = c(NA, NA))
  row.names(inputDf)=c("Geroprotective genes","Background")
  inputDf$Have_this_KO[1]=dim(tmp)[1]
  inputDf$Not_have_this_KO[1]=nrow(Geroprotective_genes_results)-dim(tmp)[1]
  inputDf$Have_this_KO[2]=length(which(Humann_gene_info$KO_number==KO))
  inputDf$Not_have_this_KO[2]=length(unique(Humann_gene_info$UniRef90))-inputDf$Have_this_KO[2]
  if (sum(inputDf==0)>0){inputDf=inputDf+1}
  if (nrow(tmp)>0){
    fisher.res <- fisher.test(inputDf)
    res.enrich[KO,1]=fisher.res$estimate
    res.enrich[KO,2]=fisher.res$p.value
    res.enrich[KO,3]=inputDf$Have_this_KO[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[KO,4]=inputDf$Not_have_this_KO[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[KO,5]=inputDf$Have_this_KO[row.names(inputDf)==c("Background")]
    res.enrich[KO,6]=inputDf$Not_have_this_KO[row.names(inputDf)==c("Background")]
    res.enrich[KO,7]=KO
  }else{
    res.enrich[KO,1]=NA
    res.enrich[KO,2]=NA
    res.enrich[KO,3]=inputDf$Have_this_KO[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[KO,4]=inputDf$Not_have_this_KO[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[KO,5]=inputDf$Have_this_KO[row.names(inputDf)==c("Background")]
    res.enrich[KO,6]=inputDf$Not_have_this_KO[row.names(inputDf)==c("Background")]
    res.enrich[KO,7]=KO
  }
}
colnames(res.enrich)=c("OddsRatio","p","Have_this_KO_interest","Not_have_this_KO_interest","Have_this_KO_B","Not_have_this_KO_B","KO")
res.enrich$fdr.BH=p.adjust(res.enrich$p,method = "BH")

#### Module enrichment ####
Humann_gene_info <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Humann_gene_info.rds")
Humann_gene_info=na.omit(Humann_gene_info[,c("UniRef90","Module")])
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
Humann_gene_info$UniRef90%>%unique(.)%>%length(.)#1787

Geroprotective_genes_results=Humann_gene_info[Humann_gene_info$UniRef90%in%c(one_a,one_b),]
res.enrich=data.frame()
for (Module in unique(na.omit(Humann_gene_info$Module))){
  print(Module)
  #Module=unique(na.omit(Humann_gene_info$Module))[1]
  #Module=c("K23779")
  tmp=Geroprotective_genes_results[which(Geroprotective_genes_results$Module==Module),]
  # the table
  inputDf <- data.frame(Have_this_Module = c(NA, NA), Not_have_this_Module = c(NA, NA))
  row.names(inputDf)=c("Geroprotective genes","Background")
  inputDf$Have_this_Module[1]=dim(tmp)[1]
  inputDf$Not_have_this_Module[1]=nrow(Geroprotective_genes_results)-dim(tmp)[1]
  inputDf$Have_this_Module[2]=length(which(Humann_gene_info$Module==Module))
  inputDf$Not_have_this_Module[2]=length(unique(Humann_gene_info$UniRef90))-inputDf$Have_this_Module[2]
  if (sum(inputDf==0)>0){inputDf=inputDf+1}
  if (nrow(tmp)>0){
    fisher.res <- fisher.test(inputDf)
    res.enrich[Module,1]=fisher.res$estimate
    res.enrich[Module,2]=fisher.res$p.value
    res.enrich[Module,3]=inputDf$Have_this_Module[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[Module,4]=inputDf$Not_have_this_Module[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[Module,5]=inputDf$Have_this_Module[row.names(inputDf)==c("Background")]
    res.enrich[Module,6]=inputDf$Not_have_this_Module[row.names(inputDf)==c("Background")]
    res.enrich[Module,7]=Module
  }else{
    res.enrich[Module,1]=NA
    res.enrich[Module,2]=NA
    res.enrich[Module,3]=inputDf$Have_this_Module[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[Module,4]=inputDf$Not_have_this_Module[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[Module,5]=inputDf$Have_this_Module[row.names(inputDf)==c("Background")]
    res.enrich[Module,6]=inputDf$Not_have_this_Module[row.names(inputDf)==c("Background")]
    res.enrich[Module,7]=Module
  }
}
colnames(res.enrich)=c("OddsRatio","p","Have_this_Module_interest","Not_have_this_Module_interest","Have_this_Module_B","Not_have_this_Module_B","Module")
res.enrich$fdr.BH=p.adjust(res.enrich$p,method = "BH")
#### RXN enrichment ####
Humann_gene_info <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Humann_gene_info.rds")
Humann_gene_info=na.omit(Humann_gene_info[,c("UniRef90","RXN")])
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
Humann_gene_info$UniRef90%>%unique(.)%>%length(.)#2656

Geroprotective_genes_results=Humann_gene_info[Humann_gene_info$UniRef90%in%c(one_a,one_b),]
res.enrich=data.frame()
for (RXN in unique(na.omit(Humann_gene_info$RXN))){
  print(RXN)
  #RXN=unique(na.omit(Humann_gene_info$RXN))[1]
  #RXN=c("K23779")
  tmp=Geroprotective_genes_results[which(Geroprotective_genes_results$RXN==RXN),]
  # the table
  inputDf <- data.frame(Have_this_RXN = c(NA, NA), Not_have_this_RXN = c(NA, NA))
  row.names(inputDf)=c("Geroprotective genes","Background")
  inputDf$Have_this_RXN[1]=dim(tmp)[1]
  inputDf$Not_have_this_RXN[1]=nrow(Geroprotective_genes_results)-dim(tmp)[1]
  inputDf$Have_this_RXN[2]=length(which(Humann_gene_info$RXN==RXN))
  inputDf$Not_have_this_RXN[2]=length(unique(Humann_gene_info$UniRef90))-inputDf$Have_this_RXN[2]
  if (sum(inputDf==0)>0){inputDf=inputDf+1}
  if (nrow(tmp)>0){
    fisher.res <- fisher.test(inputDf)
    res.enrich[RXN,1]=fisher.res$estimate
    res.enrich[RXN,2]=fisher.res$p.value
    res.enrich[RXN,3]=inputDf$Have_this_RXN[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[RXN,4]=inputDf$Not_have_this_RXN[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[RXN,5]=inputDf$Have_this_RXN[row.names(inputDf)==c("Background")]
    res.enrich[RXN,6]=inputDf$Not_have_this_RXN[row.names(inputDf)==c("Background")]
    res.enrich[RXN,7]=RXN
  }else{
    res.enrich[RXN,1]=NA
    res.enrich[RXN,2]=NA
    res.enrich[RXN,3]=inputDf$Have_this_RXN[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[RXN,4]=inputDf$Not_have_this_RXN[row.names(inputDf)==c("Geroprotective genes")]
    res.enrich[RXN,5]=inputDf$Have_this_RXN[row.names(inputDf)==c("Background")]
    res.enrich[RXN,6]=inputDf$Not_have_this_RXN[row.names(inputDf)==c("Background")]
    res.enrich[RXN,7]=RXN
  }
}
colnames(res.enrich)=c("OddsRatio","p","Have_this_RXN_interest","Not_have_this_RXN_interest","Have_this_RXN_B","Not_have_this_RXN_B","RXN")
res.enrich$fdr.BH=p.adjust(res.enrich$p,method = "BH")

#### PWY enrichment ####
Humann_gene_info <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Humann_gene_info.rds")
Humann_gene_info=na.omit(Humann_gene_info[,c("UniRef90","PWY")])
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
Humann_gene_info$UniRef90%>%unique(.)%>%length(.)#1516

Geroprotective_genes_results=Humann_gene_info[Humann_gene_info$UniRef90%in%c(one_a,one_b),]
res.enrich=data.frame()
for (PWY in unique(na.omit(Humann_gene_info$PWY))){
  print(PWY)
  #PWY=unique(na.omit(Humann_gene_info$PWY))[1]
  #PWY=c("K23779")
  tmp=Geroprotective_genes_results[which(Geroprotective_genes_results$PWY==PWY),]
  # the table
  inputDf <- data.frame(Have_this_PWY = c(NA, NA), Not_have_this_PWY = c(NA, NA))
  row.names(inputDf)=c("significant-negative","Background")
  inputDf$Have_this_PWY[1]=dim(tmp)[1]
  inputDf$Not_have_this_PWY[1]=nrow(Geroprotective_genes_results)-dim(tmp)[1]
  inputDf$Have_this_PWY[2]=length(which(Humann_gene_info$PWY==PWY))
  inputDf$Not_have_this_PWY[2]=length(unique(Humann_gene_info$UniRef90))-inputDf$Have_this_PWY[2]
  if (sum(inputDf==0)>0){inputDf=inputDf+1}
  if (nrow(tmp)>0){
    fisher.res <- fisher.test(inputDf)
    res.enrich[PWY,1]=fisher.res$estimate
    res.enrich[PWY,2]=fisher.res$p.value
    res.enrich[PWY,3]=inputDf$Have_this_PWY[row.names(inputDf)==c("significant-negative")]
    res.enrich[PWY,4]=inputDf$Not_have_this_PWY[row.names(inputDf)==c("significant-negative")]
    res.enrich[PWY,5]=inputDf$Have_this_PWY[row.names(inputDf)==c("Background")]
    res.enrich[PWY,6]=inputDf$Not_have_this_PWY[row.names(inputDf)==c("Background")]
    res.enrich[PWY,7]=PWY
  }else{
    res.enrich[PWY,1]=NA
    res.enrich[PWY,2]=NA
    res.enrich[PWY,3]=inputDf$Have_this_PWY[row.names(inputDf)==c("significant-negative")]
    res.enrich[PWY,4]=inputDf$Not_have_this_PWY[row.names(inputDf)==c("significant-negative")]
    res.enrich[PWY,5]=inputDf$Have_this_PWY[row.names(inputDf)==c("Background")]
    res.enrich[PWY,6]=inputDf$Not_have_this_PWY[row.names(inputDf)==c("Background")]
    res.enrich[PWY,7]=PWY
  }
}
colnames(res.enrich)=c("OddsRatio","p","Have_this_PWY_interest","Not_have_this_PWY_interest","Have_this_PWY_B","Not_have_this_PWY_B","PWY")
res.enrich$fdr.BH=p.adjust(res.enrich$p,method = "bonferroni")
PWY_name <- read.delim("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/map_metacyc-pwy_name.txt",
                       sep="\t", row.names=NULL, header=FALSE, check.names=FALSE, fill=TRUE, quote="")
res.enrich=merge(res.enrich,PWY_name,by.x = "PWY",by.y = "V1",all.x = T)

