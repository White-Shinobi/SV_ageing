source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part1_functions.R"))
# old version in NM_revision folder
#### 1.Humann ~ age and phenotypes association results ####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results.rds")
unique(Final_results$Y)%>%length(.)#40565 gene families
UniRef90_not_list <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/UniRef90_not_list.rds")
length(UniRef90_not_list)#69966
Final_results=Final_results[which(Final_results$Y%in%UniRef90_not_list),]
unique(Final_results$Y)%>%length(.)#40,280 genes
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)#39791 genes
fdr_threshold=0.05/(length(unique(Final_results$Y))*1)#1.256566e-06
Final_results$significance=NA
Final_results$significance[Final_results$p_phenotype<fdr_threshold&Final_results$Beta_phenotype>0]=c("significant-positive")
Final_results$significance[Final_results$p_phenotype<fdr_threshold&Final_results$Beta_phenotype<0]=c("significant-negative")
Final_results$significance[Final_results$p_phenotype>fdr_threshold]=c("non-significant")
Final_results$significance=as.factor(Final_results$significance)
Final_results$significance=factor(Final_results$significance,levels = c("significant-positive","significant-negative","non-significant"))
#Uniref90_name <- read.delim("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/map_uniref90_name.txt",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
#saveRDS(Uniref90_name, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Uniref90_name.rds")
Uniref90_name <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/Uniref90_name.rds")
Final_results=merge(Final_results,Uniref90_name,by.x = "Y",by.y = "V1",all.x = T)
for (i in c("Y","X")){`Final_results`[,i]=as.character(`Final_results`[,i])}
colnames(Final_results)[grep("V2",colnames(Final_results))]=c("UniRef90_description")
Final_results$UniRef90_description[is.na(Final_results$UniRef90_description)]=c("Unknown function")
#### 1.1 Table SX ####
Final_results=Final_results[Final_results$X==c("Age"),]
saveRDS(Final_results,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableSX1",".rds"))
write.table(Final_results, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#### 2. volcano plot ####
library(ggplot2) ##绘图使用
library(ggprism) ##设置主题私用
library(ggrepel) ##给火山图加标签使用
ggplot(Final_results[Final_results$X==c("Age"),], aes(x =Beta_phenotype, y= -log10(p_phenotype), color=significance)) +
  geom_point(alpha=1, size=2.5) +
  scale_color_manual(values=c('#E86D6C','#1D401F','gray')) +
  xlim(c(min(Final_results$Beta_phenotype), max(Final_results$Beta_phenotype))) +  ##调整x轴的取值范围，可以根据max(abs(Humann4.0_DMP.rds_with_zero$Beta))，获得差异基因最大值是多少，再进行取舍
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05/(20058*1)), lty=4,col="black",lwd=0.8) + 
  labs(x="Beta effect", y="-log10(P)") +
  scale_x_continuous(limits = c(-0.8, 0.8),n.breaks = 8) +
  ggtitle("Replication in LLD cohort") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)
#### 2.2 Numbers ####
sum(Final_results$significance==c("significant-negative"))#6242
sum(Final_results$significance==c("significant-positive"))#4244
(6242+4244)/39791

#### 3. Humann gene annotation file ####
#### Uniref90 -> GO ####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results.rds")
UniRef90_not_list <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/UniRef90_not_list.rds")
Final_results=Final_results[which(Final_results$Y%in%UniRef90_not_list),]
unique(Final_results$Y)%>%length(.)
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)
Humann_gene_info=data.frame(UniRef90=unique(Final_results$Y))
Unire90_GO <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/Unire90_GO.rds")
Humann_gene_info=merge(Humann_gene_info,Unire90_GO,by.x = "UniRef90",by.y = "Uniref90",all.x = T)
Humann_gene_info=na.omit(Humann_gene_info)
length(unique(Humann_gene_info$UniRef90))
GO_name <- read.delim("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/map_go_name.txt",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
Humann_gene_info=merge(Humann_gene_info,GO_name,by.x = "GO",by.y = "V1",all.x = T)
saveRDS(Humann_gene_info, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")

#### 3.3 Table SX ####
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
write.table(Humann_gene_info, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
length(unique(Humann_gene_info$UniRef90))

#### 4. Enrichment analysis ####
#### 4.1 GO enrichment analysis using clusterProfiler ####
library(clusterProfiler)
#### 4.1.1 significant-negative genes ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[which(humann_results$significance==c("significant-negative")),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.05, 
                    maxGSSize = 500)
GO_reasult_negative <- data.frame(go_rich@result)
saveRDS(GO_reasult_negative, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_negative.rds")
write.table(GO_reasult_negative, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
barplot(go_rich)  
dotplot(go_rich,  showCategory = 6)  
cnetplot(go_rich) 
emapplot(go_rich) 
heatplot(go_rich) 
#### 4.1.2 significant-positive genes ####
#GO annotation
Humann_gene_info <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/Humann_gene_info_GO.rds")
Humann_gene_info <- Humann_gene_info[!duplicated(Humann_gene_info), ]
length(unique(Humann_gene_info$UniRef90))

# interest gene list
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
humann_results_sig=humann_results[which(humann_results$significance==c("significant-positive")),]
humann_results_sig=humann_results_sig[humann_results_sig$Y%in%Humann_gene_info$UniRef90,]
gene_list<- humann_results_sig$Y
length(gene_list)
# Enrichment analysis
go_rich <- enricher(gene = gene_list,
                    TERM2GENE = Humann_gene_info[,c("GO","UniRef90")], 
                    TERM2NAME = Humann_gene_info[,c("GO","V2")], 
                    pvalueCutoff = 1, 
                    pAdjustMethod = 'bonferroni', 
                    qvalueCutoff = 0.2, 
                    maxGSSize = 500)
GO_reasult_positive <- data.frame(go_rich@result)
barplot(go_rich, showCategory = 21)  
dotplot(go_rich)  
cnetplot(go_rich) 
emapplot(go_rich) 
heatplot(go_rich) 
saveRDS(GO_reasult_positive, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/GO_reasult_positive.rds")
write.table(GO_reasult_positive, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 5. Figure 2C and D ####
Humann4.0_DMP <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Humann4.0_DMP.rds")
sample_number <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/Running_list/sampleCohorts.id.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
sample_number$Cohort_2=NA
sample_number$Cohort_2[sample_number$Cohort==c("300OB")]=c("300ob")
sample_number$Cohort_2[sample_number$Cohort==c("300TZFG")]=c("300tzfg")
sample_number$Cohort_2[sample_number$Cohort==c("500FG_FSK")]=c("500fg_fsk")
sample_number$Cohort_2[sample_number$Cohort==c("LLD1")]=c("lld1")
sample_number$Cohort_2[sample_number$Cohort==c("DMP")]=c("dmp")
load("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
load("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/readcounts_rawdata/full_read.RData")
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration", "Age", "Sex","BristolType")], full_read[,c("V1", "log10_counts"), drop=FALSE], by.x = "row.names", by.y = "V1", all=FALSE) %>% `rownames<-`(.[, "Row.names"]) %>% dplyr::select(-"Row.names")
humann_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableSX1.rds")
#### Negative ####
#### 5.1 UniRef90_A1A399 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_A1A399"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p1=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#1D3F21")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_A1A399\nBeta-galactosidase BgaB")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)

#### 5.2 UniRef90_A0A076JIZ1 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_A0A076JIZ1"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p2=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#1D3F21")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_A0A076JIZ1\nLacI family transcriptional regulator")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)


#### 5.3 UniRef90_A0A174CI30 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_A0A174CI30"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p3=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#1D3F21")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_A0A174CI30\nBeta-glucosidase")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)


#### Positive ####
#### 6.1 UniRef90_P76052 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_P76052"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p4=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#DD6A69")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_P76052\nP-aminobenzoyl-glutamate hydrolase subunit B")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)
#### 6.2 UniRef90_A0A379CI55 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_A0A379CI55"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p5=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#DD6A69")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_A0A379CI55\nmRNA interferase")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)





#### 6.3 UniRef90_A0A1C5M8F6 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Humann4.0_DMP), row.names(pheno_cov)))
covar=c( "DNA.Concentration","Sex", "log10_counts")
lm_input<-data.frame(Y_genes = Humann4.0_DMP[sample_name,"UniRef90_A0A1C5M8F6"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
if (length(unique(lm_input$X_phenotype)) > 2){
  lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
  for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
}
lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:5]))
lm_res$coefficients
lm_input$Y_Humann4.0_DMP_resid=resid(lm(Y_genes~Sex+DNA.Concentration+log10_counts,data=lm_input))
p6=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_Humann4.0_DMP_resid`))+
  geom_point(color="#DD6A69")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_Humann4.0_DMP_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of UniRef90_A0A1C5M8F6\nHTH-type transcriptional regulator immR")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,colour="black"),
        axis.text.x = element_text(face="plain",size=12,colour="black"),
        axis.title.y = element_text(face="bold",size=12,colour="black"),
        axis.title.x = element_text(face="bold",size=12,colour="black"),)

library(patchwork)
(p1+p2+p3)/(p5+p6+p4)

