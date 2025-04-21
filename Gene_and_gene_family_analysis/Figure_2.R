#### Figure 2 ####
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
all_sv_info_anno_old=all_sv_info_anno_old[which(all_sv_info_anno_old$SV_Name%in%meta_age_significant$Y),]
load("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/info.RData")
#### Figure 2B and numbers ####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred.rds")
#Final_results=Final_results[-grep("UniRef90",Final_results$Y),]
Progenome1_Bakta_genes_table <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Progenome1_Bakta_genes_table.rds")
MAG_genes_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/MAG_genes_age_significant.rds")
Final_results=Final_results[Final_results$Y%in%c(Progenome1_Bakta_genes_table$GeneName,as.character(MAG_genes_age_significant$Y)),]
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
shortbred_with_zero=Final_results
for (i in c("Y","X")){`shortbred_with_zero`[,i]=as.character(`shortbred_with_zero`[,i])}
length(unique(shortbred_with_zero$Y))#1141 genes 
length(unique(shortbred_with_zero$Y[grep("UniRef",shortbred_with_zero$Y)]))#30 MAG genes
shortbred_fdr=0.05/(1141*1)# one phenotype: age
shortbred_with_zero=shortbred_with_zero[shortbred_with_zero$X==c("Age"),]
shortbred_with_zero$replication=NA
shortbred_with_zero$replication[shortbred_with_zero$p<shortbred_fdr&shortbred_with_zero$Beta>0&shortbred_with_zero$X==c("Age")]=c("significant-positive")
shortbred_with_zero$replication[shortbred_with_zero$p<shortbred_fdr&shortbred_with_zero$Beta<0&shortbred_with_zero$X==c("Age")]=c("significant-negative")
shortbred_with_zero$replication[shortbred_with_zero$p>shortbred_fdr]=c("non-significant")
shortbred_with_zero$replication=as.factor(shortbred_with_zero$replication)
shortbred_with_zero$replication=factor(shortbred_with_zero$replication,levels = c("significant-positive","significant-negative","non-significant"))
shortbred_with_zero=shortbred_with_zero[shortbred_with_zero$Y%in%c(Progenome1_Bakta_genes_table$GeneName),]
library(ggplot2) 
library(ggprism) 
library(ggrepel) 
ggplot(shortbred_with_zero, aes(x =Beta_phenotype, y= -log10(p_phenotype), color=replication)) +
  geom_point(alpha=1, size=2.5) +
  scale_color_manual(values=c('#DF4428','#332847','gray')) +
  xlim(c(min(shortbred_with_zero$Beta), max(shortbred_with_zero$Beta))) +  ##调整x轴的取值范围，可以根据max(abs(shortbred_with_zero$Beta))，获得差异基因最大值是多少，再进行取舍
  geom_vline(xintercept=c(0),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(0.05/(1141*1)), lty=4,col="black",lwd=0.8) + 
  labs(x="Beta effect", y="-log10(P)") +
  scale_x_continuous(limits = c(-0.8, 0.8),n.breaks = 8) +
  ggtitle("Replication in LLD cohort") + 
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)
shortbred_with_zero=merge(shortbred_with_zero,Progenome1_Bakta_genes_table[,c("GeneName","Product.mixed")],by.x = "Y",by.y = "GeneName",all.x = T)
age_genes_negative=unique(shortbred_with_zero$Y[shortbred_with_zero$replication==c("significant-negative")])
length(age_genes_negative)
age_genes_positive=unique(shortbred_with_zero$Y[shortbred_with_zero$replication==c("significant-positive")])
length(age_genes_positive)
sum(length(age_genes_negative),length(age_genes_positive))
sum(length(age_genes_negative),length(age_genes_positive))/1111
#### Figure 2A: barplot compare how many genes still significant in Shortbred ####
plot_2A=shortbred_with_zero
Gene_SV_pair=NULL
for (i in unique(plot_2A$Y)){
  #i=Progenome1_Bakta_genes[1]
  #i=c("OOIAGD_13255")
  print(i)
  if (length(grep(i,all_sv_info_anno_old$Gene))>0){
    SV_tmp=all_sv_info_anno_old$SV_Name[grep(i,all_sv_info_anno_old$Gene)]
    tmp=data.frame(Gene=i,SV=SV_tmp)
  }else{
    SV_tmp=all_sv_info_anno_old$SV_Name[grep(i,all_sv_info_anno_old$Gene.bakta)]
    tmp=data.frame(Gene=i,SV=SV_tmp)
  }
  Gene_SV_pair=rbind(Gene_SV_pair,tmp)
}
plot_2A=merge(plot_2A,Gene_SV_pair,by.x = "Y",by.y = "Gene", all=T)
table1=table(plot_2A$SV[plot_2A$p>shortbred_fdr])%>%as.data.frame(.)
table2=table(plot_2A$SV[plot_2A$p<shortbred_fdr])%>%as.data.frame(.)
num_all=merge(table1,table2,by.x = "Var1",by.y = "Var1",all=T);num_all[is.na(num_all)]=0;num_all$sum=num_all$Freq.x+num_all$Freq.y
colnames(num_all)[2:3]=c("non_significant","significant")
num_all$species=sapply(strsplit(as.character(num_all$Var1),"\\:"),"[",1)
num_all$species_short=info$Short_name[match(num_all$species,info$organism)]
num_all$Var1_short=sapply(strsplit(as.character(num_all$Var1),"\\:"),"[",2)
num_all$Var1_short=paste(num_all$species_short,":",num_all$Var1_short)
df_long <- pivot_longer(num_all, cols = c(non_significant, significant), names_to = "variable", values_to = "value")
df_long$variable=as.factor(df_long$variable)
df_long$Var1_short=factor(df_long$Var1_short,levels=as.character(num_all$Var1_short[order(num_all$significant,decreasing = T)]))
ggplot(df_long, aes(x = Var1_short, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("non_significant" = "#397460", "significant" = "#E2BF9A")) +
  labs(x = "Category", y = "Sum of Values", fill = "Variable") +
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=10),
        axis.text.x = element_text(face="plain",size=10,angle=90,hjust = 1))
sum(shortbred_with_zero$replication==c("significant-positive")|shortbred_with_zero$replication==c("significant-negative"))/1111#757
#### Table S10 ####
write.table(df_long, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#### Table S9 ####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred.rds")
#Final_results=Final_results[-grep("UniRef90",Final_results$Y),]
Progenome1_Bakta_genes_table <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Progenome1_Bakta_genes_table.rds")
MAG_genes_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/MAG_genes_age_significant.rds")
Final_results=Final_results[Final_results$Y%in%c(Progenome1_Bakta_genes_table$GeneName,as.character(MAG_genes_age_significant$Y)),]
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
shortbred_with_zero=Final_results
for (i in c("Y","X")){`shortbred_with_zero`[,i]=as.character(`shortbred_with_zero`[,i])}
length(unique(shortbred_with_zero$Y))#1141
shortbred_fdr=0.05/(1141*1)
shortbred_with_zero$replication=NA
shortbred_with_zero$replication[shortbred_with_zero$p_phenotype<shortbred_fdr&shortbred_with_zero$Beta_phenotype>0]=c("significant-positive")
shortbred_with_zero$replication[shortbred_with_zero$p_phenotype<shortbred_fdr&shortbred_with_zero$Beta_phenotype<0]=c("significant-negative")
shortbred_with_zero$replication[shortbred_with_zero$p_phenotype>shortbred_fdr]=c("non-significant")
shortbred_with_zero$replication=as.factor(shortbred_with_zero$replication)
shortbred_with_zero$replication=factor(shortbred_with_zero$replication,levels = c("significant-positive","significant-negative","non-significant"))
shortbred_with_zero=merge(shortbred_with_zero,Progenome1_Bakta_genes_table[,c("GeneName","Product.mixed")],by.x = "Y",by.y = "GeneName",all.x = T)
shortbred_with_zero=shortbred_with_zero[which(shortbred_with_zero$X==c("Age")),]
shortbred_with_zero=merge(shortbred_with_zero,MAG_genes_age_significant[,c("Y","product")],by.x = "Y",by.y = "Y",all.x = T)
shortbred_with_zero$Product.mixed[is.na(shortbred_with_zero$Product.mixed)]=shortbred_with_zero$product[is.na(shortbred_with_zero$Product.mixed)]
shortbred_with_zero$product=NULL
shortbred_with_zero$fdr.p=NULL
shortbred_with_zero$origins=NA
shortbred_with_zero$origins[grep("UniRef",shortbred_with_zero$Y)]=c("From_MAGs")
shortbred_with_zero$origins[is.na(shortbred_with_zero$origins)]=c("From_SVs")
saveRDS(shortbred_with_zero,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableS9",".rds"))
write.table(shortbred_with_zero, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### Figure 2C and D ####
Shortbred <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Shortbred.rds")
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

#### HAKNHM_10150 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Shortbred), row.names(pheno_cov)))
covar=c("DNA.Concentration", "Sex")#, "log10_counts"
lm_input<-data.frame(Y_genes = Shortbred[sample_name,"HAKNHM_10150"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
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

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:4]))
lm_res$coefficients
lm_input$Y_shortbred_resid=resid(lm(Y_genes~Sex+DNA.Concentration,data=lm_input))
p1=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_shortbred_resid`))+
  geom_point(color="#DF4428")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_shortbred_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of HAKNHM_10150\nMazF family transcriptional regulator")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=16,colour="black"),
        axis.text.x = element_text(face="plain",size=16,colour="black"),
        axis.title.y = element_text(face="bold",size=16,colour="black"),
        axis.title.x = element_text(face="bold",size=16,colour="black"),)

#### RHOM_15565 ####
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
sample_name = Reduce(intersect, list(row.names(pheno_cohort), row.names(Shortbred), row.names(pheno_cov)))
covar=c("DNA.Concentration", "Sex")#, "log10_counts"
lm_input<-data.frame(Y_genes = Shortbred[sample_name,"RHOM_15565"], X_phenotype = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar],Raw_age = pheno_cohort[sample_name, "Age"]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
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

lm_res <- summary(lm(Y_genes~.,data = lm_input[,1:4]))
lm_res$coefficients
lm_input$Y_shortbred_resid=resid(lm(Y_genes~Sex+DNA.Concentration,data=lm_input))
p2=ggplot(data = lm_input, mapping = aes(x =`Raw_age`, y =`Y_shortbred_resid`))+
  geom_point(color="#332847")+
  geom_smooth(data =lm_input, mapping = aes(x =`Raw_age`,y = `Y_shortbred_resid`), method = "lm", formula = y ~ x,color="black")+
  labs(x = "Age", y = "Residuals of RHOM_15565\nResponse transcription regulator")+
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=16,colour="black"),
        axis.text.x = element_text(face="plain",size=16,colour="black"),
        axis.title.y = element_text(face="bold",size=16,colour="black"),
        axis.title.x = element_text(face="bold",size=16,colour="black"),)
p1/p2


#### Figure 2EandF ####
Progenome1_Bakta_genes_table <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Progenome1_Bakta_genes_table.rds")
Progenome1_Bakta_genes_table_1=Progenome1_Bakta_genes_table[Progenome1_Bakta_genes_table$GeneName%in%age_genes_positive,]
test=as.data.frame(table(Progenome1_Bakta_genes_table_1$KOmodule_textinfo))
Progenome1_Bakta_genes_table_1$KOmodule_textinfo[Progenome1_Bakta_genes_table_1$KOmodule_textinfo%in%as.character(test$Var1[test$Freq<10])]=c("Others")
library(ggpie)
p1=ggdonut(data = Progenome1_Bakta_genes_table_1, group_key = "KOmodule_textinfo", count_type = "full",
           label_info = "all",label_split = NULL,
           label_size = 4, label_pos = "in")
Progenome1_Bakta_genes_table_2=Progenome1_Bakta_genes_table[Progenome1_Bakta_genes_table$GeneName%in%age_genes_negative,]
test=as.data.frame(table(Progenome1_Bakta_genes_table_2$KOmodule_textinfo))
Progenome1_Bakta_genes_table_2$KOmodule_textinfo[Progenome1_Bakta_genes_table_2$KOmodule_textinfo%in%as.character(test$Var1[test$Freq<10])]=c("Others")
p2=ggdonut(data = Progenome1_Bakta_genes_table_2, group_key = "KOmodule_textinfo", count_type = "full",
           label_info = "all",label_split = NULL,
           label_size = 4, label_pos = "in")
library(patchwork)
p1+p2


#### 16 that remained significantly associated with host age ####
TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")
tmp=TableS9$Y[which(!TableS9$replication==c("non-significant")&TableS9$X==c("Age"))]
tmp%>%length(.)#772
tmp[grep("UniRef",tmp)]%>%length(.)#15
#### Table S15 ####
mag_significant_age_genes <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/mag_significant_age_genes.rds")
Progenome1_Bakta_genes_table <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Progenome1_Bakta_genes_table.rds")
write.table(Progenome1_Bakta_genes_table[, c("GeneName", "Uniref90")], "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)







