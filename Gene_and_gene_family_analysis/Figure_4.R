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
#### 2.number ####
TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")
TableS9_sig=TableS9[which(!TableS9$replication==c("non-significant")),]
#### 3.Table S17:phenotype info ####
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS16=Pheno_info[Pheno_info$UnifiedName%in%TableS18$X,]
write.table(TableS16, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
TableS9$X[which(!TableS9$X%in%Pheno_info$UnifiedName)]%>%unique(.)
#### 4.Table S18: Association results between 833 age-associated bacterial genes with host health realted metrics. ####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred.rds")
#Final_results=Final_results[-grep("UniRef90",Final_results$Y),]
Progenome1_Bakta_genes_table <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Progenome1_Bakta_genes_table.rds")
MAG_genes_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/MAG_genes_age_significant.rds")
Final_results=Final_results[Final_results$Y%in%c(Progenome1_Bakta_genes_table$GeneName,as.character(MAG_genes_age_significant$Y)),]
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
shortbred_with_zero=Final_results
for (i in c("Y","X")){`shortbred_with_zero`[,i]=as.character(`shortbred_with_zero`[,i])}
length(unique(shortbred_with_zero$Y))
length(unique(shortbred_with_zero$X))
shortbred_with_zero=shortbred_with_zero[which(shortbred_with_zero$Y%in%TableS9_sig$Y),]
shortbred_with_zero=shortbred_with_zero[which(!shortbred_with_zero$X==c("Age")),]
length(unique(shortbred_with_zero$Y))
length(unique(shortbred_with_zero$X))
shortbred_fdr=0.05/(833*111)#5.407568e-07
shortbred_with_zero$Significance=NA
shortbred_with_zero$Significance[which(shortbred_with_zero$p_phenotype<shortbred_fdr&shortbred_with_zero$Beta_phenotype>0)]=c("significant-positive")
shortbred_with_zero$Significance[which(shortbred_with_zero$p_phenotype<shortbred_fdr&shortbred_with_zero$Beta_phenotype<0)]=c("significant-negative")
shortbred_with_zero$Significance[which(shortbred_with_zero$p_phenotype>shortbred_fdr)]=c("non-significant")
shortbred_with_zero$Significance=as.factor(shortbred_with_zero$Significance)
shortbred_with_zero$Significance=factor(shortbred_with_zero$Significance,levels = c("significant-positive","significant-negative","non-significant"))
shortbred_with_zero=merge(shortbred_with_zero,Progenome1_Bakta_genes_table[,c("GeneName","Product.mixed")],by.x = "Y",by.y = "GeneName",all.x = T)
shortbred_with_zero=merge(shortbred_with_zero,MAG_genes_age_significant[,c("Y","product")],by.x = "Y",by.y = "Y",all.x = T)
shortbred_with_zero$Product.mixed[is.na(shortbred_with_zero$Product.mixed)]=shortbred_with_zero$product[is.na(shortbred_with_zero$Product.mixed)]
shortbred_with_zero$product=NULL
shortbred_with_zero$fdr.p=NULL
shortbred_with_zero$origins=NA
shortbred_with_zero$origins[grep("UniRef",shortbred_with_zero$Y)]=c("From_MAGs")
shortbred_with_zero$origins[is.na(shortbred_with_zero$origins)]=c("From_SVs")
length(unique(shortbred_with_zero$Y))
length(unique(shortbred_with_zero$X))
saveRDS(shortbred_with_zero,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableS18",".rds"))
write.table(shortbred_with_zero, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#### 5.Numbers ####
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS18_sig=TableS18[which(!TableS18$Significance==c("non-significant")),]
length(unique(TableS18_sig$Y))
length(unique(TableS18_sig$X))
length(unique(TableS18_sig$Y))/833
age_and_other_phenos_genes=unique(TableS18_sig$Y)
grep("UniRef",age_and_other_phenos_genes)%>%length(.) 
saveRDS(age_and_other_phenos_genes,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","age_and_other_phenos_genes",".rds"))
for (i in c("Y","X")){`TableS18_sig`[,i]=as.character(`TableS18_sig`[,i])}
TableS18_sig$X[TableS18_sig$X==c("FI41_FU")]=c("FI41")

#### 6.Figure S6 ####
# barplot for phenotype association number
tmp <- as.data.frame(table(TableS18_sig[,c("X"),drop=F]));
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
#### 7.1 Figure 4A ####
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS18_sig=TableS18[which(!TableS18$Significance==c("non-significant")),]
length(unique(TableS18_sig$Y))
length(unique(TableS18_sig$X))
TableS18_sig$X[TableS18_sig$X==c("FI41_FU")]=c("FI41")

TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")

tmp1=TableS9[which(!TableS9$replication==c("non-significant")),]
tmp2=tmp1$Y[tmp1$Beta_phenotype<0]
tmp3=tmp1$Y[tmp1$Beta_phenotype>0]
length(tmp2)
length(tmp3)

tmp4=TableS18_sig[TableS18_sig$X==c("FI41"),]
tmp5=tmp4$Y[tmp4$Beta_phenotype<0]
tmp6=tmp4$Y[tmp4$Beta_phenotype>0]
length(tmp5)
length(tmp6)

#Group1
#Group1a
intersect(tmp3,tmp5)%>%unique(.)%>%length(.)
one_a=intersect(tmp3,tmp5)%>%unique(.)
one_a[grep("UniRef",one_a)]%>%length(.)
saveRDS(one_a, paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","one_a_shortbred",".rds"))

#Group1b
intersect(tmp2,tmp6)%>%unique(.)%>%length(.)
one_b=intersect(tmp2,tmp6)%>%unique(.)
one_b[grep("UniRef",one_b)]%>%length(.)

#Group2
#Group2a
intersect(tmp3,tmp6)%>%unique(.)%>%length(.)
two_a=intersect(tmp3,tmp6)%>%unique(.)

#Group2b
intersect(tmp2,tmp5)%>%unique(.)%>%length(.)
two_b=intersect(tmp2,tmp5)%>%unique(.)
two_b[grep("UniRef",two_b)]%>%length(.)

#### 7.2 Table S21. The genes classified as the geroprotective pattern and the geropathogenic pattern. ####
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")

table1=TableS18[TableS18$Y%in%c(one_a,one_b),c(1:5,16)]
table1$type=c("geroprotective")
table2=TableS18[TableS18$Y%in%c(two_a,two_b),c(1:5,16)]
table2$type=c("geropathogenic")

table3=TableS9[TableS9$Y%in%c(one_a,one_b),c(1:5,16)]
table3$type=c("geroprotective")
table4=TableS9[TableS9$Y%in%c(two_a,two_b),c(1:5,16)]
table4$type=c("geropathogenic")

all_table=rbind(table1,table2,table3,table4)
all_table=all_table[all_table$X%in%c("Age","FI41_FU"),]

saveRDS(all_table,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","TableS_gene_type",".rds"))
write.table(all_table, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)

#### 8.1 Figure 4B : LLD replication ####
#DMP
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS18_sig=TableS18[which(!TableS18$Significance==c("non-significant")),]
length(unique(TableS18_sig$Y))
length(unique(TableS18_sig$X))
for (i in c("Y","X")){`TableS18_sig`[,i]=as.character(`TableS18_sig`[,i])}
TableS18_sig$X[grep("FI41_",TableS18_sig$X)]=c("FI41")
TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")
TableS9_sig=TableS9[which(!TableS9$replication==c("non-significant")),]
colnames(TableS9_sig)[15]=c("Significance")
TableS18_sig=rbind(TableS18_sig,TableS9_sig)
#LLD
Final_results_shortbred_LLD <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred_LLD.rds")
Final_results_shortbred_LLD$Beta_IDage=NULL;Final_results_shortbred_LLD$SE_IDage=NULL;Final_results_shortbred_LLD$p_IDage=NULL
for (i in c("Y","X")){`Final_results_shortbred_LLD`[,i]=as.character(`Final_results_shortbred_LLD`[,i])}
Final_results_shortbred_LLD$X[grep("FI41_",Final_results_shortbred_LLD$X)]=c("FI41")
replication_result=merge(TableS18_sig, Final_results_shortbred_LLD, by = c("X", "Y"), all.x = TRUE)
colnames(replication_result)=gsub("\\.x","_DMP",colnames(replication_result))
colnames(replication_result)=gsub("\\.y","_LLD",colnames(replication_result))
replication_result$replicated_in_LLD=NA
replication_result$replicated_in_LLD[replication_result$p_phenotype_LLD<0.05&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0]=c("Yes")
replication_result$replicated_in_LLD[replication_result$p_phenotype_LLD>0.05]=c("No")
replication_result$replicated_in_LLD[replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD<0]=c("No")
replication_result=replication_result[which(!is.na(replication_result$Y)),]
sum(replication_result$X==c("FI41")&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0)
#(228/281)
sum(replication_result$X==c("FI41")&replication_result$Beta_phenotype_DMP*replication_result$Beta_phenotype_LLD>0&replication_result$p_phenotype_LLD<0.05)
#(96/281)
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
  scale_color_manual(values = c("Yes" = "#EA8A33", "No" = "#73BAD3"))+
  #geom_text(x = 0.1, y = -0.4, label = paste0("R=",round(cor_test$estimate,2),", P = "," 0.001"), color = "black", size = 5)+   # Add text annotation for x = 0
  theme(axis.text.y = element_text(face="plain",size=16,colour = "black"),
        axis.text.x = element_text(face="plain",size=16,colour = "black"),
        axis.title.y = element_text(face="bold",size=16),
        axis.title.x = element_text(face="bold",size=16),)
p+xlim(-0.5,0.3)+ylim(-0.5,0.3)
#### 8.2 Table S19. LLD replication results for the genes both associated with age and health-realted metrices discovered in DMP cohort. ####
saveRDS(replication_result,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","replication_result_LLD",".rds"))
write.table(replication_result, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#### 7. Figure 4C ####
TableS18 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS18.rds")
TableS18_sig=TableS18[which(!TableS18$Significance==c("non-significant")),]
length(unique(TableS18_sig$Y))#430
length(unique(TableS18_sig$X))#57 phenotypes
for (i in c("Y","X")){`TableS18_sig`[,i]=as.character(`TableS18_sig`[,i])}
TableS18_sig$X[grep("FI41_",TableS18_sig$X)]=c("FI41")
TableS9 <- readRDS("~/Documents/SV_MWAS/NM_reviosion/TableS9.rds")
Figure_4C_plot=rbind(TableS18_sig[which(TableS18_sig$Y%in%TableS18_sig$Y[TableS18_sig$X==c("FI41")]),c("X","Y","Beta_phenotype")],
                     TableS9[which(TableS9$Y%in%TableS18_sig$Y[TableS18_sig$X==c("FI41")]),c("X","Y","Beta_phenotype")])
  
#### 7.1 R loading ####
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
#### 7.2 Plot matrix ####
plot_matrix <- dcast(Figure_4C_plot, X ~ Y, value.var = "Beta_phenotype");row.names(plot_matrix)=plot_matrix$X;plot_matrix=plot_matrix[,-1]
plot_matrix = plot_matrix[rowSums(!is.na(plot_matrix))>1,]
plot_matrix = plot_matrix[row.names(plot_matrix)%in%c("Age","FI41",
                                                      "BMI","Glucose","HbA1c","Hemoglobin",
                                                      "HeartRate","SBP","DBP",
                                                      "TG","HDL","LDL"),]
plot_matrix=as.data.frame(t(plot_matrix))
plot_matrix=plot_matrix[,c("Age","FI41","BMI","Glucose","HbA1c","Hemoglobin","HeartRate","SBP","DBP","TG","HDL")]
plot_matrix=as.data.frame(t(plot_matrix))
plot_matrix[is.na(plot_matrix)]=0
annotation=data.frame(gene=colnames(plot_matrix),num=1:length(colnames(plot_matrix)))
write.table(annotation, paste0("/Users/helloduck/Desktop/","table1",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
colnames(plot_matrix)=1:length(colnames(plot_matrix))
data <- as.matrix(plot_matrix)
cir1<-t(data)
mycol2=colorRamp2(c(min(data), 0, max(data)),c("#73BAD3","white","#E53F24")) 
#### 7.3 gene_groups ####
gene_groups <- data.frame(annotation,group = NA)
gene_groups$group[gene_groups$gene%in%c(one_a,one_b)]=c("Geroprotective genes")
gene_groups$group[gene_groups$gene%in%c(two_a,two_b)]=c("Geropathogenic genes")
gene_groups$group=factor(gene_groups$group,levels=c("Geroprotective genes","Geropathogenic genes"))
gene_groups$gene=NULL
row.names(gene_groups)=gene_groups$num;gene_groups$num=NULL
gene_groups <- as.matrix(gene_groups)
#### 7.4 plot for value ####
circos.clear()
circos.par(start.degree = 90,gap.after=c(2,30)) # circos.par()调整圆环首尾间的距离，数值越大，距离越宽;让分裂的一个口大一点，可以添加行信息
circos.heatmap(t(data),col=mycol2,split = gene_groups,
               cell.border = "white",
               #dend.side="inside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               rownames.side="outside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.4, # 轨道的高度，数值越大圆环越粗
               rownames.col="black",
               bg.border="black", # 背景边缘颜色
               # # 用行注释分裂热图
               #show.sector.labels = T,
               rownames.cex=0.4, # 字体大小
               rownames.font=0.8, # 字体粗细
               cluster=F, # cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               #dend.track.height=0.18,#调整行聚类树的高度
               #dend.callback=function(dend,m,si) {
               #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
               #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
               #color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               #}
)
#图例与列名设置
lg=Legend(title="Legend",col_fun=mycol2,direction = c("horizontal"))
grid.draw(lg)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 2) { # 2 is the last sector
    cn = colnames(t(data))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                (1:n)*1.2 - 0.7, cn, 
                cex = 0.5, adj = c(0, -0.5), facing = "outside")
  }
}, bg.border = NA)
#### 7.5 the box border for LLD replication ####
plot_matrix <- dcast(Figure_4C_plot, X ~ Y, value.var = "Beta_phenotype");row.names(plot_matrix)=plot_matrix$X;plot_matrix=plot_matrix[,-1]
plot_matrix = plot_matrix[rowSums(!is.na(plot_matrix))>1,]
plot_matrix = plot_matrix[row.names(plot_matrix)%in%c("Age","FI41",
                                                      "BMI","Glucose","HbA1c","Hemoglobin",
                                                      "HeartRate","SBP","DBP",
                                                      "TG","HDL","LDL"),]
plot_matrix=as.data.frame(t(plot_matrix))
plot_matrix=plot_matrix[,c("Age","FI41","BMI","Glucose","HbA1c","Hemoglobin","HeartRate","SBP","DBP","TG","HDL")]
plot_matrix=as.data.frame(plot_matrix)
plot_matrix[is.na(plot_matrix)]=0
rep_test <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_LLD.rds")
rep_test=rep_test[which(rep_test$replicated_in_LLD==c("Yes")),]
for (i in colnames(plot_matrix)){
  print(i)
  #i= colnames(plot_matrix)[1]
  for (j in 1:nrow(plot_matrix)){
    print(j)
    #j=1
    if (row.names(plot_matrix)[j]%in%rep_test$Y[which(rep_test$X==i)]){plot_matrix[j,i]=1}else{plot_matrix[j,i]=0}
  }
}
plot_matrix=as.data.frame(t(plot_matrix))
colnames(plot_matrix)=1:length(colnames(plot_matrix))
data <- as.matrix(plot_matrix)

#### 7.6 plot for replication ####
circos.clear()
circos.par(start.degree = 90,gap.after=c(2,30)) # circos.par()调整圆环首尾间的距离，数值越大，距离越宽;让分裂的一个口大一点，可以添加行信息
circos.heatmap(t(data),col=mycol2,split = gene_groups,
               cell.border = "white",
               #dend.side="inside",# dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               rownames.side="outside",# rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               track.height = 0.4, # 轨道的高度，数值越大圆环越粗
               rownames.col="black",
               bg.border="black", # 背景边缘颜色
               # # 用行注释分裂热图
               #show.sector.labels = T,
               rownames.cex=0.4, # 字体大小
               rownames.font=0.8, # 字体粗细
               cluster=F, # cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               #dend.track.height=0.18,#调整行聚类树的高度
               #dend.callback=function(dend,m,si) {
               #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，
               #或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
               #color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
               #}
)
#图例与列名设置
lg=Legend(title="Legend",col_fun=mycol2,direction = c("horizontal"))
grid.draw(lg)
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 2) { # 2 is the last sector
    cn = colnames(t(data))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                (1:n)*1.2 - 0.7, cn, 
                cex = 0.5, adj = c(0, -0.5), facing = "outside")
  }
}, bg.border = NA)

#### 8.0 Forest plots (Figure 4 and E) ####
rep_LLD <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_LLD.rds")
rep_LLD=rep_LLD[which(rep_LLD$replicated_in_LLD==c("Yes")),]
rep_LLD$group <- NA
rep_LLD$group[rep_LLD$Y%in%c(one_a,one_b)]=c("Geroprotective genes")
rep_LLD$group[rep_LLD$Y%in%c(two_a,two_b)]=c("Geropathogenic genes")
rep_LLD$group=as.factor(rep_LLD$group)
rep_LLD$conf.low=rep_LLD$Beta_phenotype_DMP-1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$conf.high=rep_LLD$Beta_phenotype_DMP+1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$model=paste0(rep_LLD$Y,"~",rep_LLD$X)
rep_LLD$Y=as.factor(rep_LLD$Y)
#geroprotective  
plot=rep_LLD[rep_LLD$X%in%c("Age","FI41")&rep_LLD$Y%in%c("HMPREF9436_03248","HMPREF9436_03250","HMPREF9436_03251",
                                                    "JPJG01000078_gene905","JPJG01000088_gene696","JPJG01000116_gene2019"),] 
plot$model <- as.factor(as.character(plot$model))
plot$Y <- as.factor(as.character(plot$Y))
plot$row_id <- as.numeric(plot$model)
library(wesanderson)
ggplot(plot,aes(x=fct_rev(model),y=Beta_phenotype_DMP,ymin=conf.low,ymax=conf.high)) +
  geom_rect(aes(xmin = row_id - 0.5, xmax = row_id + 0.5, ymin = -Inf, ymax = Inf), 
            fill = rep(c("gray", "white"), length.out = nrow(plot)), 
            alpha = 1, inherit.aes = FALSE)+
  geom_pointrange(aes(color = Y), size = .5) +
  scale_color_manual(values = wes_palette(6, name = "AsteroidCity2", type = "discrete"))+
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
rep_LLD <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/replication_result_LLD.rds")
rep_LLD$group <- NA
rep_LLD$group[rep_LLD$Y%in%c(one_a,one_b)]=c("Geroprotective genes")
rep_LLD$group[rep_LLD$Y%in%c(two_a,two_b)]=c("Geropathogenic genes")
rep_LLD$group=as.factor(rep_LLD$group)
rep_LLD$conf.low=rep_LLD$Beta_phenotype_DMP-1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$conf.high=rep_LLD$Beta_phenotype_DMP+1.96*rep_LLD$SE_phenotype_DMP
rep_LLD$model=paste0(rep_LLD$Y,"~",rep_LLD$X)
rep_LLD$Y=as.factor(rep_LLD$Y)
plot=rep_LLD[rep_LLD$X%in%c("Age","FI41")&rep_LLD$Y%in%c("UniRef90_A0A1Q6QCA2","UniRef90_A0A1Q6QC78","UniRef90_A0A1Q6QC66","UniRef90_A0A1Q6QC59",
                                                    "RHOM_15555","HMPREF9436_01363"),]
plot$model <- as.factor(as.character(plot$model))
plot$Y <- as.factor(as.character(plot$Y))
plot$row_id <- as.numeric(plot$model)
library(wesanderson)
ggplot(plot,aes(x=fct_rev(model),y=Beta_phenotype_DMP,ymin=conf.low,ymax=conf.high,fill=Y)) +
  geom_rect(aes(xmin = row_id - 0.5, xmax = row_id + 0.5, ymin = -Inf, ymax = Inf), 
            fill = rep(c("gray", "white"), length.out = nrow(plot)), 
            alpha = 1, inherit.aes = FALSE)+
  geom_pointrange(aes(color = Y),  size = .5) +
  scale_color_manual(values = wes_palette(7, name = "Rushmore1", type = "continuous"))+
  geom_hline(yintercept = 0, color = "steelblue") +  
  xlab(" ") +
  ylab("Coefficient (95% Confidence Interval)") +
  labs(title ="Linear Regression Models Estimating the Effects of Vehicle Weight \n on Fuel Efficiency") +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12))+
  coord_flip()

#### Figure S9: biotin gene levels in DMP and LLD cohort ####
Shortbred_DMP <- readRDS("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Shortbred.rds")
Shortbred_LLD <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Shortbred_LLD.rds")

tmp1=data.frame(Shortbred_DMP[,c("UniRef90_A0A1Q6QC59","UniRef90_A0A1Q6QCA2","UniRef90_A0A1Q6QC78","UniRef90_A0A1Q6QC66")],cohort=c("DMP"))
tmp2=data.frame(Shortbred_LLD[,c("UniRef90_A0A1Q6QC59","UniRef90_A0A1Q6QCA2","UniRef90_A0A1Q6QC78","UniRef90_A0A1Q6QC66")],cohort=c("LLD1"))
all_biotin=rbind(tmp1,tmp2);all_biotin$cohort=as.factor(all_biotin$cohort)
all_biotin$UniRef90_A0A1Q6QC59[all_biotin$UniRef90_A0A1Q6QC59==0]=min(all_biotin$UniRef90_A0A1Q6QC59[all_biotin$UniRef90_A0A1Q6QC59>0],na.rm = T)/2
all_biotin$UniRef90_A0A1Q6QC59=log2(all_biotin$UniRef90_A0A1Q6QC59)

all_biotin$UniRef90_A0A1Q6QCA2[all_biotin$UniRef90_A0A1Q6QCA2==0]=min(all_biotin$UniRef90_A0A1Q6QCA2[all_biotin$UniRef90_A0A1Q6QCA2>0],na.rm = T)/2
all_biotin$UniRef90_A0A1Q6QCA2=log2(all_biotin$UniRef90_A0A1Q6QCA2)

all_biotin$UniRef90_A0A1Q6QC78[all_biotin$UniRef90_A0A1Q6QC78==0]=min(all_biotin$UniRef90_A0A1Q6QC78[all_biotin$UniRef90_A0A1Q6QC78>0],na.rm = T)/2
all_biotin$UniRef90_A0A1Q6QC78=log2(all_biotin$UniRef90_A0A1Q6QC78)

all_biotin$UniRef90_A0A1Q6QC66[all_biotin$UniRef90_A0A1Q6QC66==0]=min(all_biotin$UniRef90_A0A1Q6QC66[all_biotin$UniRef90_A0A1Q6QC66>0],na.rm = T)/2
all_biotin$UniRef90_A0A1Q6QC66=log2(all_biotin$UniRef90_A0A1Q6QC66)

library(PupillometryR)
p1=ggplot(all_biotin, aes(x = cohort, y = UniRef90_A0A1Q6QC59)) + 
  #geom_boxplot(aes(fill = cohort),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  geom_point(aes(color = cohort),position = position_jitter(width = 0.13), size =1,alpha = 0.5)+
  geom_flat_violin(aes(fill = cohort,color = cohort),trim = T,alpha=0.4)+
  labs(x = "", y = colnames(all_biotin)[1]) +
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))+
  scale_fill_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))
p2=ggplot(all_biotin, aes(x = cohort, y = UniRef90_A0A1Q6QCA2)) + 
  #geom_boxplot(aes(fill = cohort),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  geom_point(aes(color = cohort),position = position_jitter(width = 0.13), size =1,alpha = 0.5)+
  geom_flat_violin(aes(fill = cohort,color = cohort),trim = T,alpha=0.4)+
  labs(x = "", y = colnames(all_biotin)[2]) +
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))+
  scale_fill_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))
p3=ggplot(all_biotin, aes(x = cohort, y = UniRef90_A0A1Q6QC78)) + 
  #geom_boxplot(aes(fill = cohort),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  geom_point(aes(color = cohort),position = position_jitter(width = 0.13), size =1,alpha = 0.5)+
  geom_flat_violin(aes(fill = cohort,color = cohort),trim = T,alpha=0.4)+
  labs(x = "", y = colnames(all_biotin)[3]) +
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))+
  scale_fill_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))
p4=ggplot(all_biotin, aes(x = cohort, y = UniRef90_A0A1Q6QC66)) + 
  #geom_boxplot(aes(fill = cohort),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  geom_point(aes(color = cohort),position = position_jitter(width = 0.13), size =1,alpha = 0.5)+
  geom_flat_violin(aes(fill = cohort,color = cohort),trim = T,alpha=0.4)+
  labs(x = "", y = colnames(all_biotin)[4]) +
  theme_bw()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))+
  scale_fill_manual(values = c("DMP" = "#ff595e", "LLD1" = "#1982c4"))
library(patchwork);p1+p2+p3+p4


p1=ggplot(all_biotin, aes(x=UniRef90_A0A1Q6QC59, fill=cohort)) + geom_histogram(alpha=.3)
p2=ggplot(all_biotin, aes(x=UniRef90_A0A1Q6QCA2, fill=cohort)) + geom_histogram(alpha=.3)
p3=ggplot(all_biotin, aes(x=UniRef90_A0A1Q6QC78, fill=cohort)) + geom_histogram(alpha=.3)
p4=ggplot(all_biotin, aes(x=UniRef90_A0A1Q6QC66, fill=cohort)) + geom_histogram(alpha=.3)
library(patchwork);p1+p2+p3+p4
wilcox.test(UniRef90_A0A1Q6QC59~cohort,all_biotin)$p.value
wilcox.test(UniRef90_A0A1Q6QCA2~cohort,all_biotin)$p.value
wilcox.test(UniRef90_A0A1Q6QC78~cohort,all_biotin)$p.value
wilcox.test(UniRef90_A0A1Q6QC66~cohort,all_biotin)$p.value





