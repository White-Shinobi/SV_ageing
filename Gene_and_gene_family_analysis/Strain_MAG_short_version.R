#source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step2/", "Part1_functions.R"))
#source(paste0("/scratch/p303998/SV_MWAS/Rdata_1216/Step1/", "Part3_functions.R"))
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step1/", "Part3_functions.R"))
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step2/", "Part1_functions.R"))

#### 0. pick out the MAGs for O.ER4 ####
# all MAGs infomation
#mag_file <- read.delim("/scratch/p303998/SV_MWAS/bins_summary.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
mag_file <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/bins_summary.csv",sep=",",row.names = NULL,header = T,check.names = F,fill = F)
# O.ER4: g__ER4;s__ER4 sp000765235; Completeness>50 & Contamination<5
mag_file = mag_file[which(mag_file$GTDBK.GenSpec==c("g__ER4;s__ER4 sp000765235")&mag_file$CheckM.Completeness>50&mag_file$CheckM.Contamination<5),]
mag_file$Bin.Folder.Path=gsub(".fa","",mag_file$Bin.File.Path)
mag_file$Gff3.File.Path=paste0(mag_file$Bin.Folder.Path,"/",sapply(strsplit(as.character(mag_file$Bin.Folder.Path), "\\/"), "[", 9),".gff3")
# Completeness>90 & Contamination<5
mag_file_90=mag_file[which(mag_file$CheckM.Completeness>90&mag_file$CheckM.Contamination<5),]
mag_file_90$Gff3.File.Path=sapply(strsplit(as.character(mag_file_90$Gff3.File.Path), "\\/"), "[", 10)

# reference MAGs infomation
# ls /scratch/p303998/SV_MWAS/Shortbred//O_ER4_NCBI_genomes/ncbi_dataset/data/OER4_NCBI_MAGs/*/*gff3 >> NCBI_OER4.tsv
# NCBI_OER4_file <- read.delim("/scratch/p303998/SV_MWAS/Shortbred/NCBI_OER4.tsv",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
NCBI_OER4_file <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/NCBI_OER4.tsv",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
NCBI_OER4_file$V2 <- sapply(strsplit(as.character(NCBI_OER4_file$V1), "\\/"), "[", 12)
NCBI_OER4_file$V3=gsub(".gff3","",NCBI_OER4_file$V2)

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
load("./Phenotype_rawdata/full_phen.RData")
full_phen=removeColsAllNa(full_phen)
full_phen=merge(full_phen,physical_score[,c("total_scor_VAL","MVPA_scor_VAL")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
full_phen=merge(full_phen,Frailty_index[,c("FI64","FI60","FI41_B","FI41_FU","FI39_B","FI39_FU")],by.x = "row.names",by.y = "row.names",all.x = T)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
Pheno_info <- read.delim("./Phenotype_info/phenInfoSumm.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
Pheno_info=Pheno_info[!Pheno_info$Class==c(""),]
Pheno_info=Pheno_info[-which(Pheno_info$UnifiedName==c("MeatFreqPerWeek")&Pheno_info$Unit==c("4 point scale")),]
load("./readcounts_rawdata/full_read.RData")
full_read$log10_counts=log10(full_read$V2)
pheno_cov=merge(full_phen[,c("DNA.Concentration","Age","Sex","BristolType")],full_read[,c("V1","log10_counts"),drop=F],by.x = "row.names",by.y = "V1",all=F)%>% `rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
MAG_all_unique <- readRDS("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/MAG_all_unique.rds")
plot_matrix <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/plot_matrix.rds");plot_matrix=t(plot_matrix)

#### 1.1 gene 0/1 matrix ####
#MAG_genes <- readRDS("/scratch/p303998/SV_MWAS/MAG_genes.rds")
# MAG_genes <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/MAG_genes.rds")
# plot_matrix <- dcast(MAG_genes, UniRef90 ~ MAGID, value.var = "Freq");row.names(plot_matrix)=plot_matrix$UniRef90;plot_matrix=plot_matrix[,-1]
#plot_matrix <- dcast(MAG_genes, UniRef50 ~ MAGID, value.var = "Freq");row.names(plot_matrix)=plot_matrix$UniRef50;plot_matrix=plot_matrix[,-1]
# plot_matrix[is.na(plot_matrix)]=0
#saveRDS(plot_matrix,paste0("/scratch/p303998/SV_MWAS/","plot_matrix",".rds"))#Uniref90
# saveRDS(plot_matrix,paste0("~/Documents/SV_MWAS/R/Raw_data/MAG/","plot_matrix",".rds"))#Uniref90
#saveRDS(plot_matrix,paste0("/scratch/p303998/SV_MWAS/","plot_matrix_uniref50",".rds"))
#Completeness>90&Contamination<5
plot_matrix_90=plot_matrix[row.names(plot_matrix)%in%c(mag_file_90$Bin.ID),]
plot_matrix_90=plot_matrix_90[,colSums(plot_matrix_90)>0.05*nrow(plot_matrix_90)&colSums(plot_matrix_90)<0.95*nrow(plot_matrix_90)] # Shell genes: 2022
plot_matrix_90=plot_matrix_90[rowSums(plot_matrix_90)>0.05*ncol(plot_matrix_90),]
plot_matrix_90=as.data.frame(plot_matrix_90)
plot_matrix_90$ID=sapply(strsplit(as.character(row.names(plot_matrix_90)), "\\."), "[", 1)
plot_matrix_90=plot_matrix_90[!duplicated(plot_matrix_90$ID),]#one duplicated MAG
plot_matrix_90$ID=NULL
row.names(plot_matrix_90)=sapply(strsplit(as.character(row.names(plot_matrix_90)), "\\."), "[", 1)
#### 1.2 jaccard_gene distance and k-maeans cluster ####
# 90
jaccard_gene_distance=as.data.frame(as.matrix(vegdist(plot_matrix_90),index="jaccard"))
# saveRDS(jaccard_gene_distance,paste0("~/Documents/SV_MWAS/R/Raw_data/MAG/","jaccard_gene_distance_90",".rds"))
jaccard_gene_distance=jaccard_gene_distance[row.names(plot_matrix_90),row.names(plot_matrix_90)]
# PCOA
pcoa_res<-cmdscale(jaccard_gene_distance, k=2, eig = T)
pcoa <- data.frame(pcoa_res$points)
# Choose best cluster number
wss <- NULL
# For 1 to 10 cluster centres
for (i in 1:10) {
  set.seed(123)
  kmeans_clusters <- 
    kmeans(pcoa, 
           centers = i, # just runs the solution for 1-10 clusters
           nstart = 20, 
           iter.max = 200, 
           algorithm = "MacQueen") 
  # Append the tot sum of squares i.e. compactness of the clusters for each run into our empty object
  wss <- rbind(wss, tibble("clusters" = i, "WSS" = kmeans_clusters$tot.withinss)) # record WSS and number of clusters
}
# Plot total within sum of squares vs. number of clusters
ggplot(wss, aes(x=clusters, y=WSS)) +
  geom_point() + # add dots to the plot
  geom_line() + # connect the dots with a line
  scale_x_continuous(breaks= pretty_breaks()) + # tick marks on whole numbers
  labs(x = "Number of Clusters", y ="Within cluster sum of squares") # add titles to the axes
pcoa_clusters <- 
  kmeans(pcoa, # our post-PCA data set
         centers = 3, # how many clusters we'd like
         nstart = 20, # how many times we repeat the process with different random initialisations
         iter.max = 200, # how many iterations to run k-means for
         algorithm = "MacQueen") # The default algorithm can struggle with close points
#### 1.3 Association between clades and age,phenotypes ####
#pcoa_with_clusters <- bind_cols(pcoa, pcoa_kmeans_cluster=pcoa_clusters$cluster)
#colnames(pcoa_with_clusters)[3]=c("pcoa_kmeans_cluster")
#pcoa_with_clusz¸ters$ID=sapply(strsplit(as.character(row.names(pcoa_with_clusters)), "\\."), "[", 1)
#pcoa_with_clusters=pcoa_with_clusters[!duplicated(pcoa_with_clusters$ID),]
#pheno_cohort=full_phen[sample_number$New_SV_ID[sample_number$Cohort==c("DMP")],]
#pcoa_with_clusters=merge(pcoa_with_clusters,pheno_cohort,by.x = "ID",by.y = "row.names",all.x=T)
#pcoa_with_clusters=pcoa_with_clusters[!is.na(pcoa_with_clusters$Age),]
#pcoa_with_clusters$pcoa_kmeans_cluster=factor(pcoa_with_clusters$pcoa_kmeans_cluster,levels = c("1","2","3"))
#pcoa_with_clusters=merge(pcoa_with_clusters,pheno_cov[,c("log10_counts"),drop=F],by.x = "ID",by.y = "row.names",all.x=T)
#saveRDS(pcoa_with_clusters,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","pcoa_with_clusters",".rds"))
pcoa_with_clusters <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/pcoa_with_clusters.rds")
# Age
# pcoa_with_clusters=pcoa_with_clusters[pcoa_with_clusters$Age>18,]
# lm_abundance_phenotype is comparing strain1 with the average of strain2&3
covar = c("DNA.Concentration", "Sex", "log10_counts")
strain_age_result=lm_abundance_phenotype_pairwise(pcoa_with_clusters[,c("Age"),drop=F],pcoa_with_clusters[,c("pcoa_kmeans_cluster"),drop=F],pcoa_with_clusters[,covar],covar)
# Other phenotypes
covar = c("Age","DNA.Concentration", "Sex", "log10_counts")
strain_otherpheno_result=lm_abundance_phenotype(pcoa_with_clusters[,setdiff(names(pcoa_with_clusters), c(covar,c("ID","X1","X2","pcoa_kmeans_cluster")))],pcoa_with_clusters[,c("pcoa_kmeans_cluster"),drop=F],pcoa_with_clusters[,covar],covar)
Stain_result_all=rbind(strain_age_result,strain_otherpheno_result)
# for (i in c("beta_1_2","se_1_2","p.value_1_2","beta_1_3","se_1_3","p.value_1_3","beta_2_3","se_2_3","p.value_2_3","N","y_uniq_N","x_uniq_N","Strain1_sample_number","Strain2_sample_number","Strain3_sample_number","Phenotype_0_sample_number","Phenotype_1_sample_number","Strain1_sample_rate","Strain2_sample_rate")){Stain_result_all[,i]=as.numeric(Stain_result_all[,i])}
for (i in c("beta_1","se_1","p.value_1","beta_2","se_2","p.value_2","beta_3","se_3","p.value_3","N","y_uniq_N","x_uniq_N","Strain1_sample_number","Strain2_sample_number","Strain3_sample_number","Phenotype_0_sample_number","Phenotype_1_sample_number","Strain1_sample_rate","Strain2_sample_rate")){Stain_result_all[,i]=as.numeric(Stain_result_all[,i])}
Stain_result_all_filter=Stain_result_all[Stain_result_all$N>80,]
Stain_result_all_filter=Stain_result_all_filter[!is.na(Stain_result_all_filter$p.value_1),]
Stain_result_all_filter=Stain_result_all_filter[-which(Stain_result_all_filter$y_uniq_N==2&Stain_result_all_filter$Phenotype_0_sample_number/Stain_result_all_filter$N>0.9),]
Stain_result_all_filter=Stain_result_all_filter[-which(Stain_result_all_filter$y_uniq_N==2&Stain_result_all_filter$Phenotype_0_sample_number/Stain_result_all_filter$N<0.1),]
Stain_result_all_filter=Stain_result_all_filter[!Stain_result_all_filter$Phenotype%in%c("FI64","FI60","FI39_B","FI39_FU","BristolType","BristolFreq","CleanReadCount","Latitude","Longitude","FattyLiverIndexT1Class","Rome3IBS.Factor"),]
Stain_result_all_filter$fdr.1=p.adjust(Stain_result_all_filter$p.value_1,method = "fdr")
Stain_result_all_filter$fdr.2=p.adjust(Stain_result_all_filter$p.value_2,method = "fdr")
Stain_result_all_filter$fdr.3=p.adjust(Stain_result_all_filter$p.value_3,method = "fdr")
saveRDS(Stain_result_all_filter,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","Stain_result_all_filter",".rds"))
saveRDS(Stain_result_all_filter,paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","Stain_result_all_filter_pairwise",".rds"))
Stain_result_all_filter_pairwise <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/Stain_result_all_filter_pairwise.rds")
Stain_result_all_filter <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/Stain_result_all_filter.rds")
write.table(Stain_result_all_filter_pairwise, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","Stain_result_all_filter_pairwise",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
write.table(Stain_result_all_filter, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","Stain_result_all_filter_50",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
#### 1.4 PCOA plot ####
varExp=(pcoa_res$eig/sum(pcoa_res$eig))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)
pcoa_with_clusters$age_group <- cut(pcoa_with_clusters$Age,
                                    breaks = c(-Inf, 18, 31, 41,51, 66, Inf),
                                    labels = c("<18", "18-30", "31-40", "41-50", "51-65","65+"),
                                    right = FALSE)
pcoa_with_clusters$pcoa_kmeans_cluster=as.factor(pcoa_with_clusters$pcoa_kmeans_cluster)
pcoa_with_clusters[,"Age"]=qtrans(pcoa_with_clusters[,"Age"])
# 1:<18
ggplot(pcoa_with_clusters,aes(X1,X2,
                              color = pcoa_with_clusters$pcoa_kmeans_cluster, 
                              #shape = as.factor(pcoa_with_clusters$pcoa_kmeans_cluster)
))+
  geom_point(size = 2,alpha = 0.6)+
  xlab(paste("PCo1=",round(xVar,digits = 2),"%",sep = ""))+
  ylab(paste("PCo2=",round(yVar,digits = 2),"%",sep = ""))+
  theme(plot.title = element_text(size=10, face="italic"),
        plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line.x =  element_line(),
        axis.line.y = element_line(),
        #legend.position = 'none',
        legend.title = element_blank(),
        legend.key = element_rect(fill = NA), 
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA),
        panel.background = element_rect(fill = NA))+ 
  scale_color_manual(values = c("1" = "#FF5575", "2" = "#FFD36A", "3" = "#6299FF"))
ggplot(pcoa_with_clusters, aes(x = pcoa_kmeans_cluster, y = Age)) + 
  geom_flat_violin(aes(fill = pcoa_kmeans_cluster),trim = T,alpha=0.4,color = "white")+
  geom_point(aes(color = pcoa_kmeans_cluster),position = position_jitter(width = 0.13), size =1.3)+
  geom_boxplot(aes(fill = pcoa_kmeans_cluster),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  labs(x = "Oscillibacter sp. ER4 subspecies", y = "Age") +
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("1" = "#FF5575", "2" = "#FFD36A", "3" = "#6299FF"))+
  scale_fill_manual(values = c("1" = "#FF5575", "2" = "#FFD36A", "3" = "#6299FF"))

#### 2. pangenome association: gene ~ age ####
# pangenome gene
#pcoa_with_clusters <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/pcoa_with_clusters.rds")
#plot_matrix <- readRDS("/scratch/p303998/SV_MWAS/plot_matrix.rds")
plot_matrix_90=plot_matrix_90[,row.names(plot_matrix_90)%in%pcoa_with_clusters$ID] #1921 DAG3 uniref genes
# phenotype
cohort = c("dmp")
pheno_cohort = full_phen[row.names(full_phen) %in% sample_number$New_SV_ID[sample_number$Cohort_2 == cohort],c("Age","BMI")]
covar=c("DNA.Concentration", "Sex", "log10_counts")
sample_name = Reduce(intersect, list(row.names(plot_matrix_90),row.names(pheno_cov),row.names(pheno_cohort)))
O_ER4_MAG_gene_age_res=glm_gene_SV(pheno_cohort[sample_name,"Age",drop=F],plot_matrix_90[sample_name,],c("DMP"),pheno_cov[sample_name,],covar) 
for (i in c("Beta","SE","p","N","y_uniq_N","x_uniq_N","y_non_zero_N","x_non_zero_N","y_non_zero_rate","x_non_zero_rate")){O_ER4_MAG_gene_age_res[,i]=as.numeric(O_ER4_MAG_gene_age_res[,i])}
#saveRDS(O_ER4_MAG_gene_age_res,paste0("/scratch/p303998/SV_MWAS/","O_ER4_MAG_gene_age_res",".rds"))
#O_ER4_MAG_gene_age_res <- readRDS("/scratch/p303998/SV_MWAS/O_ER4_MAG_gene_age_res.rds")
O_ER4_MAG_gene_age_res=merge(O_ER4_MAG_gene_age_res,MAG_all_unique,by.x = "Y",by.y = "UniRef90",all.x =T)
O_ER4_MAG_gene_age_res=O_ER4_MAG_gene_age_res[!duplicated(O_ER4_MAG_gene_age_res$Y),]
saveRDS(O_ER4_MAG_gene_age_res,paste0("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/","O_ER4_MAG_gene_age_res",".rds"))
O_ER4_MAG_gene_age_res$p.bonforroni=p.adjust(O_ER4_MAG_gene_age_res$p,method="bonferroni")
sum(is.na(O_ER4_MAG_gene_age_res$p))
res_significant=O_ER4_MAG_gene_age_res[O_ER4_MAG_gene_age_res$p<(0.05/(length(unique(O_ER4_MAG_gene_age_res$Y))*length(unique(O_ER4_MAG_gene_age_res$X)))),]
#saveRDS(res_significant,paste0("/scratch/p303998/SV_MWAS/","MAG_genes_age_significant",".rds"))
saveRDS(res_significant,paste0("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/","MAG_genes_age_significant",".rds"))
saveRDS(as.character(O_ER4_MAG_gene_age_res$Y[O_ER4_MAG_gene_age_res$p.bonforroni<0.05]),paste0("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/MAG/","mag_significant_age_genes",".rds"))


#### 2.1 heatmap ####
pcoa_with_clusters <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/pcoa_with_clusters.rds")
plot_matrix <- readRDS("/scratch/p303998/SV_MWAS/plot_matrix.rds")
plot_matrix_90=plot_matrix[,colnames(plot_matrix)%in%c(mag_file_90$Bin.ID)] #2023 DAG3 MAGs
plot_matrix_90=plot_matrix_90[rowSums(plot_matrix_90)>0.05*ncol(plot_matrix_90)&rowSums(plot_matrix_90)<0.95*ncol(plot_matrix_90),] # Shell genes: 2022
plot_matrix_90=plot_matrix_90[,colSums(plot_matrix_90)>0.05*nrow(plot_matrix_90)]
plot_matrix_90=as.data.frame(plot_matrix_90)
pangenome_mag <- readRDS("/scratch/p303998/SV_MWAS/pangenome_mag.rds")
plot_matrix_90=plot_matrix_90[,pangenome_mag]#one duplicated MAG
colnames(plot_matrix_90)=sapply(strsplit(as.character(colnames(plot_matrix_90)), "\\."), "[", 1)
plot_matrix_90=plot_matrix_90[,colnames(plot_matrix_90)%in%pcoa_with_clusters$ID] #1921 DAG3 MAGs
test=MAG_all_unique[MAG_all_unique$UniRef90%in%row.names(plot_matrix_90),];test=test[!duplicated(test$UniRef90),];write.table(test, paste0("/scratch/p303998/SV_MWAS/Rdata_0828/","genes",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
nrow(plot_matrix_90)-sum(rowSums(plot_matrix_90)==0)
#plot_matrix_90=plot_matrix_90[1:10,]
# zoom-in part
MAG_genes_age_significant <- readRDS("/scratch/p303998/SV_MWAS/MAG_genes_age_significant.rds")
biotin_genes=c(MAG_genes_age_significant$Y[grep("iotin",MAG_genes_age_significant$product)],MAG_genes_age_significant$Y[grep("adenosylmethionine--8-amino-7-oxononanoate transaminase",MAG_genes_age_significant$product)])%>%as.character(.)
phenylalanine_genes=c(MAG_genes_age_significant$Y[grep("Phenylalanine",MAG_genes_age_significant$product)],MAG_genes_age_significant$Y[grep("hydroxybutyryl-CoA",MAG_genes_age_significant$product)])%>%as.character(.)
plot_matrix_90=plot_matrix_90[c(biotin_genes,phenylalanine_genes),]

#plot_matrix_90=plot_matrix_90[,colnames(plot_matrix_90)%in%MAG_genes_age_significant$Y]
info_gene=data.frame(genes=row.names(plot_matrix_90),direction=NA)
info_gene$direction[info_gene$genes%in%MAG_genes_age_significant$Y[MAG_genes_age_significant$Beta>0]]=c("Positive")
info_gene$direction[info_gene$genes%in%MAG_genes_age_significant$Y[MAG_genes_age_significant$Beta<0]]=c("Negative")

# gene annotation
genes_groups <- data.frame(genes = row.names(plot_matrix_90),
                           direction = info_gene$direction[match(row.names(plot_matrix_90),info_gene$genes)])
genes_groups$direction[is.na(genes_groups$direction)]=c("not associated")
genes_groups$direction=factor(genes_groups$direction,levels = c("not associated","Positive","Negative"))
#metabolites_group=arrange(metabolites_group,age_direction,SuperClass)
#color: 34 kinds of species in total, 34 colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#set.seed(2004); direction_colors = sample(col_vector, length(unique(genes_groups$direction)), replace = F)
direction_colors=c("#cdd7d6","#235789","#c1292e")
names(direction_colors) <- unique(genes_groups$direction)
colours1 <- list(direction=direction_colors) # if species also has specific color: can add "species=gene_block_colors"; the generation method is the same
ha <- HeatmapAnnotation(df = genes_groups[,2,drop=F], show_legend = TRUE,which='row',col = colours1)
# MAG --> strain annotation
strain_groups <- data.frame(sample = colnames(plot_matrix_90),
                            strain = pcoa_with_clusters$pcoa_kmeans_cluster[match(colnames(plot_matrix_90),pcoa_with_clusters$ID)])
strain_groups$strain=factor(strain_groups$strain,levels = c("1","2","3"))
#strain_groups=arrange(strain_groups,strain)
strain_colors=c("#FF5575","#FFD36A","#6299FF")
names(strain_colors) <- unique(strain_groups$strain)
colours2 <- list(strain=strain_colors) # if species also has specific color: can add "species=gene_block_colors"; the generation method is the same
ha2 <- HeatmapAnnotation(df = strain_groups[,2,drop=F], show_legend = TRUE,which='column',col = colours2)

cluster_within = cluster_within_group(plot_matrix_90, 
                                      as.character(pcoa_with_clusters$pcoa_kmeans_cluster[match(colnames(plot_matrix_90),pcoa_with_clusters$ID)]))
CairoPNG("Rplot53.png", width = 700, height = 1000)
Heatmap(plot_matrix_90,
        name = "Beta",
        #row_km = 3, border = TRUE,
        #clustering_distance_rows = "spearman",
        #clustering_method_rows = "complete",
        #rect_gp = gpar(col = "white", lwd = 1.5),
        cluster_rows = T,  # 聚类Y轴
        show_row_names = T,
        cluster_columns = cluster_within,  # 聚类X轴
        show_column_names = F,
        top_annotation=ha2,
        #right_annotation = ha,
        col = colorRamp2(c(min(plot_matrix_90),0, max(plot_matrix_90)), c("#26547c", "white","#ef476f")),
        use_raster=F)
dev.off()

#### 2.2 volcano plot for MAG genes associated with age ####
O_ER4_MAG_gene_age_res <- readRDS("/scratch/p303998/SV_MWAS/O_ER4_MAG_gene_age_res.rds")
for (i in c("Y","X")){O_ER4_MAG_gene_age_res[,i]=as.character(O_ER4_MAG_gene_age_res[,i])}
length(unique(O_ER4_MAG_gene_age_res$Y))#1995
length(unique(O_ER4_MAG_gene_age_res$X))#1 phenotypes
sum(is.na(O_ER4_MAG_gene_age_res$p))
O_ER4_MAG_gene_age_res <- O_ER4_MAG_gene_age_res %>%
  mutate(significant = case_when(
    p < 0.05/(1995*1) & Beta >0 ~ "Positive",
    p < 0.05/(1995*1) & Beta <0 ~ "Negative",
    p > 0.05/(1995*1) ~ "non-significant"))
O_ER4_MAG_gene_age_res$significant=as.factor(O_ER4_MAG_gene_age_res$significant)
O_ER4_MAG_gene_age_res$annotation=NA
O_ER4_MAG_gene_age_res$annotation[grep("iotin",O_ER4_MAG_gene_age_res$product)]=O_ER4_MAG_gene_age_res$product[grep("iotin",O_ER4_MAG_gene_age_res$product)]
O_ER4_MAG_gene_age_res$annotation[grep("adenosylmethionine--8-amino-7-oxononanoate transaminase",O_ER4_MAG_gene_age_res$product)]=O_ER4_MAG_gene_age_res$product[grep("adenosylmethionine--8-amino-7-oxononanoate transaminase",O_ER4_MAG_gene_age_res$product)]
O_ER4_MAG_gene_age_res$annotation[O_ER4_MAG_gene_age_res$p>0.05/(1995*1)]=NA
ggplot(O_ER4_MAG_gene_age_res, aes(x =Beta, y= -log10(p), color=significant)) +
  geom_point(alpha=1, size=2.5) +
  scale_color_manual(values = c("non-significant" = "#cdd7d6", "Negative" = "#235789", "Positive" = "#c1292e")) +
  geom_vline(xintercept=c(-0.2,0.2),lty=4,col="gray",lwd=0.5)+
  geom_hline(yintercept = -log10(0.05/(1995*1)), lty=4,col="gray50",lwd=0.5) + 
  labs(x="Beta effect", y="-log10(P)") +
  theme(plot.title = element_text(hjust = 0.5),legend.position="right",legend.title = element_blank())+
  theme_prism(border = T)+
  geom_label_repel(
    aes(label = annotation),
    data=O_ER4_MAG_gene_age_res[!is.na(O_ER4_MAG_gene_age_res$annotation),],
    size = 4,  # Adjust label size for better visibility
    box.padding = 0.35,  # Adjust padding inside the label box (if needed)
    point.padding = 0.5,  # Adjust space between the point and the label
    segment.color = '#235789',
    min.segment.length = 0# Add segment lines from points to labels
  )

#### 2.3 biotin synthase BioB ####
lm_input<-data.frame(Y_gene = plot_matrix_90[sample_name,"UniRef90_A0A1Q6QCA2"], X_SV = pheno_cohort[sample_name,"Age"],pheno_cov[sample_name,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
print(lm_input$Y_gene)
# if (length(unique(lm_input$X_SV)) == 2){lm_input$X_SV=factor(lm_input$X_SV,levels=c("0","1"))}
# if (length(unique(lm_input$X_SV)) > 2){for (a in c("X_SV")){lm_input[,a]=qtrans(lm_input[,a])}}
if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
if ("DNA.Concentration" %in% names(lm_input)){
  for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
if ("log10_counts" %in% names(lm_input)){
  for (d in c("log10_counts")){lm_input[,d]=qtrans(lm_input[,d])}}
lm_res <- summary(glm(Y_gene~.,data = lm_input,family = "binomial"))
lm_input$Y_gene=as.factor(lm_input$Y_gene)
ggplot(lm_input, aes(x = Y_gene, y = X_SV)) +
  geom_flat_violin(aes(fill = Y_gene),trim = T,alpha=0.4,color = "white")+
  geom_point(aes(color = Y_gene),position = position_jitter(width = 0.13), size =1.3)+
  geom_boxplot(aes(fill = Y_gene),outlier.shape = NA, alpha = 0.5, width = .13,color = "black") +
  labs(x = "UniRef90_A0A1Q6QCA2\nbiotin synthase BioB", y = "Age (invrank-transformed)") +
  theme_classic()+
  theme(axis.text.y = element_text(face="plain",size=12,color="black"),
        axis.text.x = element_text(face="plain",size=12,color="black"),
        axis.title.x = element_text(color="black", size=14, face="plain"),
        axis.title.y = element_text(color="black", size=14, face="plain"))+
  scale_color_manual(values = c("0" = "#26547c", "1" = "#ef476f")) +
  scale_fill_manual(values = c("0" = "#26547c", "1" = "#ef476f")) 

#### 3. pangenome association: genes ~ strains ####
#pcoa_with_clusters <- readRDS("/scratch/p303998/SV_MWAS/Rdata_1216/pcoa_with_clusters.rds")
pcoa_with_clusters <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/pcoa_with_clusters.rds")
row.names(pcoa_with_clusters)=pcoa_with_clusters$ID
#plot_matrix <- readRDS("/scratch/p303998/SV_MWAS/plot_matrix.rds")
plot_matrix <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/plot_matrix.rds")
plot_matrix_90=plot_matrix[,colnames(plot_matrix)%in%c(mag_file_90$Bin.ID)] #2023 DAG3 MAGs
plot_matrix_90=plot_matrix_90[rowSums(plot_matrix_90)>0.05*ncol(plot_matrix_90)&rowSums(plot_matrix_90)<0.95*ncol(plot_matrix_90),] # Shell genes: 2022
plot_matrix_90=plot_matrix_90[,colSums(plot_matrix_90)>0.05*nrow(plot_matrix_90)]
plot_matrix_90=as.data.frame(plot_matrix_90)
#pangenome_mag <- readRDS("/scratch/p303998/SV_MWAS/pangenome_mag.rds")
pangenome_mag <- readRDS("~/Documents/SV_MWAS/R/Raw_data/MAG/pangenome_mag.rds")
plot_matrix_90=plot_matrix_90[,pangenome_mag]#one duplicated MAG
colnames(plot_matrix_90)=sapply(strsplit(as.character(colnames(plot_matrix_90)), "\\."), "[", 1)
plot_matrix_90=plot_matrix_90[,colnames(plot_matrix_90)%in%pcoa_with_clusters$ID] #1921 DAG3 MAGs
plot_matrix_90=t(plot_matrix_90)
# phenotype
cohort = c("dmp")
covar=c("DNA.Concentration", "Sex", "log10_counts")
sample_name = row.names(plot_matrix_90)
O_ER4_MAG_gene_strain_res=lm_gene_strain_pairwise(plot_matrix_90[sample_name,],pcoa_with_clusters[sample_name,"pcoa_kmeans_cluster",drop=F],pheno_cov[sample_name,],covar) 
saveRDS(O_ER4_MAG_gene_strain_res,paste0("/scratch/p303998/SV_MWAS/","O_ER4_MAG_gene_strain_res",".rds"))
O_ER4_MAG_gene_strain_res <- readRDS("/scratch/p303998/SV_MWAS/O_ER4_MAG_gene_strain_res.rds")
O_ER4_MAG_gene_strain_res$fdr.1=p.adjust(O_ER4_MAG_gene_strain_res$p.value_1,method = "fdr")
O_ER4_MAG_gene_strain_res$fdr.2=p.adjust(O_ER4_MAG_gene_strain_res$p.value_2,method = "fdr")
O_ER4_MAG_gene_strain_res$fdr.3=p.adjust(O_ER4_MAG_gene_strain_res$p.value_3,method = "fdr")
write.table(O_ER4_MAG_gene_strain_res, paste0("/scratch/p303998/SV_MWAS/Rdata_1216/","O_ER4_MAG_gene_strain_res",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
test=merge(plot_matrix_90[sample_name,],pcoa_with_clusters[sample_name,"pcoa_kmeans_cluster",drop=F],by.x = "row.names",by.y = "row.names",all=T) %>%`rownames<-`(.[,'Row.names']) %>% dplyr::select(-'Row.names')
for (i in 1:1995){
  test1=as.data.frame(table(test[,c(i,1996)]))
}