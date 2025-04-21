#### 1. load data (generated before) ####
source(paste0("/Users/helloduck/Documents/SV_MWAS/R/R_scripts/Step2/", "Part1_functions.R"))
setwd("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/")
load("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/svAnnoDb.RData")
load("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/info.RData")
test=svAnnoDb$gene.df.anno
test2=svAnnoDb$bakta.full;test2 <- test2 %>%   mutate(RefSeq = str_extract(DbXref, "WP_\\d+\\.\\d+"))
tmp_1=svAnnoDb$gene.df.anno[,c("Taxid","Start","End","GeneName","Product.mixed","KoNumber.mixed")];tmp_1$Annotation_tool=c("Progenome1")
tmp_2=svAnnoDb$bakta.full[,c("Species","Start","End","GeneID","Product","AnnotationTool","Attributes")]
tmp_2 <- tmp_2 %>% mutate(Uniref90 = str_extract(Attributes, "UniRef90_[^,;]+"))
tmp_2 <- tmp_2 %>% mutate(KEGG = str_extract(Attributes, "KEGG[^,,]+"))
tmp_2$Attributes=NULL
colnames(tmp_2)[1]=colnames(tmp_1)[1]
#### 2. generate gene annotation file based on two methods  ####
all_annotation <- merge(tmp_1,tmp_2,by=c("Start","End","Taxid"),all = T)
all_annotation$GeneName[is.na(all_annotation$GeneName)]=all_annotation$GeneID[is.na(all_annotation$GeneName)]
all_annotation$Product.mixed[is.na(all_annotation$Product.mixed)]=all_annotation$Product[is.na(all_annotation$Product.mixed)]
all_annotation$Annotation_tool[is.na(all_annotation$Annotation_tool)]=all_annotation$AnnotationTool[is.na(all_annotation$Annotation_tool)]
all_annotation$KoNumber.mixed[is.na(all_annotation$KoNumber.mixed)]=all_annotation$KEGG[is.na(all_annotation$KoNumber.mixed)]
all_annotation$GeneID=NULL;all_annotation$Product=NULL;all_annotation$AnnotationTool=NULL;all_annotation$KEGG=NULL
saveRDS(all_annotation, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/all_annotation.rds")# 401,813
all_annotation <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_annotation.rds")
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
# step 1: get all genes predicted from progenome1 annotation files
gene_progenome1=strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)
# step 2: get all genes predicted from Bakta and other methods
gene_bakta=strsplit(as.character(all_sv_info_anno_old$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)
# step 3: combine them together and remove the replicated ones
gene_all=c(gene_progenome1,gene_bakta)%>%unique(.)
gene_all=gene_all[gene_all%in%all_annotation$GeneName]
length(gene_all)#113,521
test=all_annotation[all_annotation$GeneName%in%gene_all,]# 80,053 genes have been annotated to Uniref90: (80053/113521 = 70%)
saveRDS(test, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/bacterial_contribution_uniref90.rds")
#### 3. Table S8: 1369 genes and their annotation extracted from 105 replicable age-associated SVs. ####
#Progenome1_gene
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)%>%length(.)
strsplit(as.character(all_sv_info_anno_old$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)%>%length(.)
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
all_sv_info_anno_old=all_sv_info_anno_old[which(all_sv_info_anno_old$SV_Name%in%meta_age_significant$Y),]
Progenome1_gene=strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)
# bakta_genes
bakta_genes=strsplit(as.character(all_sv_info_anno_old$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)
bakta_genes=bakta_genes[bakta_genes%in%all_annotation$GeneName]#bakta_genes to run
bakta_genes_gff=test2[which(test2$GeneID%in%bakta_genes),]
bakta_genes_gff=bakta_genes_gff[which(bakta_genes_gff$Feature==c("CDS")),]#334
bakta_genes=bakta_genes_gff$GeneID
#Progenome1+Bakta genes
Progenome1_Bakta_genes=c(Progenome1_gene,bakta_genes)#1369
Progenome1_Bakta_genes_table=all_annotation[which(all_annotation$GeneName%in%Progenome1_Bakta_genes),]
Progenome1_Bakta_genes_table$SV=NA
for (i in Progenome1_Bakta_genes){
  #i=Progenome1_Bakta_genes[1]
  #i=c("OOIAGD_13255")
  print(i)
  if (length(grep(i,all_sv_info_anno_old$Gene))>0){
    SV_tmp=paste(all_sv_info_anno_old$SV_Name[grep(i,all_sv_info_anno_old$Gene)],collapse=";")
  }else{
    SV_tmp=paste(all_sv_info_anno_old$SV_Name[grep(i,all_sv_info_anno_old$Gene.bakta)],collapse=";")
  }
  Progenome1_Bakta_genes_table$SV[grep(i,Progenome1_Bakta_genes_table$GeneName)]=SV_tmp
}
# Double check
for (i in Progenome1_Bakta_genes){
  #i=c("OOIAGD_13255)
  print(i)
  if (length(grep(i,all_sv_info_anno_old$Gene))>0){
    print(grep(i,all_sv_info_anno_old$Gene))
  }else{
    print(grep(i,all_sv_info_anno_old$Gene.bakta))
  }
}  
# Shortbred marker number
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred.rds")
Final_results=Final_results[-grep("UniRef90",Final_results$Y),]
Final_results=Final_results[Final_results$Y%in%Progenome1_Bakta_genes_table$GeneName,]
unique(Final_results$Y)%>%length(.)

Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)#!!! 1111 genes (shortbred) !!!!
# Function
Progenome1_Bakta_genes_table$Product.mixed[grep("hypothetical protein",Progenome1_Bakta_genes_table$Product.mixed)]=c("Unknown function")
Progenome1_Bakta_genes_table$Product.mixed[is.na(Progenome1_Bakta_genes_table$Product.mixed)]=c("Unknown function")
Progenome1_Bakta_genes_table$KoNumber=Progenome1_Bakta_genes_table$KoNumber.mixed
Progenome1_Bakta_genes_table$KoNumber[grep(";",Progenome1_Bakta_genes_table$KoNumber.mixed)]=sapply(strsplit(as.character(Progenome1_Bakta_genes_table$KoNumber.mixed[grep(";",Progenome1_Bakta_genes_table$KoNumber.mixed)]),"\\;"),"[",1)
Progenome1_Bakta_genes_table$KoNumber[grep("=",Progenome1_Bakta_genes_table$KoNumber.mixed)]=sapply(strsplit(as.character(Progenome1_Bakta_genes_table$KoNumber.mixed[grep("=",Progenome1_Bakta_genes_table$KoNumber.mixed)]),"\\="),"[",2)
Progenome1_Bakta_genes_table$KoNumber=gsub("KEGG:","",Progenome1_Bakta_genes_table$KoNumber)
Progenome1_Bakta_genes_table$KoNumber[grep(";",Progenome1_Bakta_genes_table$KoNumber)]=sapply(strsplit(as.character(Progenome1_Bakta_genes_table$KoNumber[grep(";",Progenome1_Bakta_genes_table$KoNumber)]),"\\;"),"[",1)
# add module information
# Detailed KO information: https://www.genome.jp/kegg-bin/download_htext?htext=ko00001&format=htext
ko00001_htext2df=function(input_file, output_file,header_num = 0) {
  # Read the input file
  input <- readLines(input_file)
  
  # Open the output file
  output <- file(output_file, "w")
  
  count <- list()
  while (header_num > 0) {
    input <- input[-1]
    header_num <- header_num - 1
  }
  
  # Process the input file
  for (line in input) {
    if (grepl("^A", line)) {
      a <- line
      a <- sub(" ", "\t", a)
      next
    }
    if (grepl("^B  ", line)) {
      b <- line
      b <- sub("B  ", "", b)
      b <- sub(" ", "\t", b)
      next
    }
    if (grepl("^C    ", line)) {
      c <- line
      c <- sub("C    ", "", c)
      c <- sub(" ", "\t", c)
      next
    }
    if (grepl("^D      ", line)) {
      d <- line
      d <- sub("D      ", "", d)
      d <- sub("  ", "\t", d)
      writeLines(paste(a, b, c, d, sep = "\t"), output)
    }
  }
  
  # Close the output file
  close(output)
}
ko00001_htext2df("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/ko00001.keg","/Users/helloduck/Documents/SV_MWAS/NM_reviosion/kegg.file")
KO_info <- read.delim("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/kegg.file",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
Progenome1_Bakta_genes_table$KoNumber[!grepl("K", Progenome1_Bakta_genes_table$KoNumber)]=NA
Progenome1_Bakta_genes_table$Module=NA
for (i in 1:nrow(Progenome1_Bakta_genes_table)){
  if (!is.na(Progenome1_Bakta_genes_table$KoNumber[i])){
    print(i)
    if (length(grep(Progenome1_Bakta_genes_table$KoNumber[i], KO_info$V7))>0){
      Progenome1_Bakta_genes_table$Module[i]=KO_info$V4[grep(Progenome1_Bakta_genes_table$KoNumber[i], KO_info$V7)]}
  }
}
Progenome1_Bakta_genes_table$Module[grep("Unknown function", Progenome1_Bakta_genes_table$Product.mixed)]=c("Unknown function")
Progenome1_Bakta_genes_table$Module[grep("Poorly characterized", Progenome1_Bakta_genes_table$Module)]=c("Unknown function")
# add Text information from: https://doi.org/10.1038/s41586-019-1065-y; SGVFinder original paper!
Progenome1_Bakta_genes_table$Eran_text=NA
#transposon
Progenome1_Bakta_genes_table$Eran_text[grep("transpos\\S*|insertion|Tra[A-Z]|Tra[0-9]|IS[0-9]|conjugate transposon", Progenome1_Bakta_genes_table$Product.mixed)]=c("Transposon")
#plasmid
Progenome1_Bakta_genes_table$Eran_text[grep("relax\\S*|conjug\\S*|mob\\S*|plasmid|type IV|chromosome partitioning|chromosome segregation’", Progenome1_Bakta_genes_table$Product.mixed)]=c("Plasmid")
#phage
Progenome1_Bakta_genes_table$Eran_text[grep("apsid|phage|tail|head|tape measure|antiterminatio", Progenome1_Bakta_genes_table$Product.mixed)]=c("Phage")
#other HGT mechanisms
Progenome1_Bakta_genes_table$Eran_text[grep("ntegrase|excision\\S*|exonuclease|recomb|toxin|restrict\\S*|resolv\\S*|topoisomerase|reverse transcrip", Progenome1_Bakta_genes_table$Product.mixed)]=c("Other HGT")
#carbohydrate
Progenome1_Bakta_genes_table$Eran_text[grep("glycosyltransferase|glycoside hydrolase|xylan|monooxygenase|rhamnos\\S*|cellulose|sialidase|\\S*ose($|\\S|\\-)|acetylglucosaminidase|cellobiose|galact\\S*|fructose|aldose|starch|mannose|mannan\\S*|glucan|lyase|glycosyltransferase|glycosidase|pectin|SusD|SusC|fructokinase|galacto\\S*|arabino\\S*", Progenome1_Bakta_genes_table$Product.mixed)]=c("Carbohydrate")
#antibiotic resistance:
Progenome1_Bakta_genes_table$Eran_text[grep("azole resistance|antibiotic resistance|TetR|tetracycline resistance|VanZ|betalactam\\S*|beta-lactam|antimicrob\\S*|lantibio\\S*", Progenome1_Bakta_genes_table$Product.mixed)]=c("antibiotic resistence")
#DNA/RNA modification and transcription factor
Progenome1_Bakta_genes_table$Eran_text[grep("RNA|DNA|rna|igma|transcriptional regulator|Transcriptional regulator", Progenome1_Bakta_genes_table$Product.mixed)]=c("Protein families: genetic information processing")
#KO module + Text info
Progenome1_Bakta_genes_table$KOmodule_textinfo=Progenome1_Bakta_genes_table$Module
Progenome1_Bakta_genes_table$KOmodule_textinfo[which(is.na(Progenome1_Bakta_genes_table$KOmodule_textinfo))]=Progenome1_Bakta_genes_table$Eran_text[which(is.na(Progenome1_Bakta_genes_table$KOmodule_textinfo))]
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results_shortbred.rds")
Final_results=Final_results[-grep("UniRef90",Final_results$Y),]
Final_results=Final_results[Final_results$Y%in%Progenome1_Bakta_genes_table$GeneName,]
unique(Final_results$Y)%>%length(.)
Progenome1_Bakta_genes_table$genes_1115_with_representative_markers=NA
Progenome1_Bakta_genes_table$genes_1115_with_representative_markers[which(Progenome1_Bakta_genes_table$GeneName%in%Final_results$Y)]=c("Yes_1115")
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)#!!! 1111 genes (shortbred) !!!!
Progenome1_Bakta_genes_table$genes_1111_after_filtering=NA
Progenome1_Bakta_genes_table$genes_1111_after_filtering[which(Progenome1_Bakta_genes_table$GeneName%in%Final_results$Y)]=c("Yes_1111")
Progenome1_Bakta_genes_table$KOmodule_textinfo[which(is.na(Progenome1_Bakta_genes_table$KOmodule_textinfo))]=Progenome1_Bakta_genes_table$Product.mixed[which(is.na(Progenome1_Bakta_genes_table$KOmodule_textinfo))]
saveRDS(Progenome1_Bakta_genes_table,paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","Progenome1_Bakta_genes_table",".rds"))
write.table(Progenome1_Bakta_genes_table, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
#### 4.1 generate uniref90 file for extracting from Humann  ####
UniRef90_list=as.character(all_annotation$Uniref90[all_annotation$GeneName%in%gene_all])%>%unique(.)
UniRef90_list=UniRef90_list[!UniRef90_list==c("NA")]
UniRef90_list=UniRef90_list[!is.na(UniRef90_list)]
length(UniRef90_list)
saveRDS(UniRef90_list, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/UniRef90_list.rds")
write.table(UniRef90_list, paste0("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/","UniRef90_list",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
#### 4.2 Download Humann4 table of these UniRef90_list from cluster ####
#### 4.3 Humann4.0 generation (DMP)  ####
Humann4.0_DMP <- read.delim("~/Documents/SV_MWAS/NM_reviosion/filtered_output.txt",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
colnames(Humann4.0_DMP)=sapply(strsplit(as.character(colnames(Humann4.0_DMP)), "_"), function(x) paste(x[1:3], collapse = "_"))
colnames(Humann4.0_DMP)= gsub("_NA","",colnames(Humann4.0_DMP))
colnames(Humann4.0_DMP)= gsub("_kneaddata_paired","",colnames(Humann4.0_DMP))
colnames(Humann4.0_DMP)= gsub("_kneaddata","",colnames(Humann4.0_DMP))
row.names(Humann4.0_DMP)=Humann4.0_DMP$"# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs";Humann4.0_DMP$"# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs"=NULL
Humann4.0_DMP=as.data.frame(t(Humann4.0_DMP))
saveRDS(Humann4.0_DMP, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Humann4.0_DMP.rds")
#### 4.4 Humann4.0 generation (LLD1) ####
Humann4.0_LLD <- read.delim("~/Documents/SV_MWAS/NM_reviosion/filtered_output_LLD.txt",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
colnames(Humann4.0_LLD)= gsub("_kneaddata_paired_merged","",colnames(Humann4.0_LLD))
row.names(Humann4.0_LLD)=Humann4.0_LLD$"# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs"
Humann4.0_LLD$"# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs"=NULL
Humann4.0_LLD=as.data.frame(t(Humann4.0_LLD))
saveRDS(Humann4.0_LLD, "/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Humann4.0_LLD.rds")

#### 5. Generate all summary numbers [file: Keynote-NM_revision_2025_02_21] ####
#### 5.1 UniRef90_age_list ####
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
filter_zero_dmp <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/filter_zero_dmp.rds")
non_NA_final <- readRDS("~/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
SV_6447=non_NA_final$SV_Name[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]
all_sv_info_anno_6447=all_sv_info_anno_old[which(all_sv_info_anno_old$SV_Name%in%SV_6447),]
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
all_sv_info_anno_105=all_sv_info_anno_6447[which(all_sv_info_anno_6447$SV_Name%in%meta_age_significant$Y),]
all_annotation <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_annotation.rds")
test1=strsplit(as.character(all_sv_info_anno_105$Gene),"\\|")%>%unlist(.)%>%unique(.)
test2=strsplit(as.character(all_sv_info_anno_105$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)
test_all=c(test1,test2)%>%unique(.)
test_all=test_all[test_all%in%all_annotation$GeneName]
UniRef90_age_list=all_annotation$Uniref90[all_annotation$GeneName%in%test_all]%>%unique(.)
UniRef90_age_list=UniRef90_age_list[!is.na(UniRef90_age_list)]
length(UniRef90_age_list)
saveRDS(UniRef90_age_list, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/UniRef90_age_list.rds")
#### 5.2 UniRef90_not_list ####
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
filter_zero_dmp <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Results/filter_zero/filter_zero_dmp.rds")
non_NA_final <- readRDS("~/Documents/SV_MWAS/NM_reviosion/non_NA_final.rds")
SV_6447=non_NA_final$SV_Name[which(non_NA_final$After_filtering==c("passed")&non_NA_final$Pass_cutoff_in_how_many_other_cohorts>0)]
all_sv_info_anno_6447=all_sv_info_anno_old[which(all_sv_info_anno_old$SV_Name%in%SV_6447),]
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
all_sv_info_anno_6342=all_sv_info_anno_6447[which(!all_sv_info_anno_6447$SV_Name%in%meta_age_significant$Y),]
all_annotation <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_annotation.rds")
test1=strsplit(as.character(all_sv_info_anno_6342$Gene),"\\|")%>%unlist(.)%>%unique(.)
test2=strsplit(as.character(all_sv_info_anno_6342$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)
test_all=c(test1,test2)%>%unique(.)
test_all=test_all[test_all%in%all_annotation$GeneName]
UniRef90_not_list=all_annotation$Uniref90[all_annotation$GeneName%in%test_all]%>%unique(.)
UniRef90_not_list=UniRef90_not_list[!is.na(UniRef90_not_list)]
length(UniRef90_not_list)
UniRef90_not_list=setdiff(UniRef90_not_list,UniRef90_age_list)
length(UniRef90_not_list)
saveRDS(UniRef90_not_list, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/UniRef90_not_list.rds")
#### 5.3 humann number check ####
Humann4.0_DMP <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Humann4.0_DMP.rds")
colnames(Humann4.0_DMP)%in%UniRef90_age_list%>%sum(.)
colnames(Humann4.0_DMP)%in%UniRef90_not_list%>%sum(.)
number_file <- read.delim("~/Documents/SV_MWAS/NM_reviosion/number.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F)
number_file=number_file[grep("UniRef90",number_file$`# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs`),,drop=F]
number_file=number_file[-grep("\\|",number_file$`# Gene Family HUMAnN v4.0.0.alpha.1 Adjusted CPMs`),,drop=F]

#### 5.4 generate protein markers for shortbred  ####
all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)%>%length(.)
strsplit(as.character(all_sv_info_anno_old$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)%>%length(.)
meta_age_significant <- readRDS("~/Documents/SV_MWAS/R/Raw_data/meta_age_significant.rds")
meta_age_significant =meta_age_significant[,1:17]
all_sv_info_anno_old=all_sv_info_anno_old[which(all_sv_info_anno_old$SV_Name%in%meta_age_significant$Y),]
strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)%>%length(.)
#Progenome1_gene
write.table(data.frame(GeneID=strsplit(as.character(all_sv_info_anno_old$Gene),"\\|")%>%unlist(.)%>%unique(.)), paste0("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/","Progenome1_gene",".tsv"),sep = "\t", quote = F,col.names = T,row.names = F)
# bakta_genes
bakta_genes=strsplit(as.character(all_sv_info_anno_old$Gene.bakta),"\\|")%>%unlist(.)%>%unique(.)
bakta_genes=bakta_genes[bakta_genes%in%all_annotation$GeneName]#bakta_genes to run
bakta_genes_gff=test2[which(test2$GeneID%in%bakta_genes),]
bakta_genes_gff=bakta_genes_gff[which(bakta_genes_gff$Feature==c("CDS")),]#334! Only CDS can be run
#### !!!! starting position is different based on different tools ####
bakta_genes_gff$Start=bakta_genes_gff$Start+1;bakta_genes_gff$End=bakta_genes_gff$End+1 #less 1121098.PRJDB570.shortbred.final.faa | grep X| wc -l ---> 22
#### !!!! ####
write.table(bakta_genes_gff[which(bakta_genes_gff$Feature==c("CDS")),1:9], paste0("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/","Bakta_gene",".tsv"),sep = "\t", quote = F,col.names = F,row.names = F)
write.table(data.frame(unique(bakta_genes_gff$Species)), paste0("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/","Species",".tsv"),sep = "\t", quote = F,col.names = ,row.names = F)

#### Table S. 40,280 bacterial gene families and their annotation extracted from 6,342 non-age-associated SVs.####
Final_results <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final_results.rds")
unique(Final_results$Y)%>%length(.)
UniRef90_not_list <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/UniRef90_not_list.rds")
length(UniRef90_not_list)
Final_results=Final_results[which(Final_results$Y%in%UniRef90_not_list),]
unique(Final_results$Y)%>%length(.)

all_sv_info_anno_old <- readRDS("~/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")##all_sv_info_anno_old.rds is the one combined of two annotations
UniRef90_info=data.frame(UniRef90=as.character(unique(Final_results$Y)),SV=NA)
for (i in as.character(unique(Final_results$Y))){
  #i=as.character(unique(Final_results$Y))[1]
  print(i)
  print(grep(i,UniRef90_info$UniRef90))
  UniRef90_info$SV[grep(i,UniRef90_info$UniRef90)]=paste0(all_sv_info_anno_old$SV_Name[grep(i,all_sv_info_anno_old$UniRef90.bakta)],collapse = ";")
}
Final_results$Beta_IDage=NULL;Final_results$SE_IDage=NULL;Final_results$p_IDage=NULL
Final_results=na.omit(Final_results)
unique(Final_results$Y)%>%length(.)#39791 genes
UniRef90_info$Included_in_downstream_analysis=NA
UniRef90_info$Included_in_downstream_analysis[which(UniRef90_info$UniRef90%in%unique(Final_results$Y))]=c("Yes")
saveRDS(UniRef90_info, paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","UniRef90_40280_info",".rds"))
write.table(UniRef90_info, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)
Unire90_GO <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/Unire90_GO.rds")
UniRef90_info=merge(UniRef90_info,Unire90_GO,by.x = "UniRef90",by.y = "Uniref90",all.x = T)
GO_name <- read.delim("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/Final R scripts/full_mapping_v4_alpha/map_go_name.txt",sep="\t",row.names = NULL,header = F,check.names = F,fill = F)
UniRef90_info=merge(UniRef90_info,GO_name,by.x = "GO",by.y = "V1",all.x = T)
colnames(UniRef90_info)[5]=c("GO_name")
Uniref90_name <- readRDS("~/Documents/SV_MWAS/NM_reviosion/Final R scripts/Uniref90_name.rds")
UniRef90_info=merge(UniRef90_info,Uniref90_name,by.x = "UniRef90",by.y = "V1",all.x = T)
UniRef90_info[is.na(UniRef90_info)]=c("Unknown")
length(unique(UniRef90_info$UniRef90))#40280
saveRDS(UniRef90_info, paste0("/Users/helloduck/Documents/SV_MWAS/NM_reviosion/","UniRef90_40280_info_annotation",".rds"))
write.table(UniRef90_info, "/Users/helloduck/Desktop/table1.tsv",sep = "\t", col.names = T, row.names = F, quote = F)


