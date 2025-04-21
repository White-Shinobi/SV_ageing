# Generation of all_sv_info_anno_old.rds (the one combined of two annotations)
#### load data (generated before)####
#setwd("/scratch/p303998/SV_MWAS/Rdata_0828/SV_annotation/")
setwd("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/")
#load("/Users/helloduck/Resilio Sync/SV_MWAS_sync/Analysis/SV_anno/02.getAnno/RData/svAnnoDb.RData")
#load("/scratch/p303998/SV_MWAS/Rdata_0828/SV_annotation/svAnnoDb.RData")
#load("/scratch/p303998/SV_MWAS/Rdata_0828/SV_info/info.RData")
load("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/svAnnoDb.RData")
load("/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/info.RData")
tmp=svAnnoDb$gene.df.anno
test=svAnnoDb$bakta.full
svAnnoDb$bakta.full <- svAnnoDb$bakta.full %>%
  mutate(Uniref90 = str_extract(Attributes, "UniRef90_[^,;]+"))
colnames(svAnnoDb$bakta.full)[1]=colnames(svAnnoDb$gene.df.anno)[1]
svAnnoDb$gene.df.anno <- merge(svAnnoDb$gene.df.anno,svAnnoDb$bakta.full[,c("Taxid","Start","End","Uniref90")],by=c("Start","End","Taxid"),all.x = T)
svAnnoDb$gene.df.anno = svAnnoDb$gene.df.anno[,c(colnames(tmp),"Uniref90")]
tmp=svAnnoDb$gene.df.anno
#### Functions ####
library(tidyverse)
library(dplyr)
library(stringr)
library(biomaRt)
library(UniprotR)
library(KEGGREST)
library(segmenTools)
# flanking = 1, shuffling = F
getSvAnno<-function(inSv, geneAnno = gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 25, flanking = 1, shuffling = F){
  # inSv <- c("1519439.PRJNA253252:2001_2003")
  # geneAnno <- svAnnoDb[["gene.df.anno"]]
  # tax <- 1
  # geneName <- 2
  # start <- 4
  # end <- 5
  # anno.col <- 25
  # flanking <- 0
  # shuffling = F
  
  spe <- str_replace_all(inSv, ":.*", "")
  sv  <- str_replace_all(inSv, ".*:", "") %>% str_split(";")%>% unlist %>% str_split("_") %>%unlist%>%as.numeric %>%matrix(byrow = T, ncol = 2)
  
  sv[,1] <- sv[,1] - flanking
  sv[,2] <- sv[,2] + flanking
  
  geneAnno_spe <- geneAnno[geneAnno[,tax]==spe, ]
  
  if(shuffling==T){
    sv<-sv-min(sv)
    sample_limit <- trunc(max(geneAnno_spe[,end])/1000)-max(sv)
    #    shuffledStart <- sample(0:sample_limit, 1)
    shuffledStart <- trunc(runif(1, 0, sample_limit))
    sv<-sv+shuffledStart
  }
  
  gene <- NULL
  for (i in 1:nrow(sv)) {
    #    i<-1
    gene_i <- geneAnno_spe[sv[i,1]*1000<=geneAnno_spe[,end] & sv[i,2]*1000>=geneAnno_spe[,start], geneName]
    gene   <- c(gene, gene_i)
  }
  
  gene<-unique(gene)
  anno_return <- geneAnno[match(gene, geneAnno[,geneName]), anno.col]
  return(anno_return)
}


keggSet<-function(query.gene, background.gene, targetSetId, 
                  kegg_info = svAnnoDb$module_info, 
                  kegg_info.id.col = 1, 
                  kegg_info.ko.col = 6, 
                  geneKoTable = svAnnoDb$gene.df.anno,
                  geneKoTable.gene.col = 2, 
                  geneKoTable.ko.col = 27){
  # query.gene <- svs.info$Gene
  # background.gene <- gene.df.anno$GeneName[!is.na(gene.df.anno$KoNumber)]
  # targetSetId <- "M00001"
  # kegg_info <- module_info
  # kegg_info.id.col <- 1
  # kegg_info.ko.col <- 6
  # geneKoTable<-gene.df.anno
  # geneKoTable.gene.col <- 2
  # geneKoTable.ko.col <- 27
  
  targetSetKo <- kegg_info[match(targetSetId, kegg_info[,kegg_info.id.col]), kegg_info.ko.col] %>% str_split(";") %>% unlist
  background.ko <- match(background.gene, geneKoTable[,geneKoTable.gene.col]) %>% geneKoTable[., geneKoTable.ko.col]
  
  query.gene <- query.gene[query.gene %in% background.gene]
  geneInTarget <- background.gene[background.ko %in% targetSetKo]
  
  if(length(geneInTarget) > 1){
    geneNotQuery <- background.gene[! background.gene %in% query.gene]
    
    geneNotQuery.inTarget    <- geneNotQuery[ geneNotQuery %in% geneInTarget]
    geneNotQuery.notInTarget <- geneNotQuery[!geneNotQuery %in% geneInTarget]
    
    geneQuery.inTarget    <- query.gene[ query.gene %in% geneInTarget]
    geneQuery.notInTarget <- query.gene[!query.gene %in% geneInTarget]
    #    geneQuery.inTarget    <- background.gene[(  background.gene %in% geneInTarget) & (background.gene %in% query.gene)]
    #    geneQuery.notInTarget <- background.gene[(! background.gene %in% geneInTarget) & (background.gene %in% query.gene)]
    
    inputDf <- data.frame(geneQuery    = c(length(geneQuery.inTarget),    length(geneQuery.notInTarget)),
                          geneNotQuery = c(length(geneNotQuery.inTarget), length(geneNotQuery.notInTarget)))
    rownames(inputDf) <- c("geneInTarget", "geneNotInTarget")
    #inputDf
    fisher.res <- fisher.test(inputDf)
    enrich.res <- c(fisher.res$estimate, fisher.res$p.value, length(query.gene) ,length(geneQuery.inTarget), length(geneInTarget), length(background.gene))
  }else{
    geneQuery.inTarget <- NA
    enrich.res <- c(NA, NA, length(query.gene), NA, length(geneInTarget), length(background.gene))
  }
  
  names(enrich.res) <- c("OddsRatio", "P", "QueryN", "QueryInSetN", "SetN", "BackgroundN")
  
  return(list(enrich.res = enrich.res, geneQuery.inTarget = paste(geneQuery.inTarget, collapse = ";")))
}


allKeggSet <- function(query.gene, background.gene, 
                       kegg_info = svAnnoDb$module_info, 
                       kegg_info.id.col = 1, 
                       kegg_info.ko.col = 6, 
                       geneKoTable = svAnnoDb$gene.df.anno,
                       geneKoTable.gene.col = 2, 
                       geneKoTable.ko.col = 27){
  
  # query.gene <- svs.info$Gene
  # background.gene <- gene.df.anno$GeneName[!is.na(gene.df.anno$KoNumber)]
  # kegg_info <- module_info
  # kegg_info.id.col <- 1
  # kegg_info.ko.col <- 6
  # geneKoTable<-gene.df.anno
  # geneKoTable.gene.col <- 2
  # geneKoTable.ko.col <- 27
  
  keggEnrich.res <- matrix(NA, nrow = nrow(kegg_info), ncol = 8) %>% as.data.frame
  colnames(keggEnrich.res) <- c("SetId", "OddsRatio", "P", "QueryN", "QueryInSetN", "SetN", "BackgroundN","QueryInSet")
  for (i in 1:nrow(kegg_info)) {
    geneInTarget <- kegg_info[i, kegg_info.ko.col] %>% str_split(";") %>% unlist 
    keggEnrich.res_i <-keggSet(query.gene, background.gene, kegg_info[i, kegg_info.id.col],
                               kegg_info, 
                               kegg_info.id.col = 1, 
                               kegg_info.ko.col = 6, 
                               geneKoTable = geneKoTable,
                               geneKoTable.gene.col = 2, 
                               geneKoTable.ko.col = 27)
    
    keggEnrich.res$SetId[i]       <- kegg_info[i, kegg_info.id.col] 
    keggEnrich.res$OddsRatio[i]   <- keggEnrich.res_i$enrich.res[1]
    keggEnrich.res$P[i]           <- keggEnrich.res_i$enrich.res[2]
    keggEnrich.res$QueryN[i]      <- keggEnrich.res_i$enrich.res[3]
    keggEnrich.res$QueryInSetN[i] <- keggEnrich.res_i$enrich.res[4]
    keggEnrich.res$SetN[i]        <- keggEnrich.res_i$enrich.res[5]
    keggEnrich.res$BackgroundN[i] <- keggEnrich.res_i$enrich.res[6]
    keggEnrich.res$QueryInSet[i]  <- keggEnrich.res_i$geneQuery.inTarget
  }
  keggEnrich.res$FDR <- p.adjust(keggEnrich.res$P, method = "fdr")
  return(keggEnrich.res)
}

#### Add annotation information for vsv ####
vsv_info_anno <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_vsgv_info_anno.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) #3899 vSVs
# test block: vsv_info_anno=vsv_info_anno[200:210,]

vsv_info_anno$Gene <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 2, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
vsv_info_anno$Product.mixed <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 25, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
vsv_info_anno$KoNumber <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 27, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
vsv_info_anno$UniRef90 <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 28, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))

vsv_info_anno$Gene.bakta <-vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 10, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
vsv_info_anno$Product.bakta <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 12, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
vsv_info_anno$UniRef50.bakta <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 21, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
vsv_info_anno$UniRef90.bakta <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 23, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
vsv_info_anno$GeneSymbol <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 16, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
vsv_info_anno$HumanProteinInteraction <- vsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 22, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))

vsv_info_anno$SV_type=c("vSV")
saveRDS(vsv_info_anno, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/vsv_info_anno.rds")
#### Add annotation information for dsv ####
dsv_info_anno <- read.delim("/Users/helloduck/Documents/SV_MWAS/R/Raw_data/SV_info/20230827_full_dsgv_info_anno.tsv",sep="\t",row.names = NULL,header = T,check.names = F,fill = F) #10350 vSVs

dsv_info_anno$Gene <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 2, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
dsv_info_anno$Product.mixed <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 25, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
dsv_info_anno$KoNumber <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 27, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))
dsv_info_anno$UniRef90 <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$gene.df.anno, tax = 1, geneName = 2, start = 4, end = 5, anno.col = 28, flanking = 1, shuffling = F) %>%
  map_chr(~paste(.,collapse = "|"))

dsv_info_anno$Gene.bakta <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 10, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
dsv_info_anno$Product.bakta <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 12, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
dsv_info_anno$UniRef50.bakta <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 21, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
dsv_info_anno$UniRef90.bakta <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 23, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
dsv_info_anno$GeneSymbol <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 16, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))
dsv_info_anno$HumanProteinInteraction <- dsv_info_anno$SV_ID %>% as.list %>% lapply(getSvAnno, geneAnno = svAnnoDb$bakta.full, tax = 1, geneName = 10, start = 4, end = 5, anno.col = 22, flanking = 1, shuffling = F) %>%
map_chr(~paste(.,collapse = "|"))

dsv_info_anno$SV_type=c("dSV")
saveRDS(dsv_info_anno, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/dsv_info_anno.rds")

#### all SV annotation ####
all_sv_info_anno=rbind(vsv_info_anno,dsv_info_anno)
saveRDS(all_sv_info_anno, "/Users/helloduck/Documents/SV_MWAS/R/SV_annotation/all_sv_info_anno_old.rds")
