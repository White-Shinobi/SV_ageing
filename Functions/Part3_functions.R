library(scales)
library(reshape2)
#library(VennDetail)
library(Hmisc)
library(tidyverse)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(qqman)
library(dplyr)
library(ggprism)
library(ggforce)
library(ggrepel)
library(ggpie)
library(ggplot2)
library(PupillometryR)
#library(ggunchained)
library(ggpubr)
library(ggsci)
library(reshape2)              
#library(introdataviz)
#library(mediation)
#library("ggalluvial")

removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
removeRowsAllNa <- function(x) {x[which(apply(x, 1, function(y) any(!is.na(y)))), ]}
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}
#clr_transformation(row:taxa,column:sample)(Alex's version)
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  #######how to deal with zero
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  #geometric mean=exp(mean(log(x)))
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 2, gm_mean)#most important:Gmean_core should be based on all kinds of spp in each sample
  
  #do transformation
  data_prepared = rbind(Gmean_core,interest_matrix)
  data_transformed = apply(data_prepared,2,function(x){
    log(x / x[1])[-1]
  })
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

lm_abundance_age<-function(y_mat, x_mat, cov_mat, covar){
  
  sample_name=Reduce(intersect, list(row.names(y_mat), row.names(x_mat), row.names(cov_mat)))
  y_mat=y_mat[sample_name,,drop=F]
  x_mat=x_mat[sample_name,,drop=F]
  cov_mat=cov_mat[sample_name,covar,drop=F]
  y_x.edge<-NULL
  for (i in c(1:ncol(y_mat))){
    # i<-5
    cat(paste(i,colnames(y_mat)[i],"\n"))
    y_mat_i<-y_mat[,i,drop=F]#grep i not species name! grep species name will cause duplicated results!
    y_x_i <- lm_abundance_age_2(y_mat_i,x_mat, cov_mat, covar)
    y_x.edge<-rbind(y_x.edge,y_x_i)
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  return(y_x.edge)
}

lm_abundance_age_2<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  my_lm<-function(y,x){
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    older_18_sample_number    <-NA
    younger_18_sample_number    <-NA
    older_18_sample_rate <-NA
    younger_18_sample_rate <-NA
    
    #Y:SV,X:SO
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar,drop=F]) %>% na.omit(.) %>%as.data.frame(.)
    # lm_input<-data.frame(Y = sp_re_metaphlan_clr_EC_lost[sample_name,"t__SGB6571"], X = EC_pheno_use[sample_name,c("EC_lost")], pheno_cov[sample_name,covar]) %>% na.omit%>%as.data.frame
    colnames(lm_input)[1:2]=c("bacterial_abundance","phenotype")
    print(lm_input[1,])
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"bacterial_abundance"]))
    x_uniq_N   <- length(unique(lm_input[,"phenotype"]))
    
    older_18_sample_number <- sum(lm_input[,"phenotype"]>18,na.rm = T)
    younger_18_sample_number <- sum(!lm_input[,"phenotype"]>18,na.rm = T)

    older_18_sample_rate<-older_18_sample_number/N
    younger_18_sample_rate<-younger_18_sample_number/N
    
    if ("DNA.Concentration" %in% names(lm_input)){
      for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
    if ("log10_counts" %in% names(lm_input)){
      for (d in c("log10_counts")){lm_input[,d]=qtrans(lm_input[,d])}}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
    if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    
    if (length(unique(lm_input$bacterial_abundance))>1&length(unique(lm_input$phenotype))>1){
      try(lm_res <- summary(lm(bacterial_abundance~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        older_18_sample_number = older_18_sample_number,
                        younger_18_sample_number = younger_18_sample_number,
                        older_18_sample_rate = older_18_sample_rate,
                        younger_18_sample_rate= younger_18_sample_rate)),silent = T)
      }else{
        try(beta    <- NA,silent = T)
        try(se      <- NA, silent = T)
        try(p.value <- NA,silent = T)
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        older_18_sample_number = older_18_sample_number,
                        younger_18_sample_number = younger_18_sample_number,
                        older_18_sample_rate = older_18_sample_rate,
                        younger_18_sample_rate= younger_18_sample_rate)),silent = T)
      }
    }
    if (length(unique(lm_input$bacterial_abundance))<2&length(unique(lm_input$phenotype))<2){
        try(beta    <- NA,silent = T)
        try(se      <- NA, silent = T)
        try(p.value <- NA,silent = T)
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        older_18_sample_number = older_18_sample_number,
                        younger_18_sample_number = younger_18_sample_number,
                        older_18_sample_rate = older_18_sample_rate,
                        younger_18_sample_rate= younger_18_sample_rate)),silent = T)
    }
    
  }
  
  
  y_x<-sapply( 
    # y_mat=y_mat_i
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("bacterial_abundance", "phenotype", 
                        "beta","se","p.value",
                        "N","y_uniq_N","x_uniq_N", "older_18_sample_number", "younger_18_sample_number","older_18_sample_rate","younger_18_sample_rate")
  
  return(y_x_edge)
}

SV_Age_interaction<-function(pheno_mat, SV_mat, cov_mat_abun, covar, info,cohort){
  # pheno_mat<-pheno_3.1.1_Continuous[sample_name, setdiff(names(pheno_3.1.1_Continuous), covar)]
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov[sample_name,]
  # covar<-covar
  # abun<-sp_re_lld_baseline_clr[sample_name,]
  # info<-running_info_cohort
  
  y_x.edge<-NULL
  for (i in unique(info$V3)){
    # i<-unique(info$V3)[1]
    cat(paste(i,"\n"))
    spe_name<-str_replace_all(i,"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    SV_mat_i<-SV_mat[,grep(spe_name,colnames(SV_mat)),drop=F]
    pheno_mat_i<-pheno_mat[,colnames(pheno_mat)%in%info$V2[grep(spe_name,info$V3)],drop=F]
    if(dim(SV_mat_i)[2]>0&dim(pheno_mat_i)[2]>0){
      if(i%in%colnames(cov_mat_abun)){
        covar_i<-c(covar,i)
        y_x_i <- SV_Age_interaction_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        # Y is from SV_mat_i and X is from pheno_mat
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- SV_Age_interaction_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

SV_Age_interaction_2<-function(pheno_mat,SV_mat,cov_mat_abun,covar,cohort){
  require(reshape2)
  require(R.utils)
  
  my_lm<-function(y,x){
    # Note! Y is from SV_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-SV_mat[,1]
    
    beta.Age <- NA
    se.Age <- NA
    p.value.Age <- NA
    beta.SV <- NA
    se.SV <- NA
    p.value.SV <- NA
    beta.SV_Age <- NA
    se.SV_Age <- NA
    p.value.SV_Age <- NA
    anova_p <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(X_SV = x, Y_phenotype = y,cov_mat_abun[,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    # lm_input<-data.frame(X_SV = sgv_full[sample_name,c("Oscillibacter sp. ER4:2264_2265")], Y_phenotype = pheno_cohort[sample_name, "BMI"],pheno_cov_abundance[sample_name,c("DNA.Concentration","Age","Sex","log10_counts","Oscillibacter sp. ER4")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    # lm_input<-data.frame(X_SV = sgv_full[sample_name,c("Oscillibacter sp. ER4:2264_2265")], Y_phenotype = pheno_cohort[sample_name, "NoDisease"],pheno_cov_abundance[sample_name,c("DNA.Concentration","Age","Sex","log10_counts","Oscillibacter sp. ER4")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    if (dim(lm_input)[2]==7){colnames(lm_input)[7]=c("abundance")}
    print(lm_input$Y_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"X_SV"]))
    x_uniq_N   <- length(unique(lm_input[,"Y_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"X_SV"]==0)
    x_non_zero_N <- sum(!lm_input[,"Y_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$Y_phenotype)) == 2){lm_input$Y_phenotype=factor(lm_input$Y_phenotype,levels=c("0","1"))}
    if (length(unique(lm_input$Y_phenotype)) > 2){
      lm_input$Y_phenotype=as.numeric(lm_input$Y_phenotype)
      for (a in c("Y_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
    }
    
    if (length(unique(lm_input$X_SV)) > 2){for (b in c("X_SV")){lm_input[,b]=qtrans(lm_input[,b])}}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
    if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
    
    if ("DNA.Concentration" %in% names(lm_input)){
      for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
    
    if (length(unique(lm_input$Y_phenotype))>2){
      try(lm_res_2 <- summary(lm(Y_phenotype~X_SV*Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input)), silent = T)
      try(lm_res_1 <- summary(lm(Y_phenotype~X_SV+Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input)), silent = T)
      try(model_2 <- lm(Y_phenotype~X_SV*Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input), silent = T)
      try(model_1 <- lm(Y_phenotype~X_SV+Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input), silent = T)
      
      try(beta.Age    <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.Age      <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.Age <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.SV    <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.SV      <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.SV <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.SV_Age    <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),1],silent = T)
      try(se.SV_Age      <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),2], silent = T)
      try(p.value.SV_Age <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),4],silent = T)
      
      try(anova_p <- anova(model_2, model_1)$`Pr(>F)`[2],silent = T)
      try(if (is.null(anova_p)){anova_p <- NA},silent = T)
      
      
      try(return(list(beta.Age = beta.Age,
                      se.Age = se.Age,
                      p.value.Age = p.value.Age,
                      beta.SV = beta.SV,
                      se.SV = se.SV,
                      p.value.SV = p.value.SV,
                      beta.SV_Age = beta.SV_Age,
                      se.SV_Age = se.SV_Age,
                      p.value.SV_Age = p.value.SV_Age,
                      anova_p = anova_p,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=cohort)),
          silent = T)
      
    }
    if (length(unique(lm_input$Y_phenotype))==2){
      try(lm_res_2 <- summary(glm(Y_phenotype~X_SV*Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input,family = "binomial")), silent = T)
      try(lm_res_1 <- summary(glm(Y_phenotype~X_SV+Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input,family = "binomial")), silent = T)
      try(model_2 <- glm(Y_phenotype~X_SV*Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input,family = "binomial"), silent = T)
      try(model_1 <- glm(Y_phenotype~X_SV+Age+DNA.Concentration+Sex+log10_counts+abundance,data = lm_input,family = "binomial"), silent = T)
      
      try(beta.Age    <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.Age      <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.Age <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.SV    <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.SV      <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.SV <- lm_res_1$coefficients[match(c("X_SV"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.SV_Age    <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),1],silent = T)
      try(se.SV_Age      <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),2], silent = T)
      try(p.value.SV_Age <- lm_res_2$coefficients[match(c("X_SV:Age"),rownames(lm_res_2$coefficients)),4],silent = T)
      
      try(anova_p <- anova(model_2, model_1,test="Chisq")$`Pr(>Chi)`[2],silent = T)
      try(if (is.null(anova_p)){anova_p <- NA},silent = T)
      
      
      try(return(list(beta.Age = beta.Age,
                      se.Age = se.Age,
                      p.value.Age = p.value.Age,
                      beta.SV = beta.SV,
                      se.SV = se.SV,
                      p.value.SV = p.value.SV,
                      beta.SV_Age = beta.SV_Age,
                      se.SV_Age = se.SV_Age,
                      p.value.SV_Age = p.value.SV_Age,
                      anova_p = anova_p,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=cohort)),
          silent = T)
      
    }
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(pheno_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(SV_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 18, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(SV_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", 
                        "Beta.Age","SE.Age", "P.Age",
                        "Beta.SV","SE.SV", "p.SV",
                        "Beta.SV_Age","SE.SV_Age", "p.SV_Age",
                        "anova_p",
                        "N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}

calculate_cohens_f2 <- function(mod_full, mod_part){
  r_full <- summary(mod_full)$r.sq
  r_sub <- summary(mod_part)$r.sq
  f2 <- (r_full - r_sub)/(1-r_full)
  return(f2)
}

shortbred_interaction<-function(shortbred_data,pheno_mat, SV_mat, cov_mat_abun, covar, proteins){
  # pheno_mat<-pheno_3.1.1_Continuous[sample_name, setdiff(names(pheno_3.1.1_Continuous), covar)]
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov[sample_name,]
  # covar<-covar
  # abun<-sp_re_lld_baseline_clr[sample_name,]
  # info<- 
  
  y_x.edge<-NULL
  for (i in unique(proteins$SV_Name)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    # spe_name<-str_replace_all(spe_name,"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    SV_mat_i<-SV_mat[,i,drop=F]
    # shortbred_data=shortbred_all
    if(sum(colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i])>0){
      name=colnames(shortbred_data)[colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i]]
      shortbred_data_i=shortbred_data[,name,drop=F]
      covar_i<-c(covar)
      y_x_i <- shortbred_interaction_2(pheno_mat, shortbred_data_i, SV_mat_i, cov_mat_abun, covar_i)
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

shortbred_interaction_2<-function(pheno_mat,shortbred_mat,SV_mat_use,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  my_lm<-function(y,x){
    
    
    beta.Age <- NA
    se.Age <- NA
    p.value.Age <- NA
    beta.gene <- NA
    se.gene <- NA
    p.value.gene <- NA
    beta.gene_Age <- NA
    se.gene_Age <- NA
    p.value.gene_Age <- NA
    anova_p <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(X_shortbred = x, Y_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #lm_input<-data.frame(X_shortbred = shortbred_all[sample_name,"RHOM_15565"], Y_phenotype = pheno_cohort[sample_name, "NoDisease"],pheno_cov_abundance[sample_name,c("DNA.Concentration","Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$X_shortbred>0,]
    #### change all genes == 0 to smallest non-zero value/10
    lm_input$X_shortbred[lm_input$X_shortbred==0]=min(lm_input$X_shortbred[lm_input$X_shortbred>0])/10
    # lm_input<-data.frame(Y = vsgv_full[sample_name,c("Paraprevotella clara YIT 11840:3240_3241")], X = dag3_continuous[sample_name,c("hsCRP")],cov_mat_abun[sample_name,c("DNA.Concentration","Age","Sex","log10_counts","Paraprevotella clara YIT 11840")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input$Y_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"X_shortbred"]))
    x_uniq_N   <- length(unique(lm_input[,"Y_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"X_shortbred"]==0)
    x_non_zero_N <- sum(!lm_input[,"Y_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$Y_phenotype)) == 2){lm_input$Y_phenotype=factor(lm_input$Y_phenotype,levels=c("0","1"))}
    if (length(unique(lm_input$Y_phenotype)) > 2){
      lm_input$Y_phenotype=as.numeric(lm_input$Y_phenotype)
      for (a in c("Y_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
    }
    
    if (length(unique(lm_input$X_shortbred)) > 2){for (b in c("X_shortbred")){lm_input[,b]=log2(lm_input[,b])}}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
    if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
    if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
    
    if ("DNA.Concentration" %in% names(lm_input)){
      for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
    
    if (length(unique(lm_input$Y_phenotype))>2){
      try(lm_res_2 <- summary(lm(Y_phenotype~X_shortbred*Age+DNA.Concentration+Sex+log10_counts,data = lm_input)), silent = T)
      try(lm_res_1 <- summary(lm(Y_phenotype~X_shortbred+Age+DNA.Concentration+Sex+log10_counts,data = lm_input)), silent = T)
      try(model_2 <- lm(Y_phenotype~X_shortbred*Age+DNA.Concentration+Sex+log10_counts,data = lm_input), silent = T)
      try(model_1 <- lm(Y_phenotype~X_shortbred+Age+DNA.Concentration+Sex+log10_counts,data = lm_input), silent = T)
      
      try(beta.Age    <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.Age      <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.Age <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.gene    <- lm_res_1$coefficients[match(c("X_shortbred"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.gene      <- lm_res_1$coefficients[match(c("X_shortbred"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.gene <- lm_res_1$coefficients[match(c("X_shortbred"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.gene_Age    <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),1],silent = T)
      try(se.gene_Age      <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),2], silent = T)
      try(p.value.gene_Age <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),4],silent = T)
      
      try(anova_p <- anova(model_2, model_1)$`Pr(>F)`[2],silent = T)
      try(if (is.null(anova_p)){anova_p <- NA},silent = T)
      
      try(f2 <- calculate_cohens_f2(model_2, model_1),silent = T)
      try(if (is.null(f2)){f2 <- c("no value")},silent = T)
      
      try(return(list(beta.Age = beta.Age,
                      se.Age = se.Age,
                      p.value.Age = p.value.Age,
                      beta.gene = beta.gene,
                      se.gene = se.gene,
                      p.value.gene = p.value.gene,
                      beta.gene_Age = beta.gene_Age,
                      se.gene_Age = se.gene_Age,
                      p.value.gene_Age = p.value.gene_Age,
                      anova_p = anova_p,
                      f2 = f2,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=cohort)),
          silent = T)
    }
    
    if (length(unique(lm_input$Y_phenotype))==2){
      try(lm_res_2 <- summary(glm(Y_phenotype~X_shortbred*Age+DNA.Concentration+Sex+log10_counts,data = lm_input,family = "binomial")), silent = T)
      try(lm_res_1 <- summary(glm(Y_phenotype~X_shortbred+Age+DNA.Concentration+Sex+log10_counts,data = lm_input,family = "binomial")), silent = T)
      try(model_2 <- glm(Y_phenotype~X_shortbred*Age+DNA.Concentration+Sex+log10_counts,data = lm_input,family = "binomial"), silent = T)
      try(model_1 <- glm(Y_phenotype~X_shortbred+Age+DNA.Concentration+Sex+log10_counts,data = lm_input,family = "binomial"), silent = T)
      
      try(beta.Age    <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.Age      <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.Age <- lm_res_1$coefficients[match(c("Age"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.gene    <- lm_res_1$coefficients[match(c("Y_phenotype1"),rownames(lm_res_1$coefficients)),1],silent = T)
      try(se.gene      <- lm_res_1$coefficients[match(c("Y_phenotype1"),rownames(lm_res_1$coefficients)),2], silent = T)
      try(p.value.gene <- lm_res_1$coefficients[match(c("Y_phenotype1"),rownames(lm_res_1$coefficients)),4],silent = T)
      
      try(beta.gene_Age    <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),1],silent = T)
      try(se.gene_Age      <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),2], silent = T)
      try(p.value.gene_Age <- lm_res_2$coefficients[match(c("X_shortbred:Age"),rownames(lm_res_2$coefficients)),4],silent = T)
      
      try(anova_p <- anova(model_2, model_1,test="Chisq")$`Pr(>Chi)`[2],silent = T)
      try(if (is.null(anova_p)){anova_p <- NA},silent = T)
      
      try(f2 <- calculate_cohens_f2(model_2, model_1),silent = T)
      try(if (length(f2)==0){f2 <- NA},silent = T)
      
      try(return(list(beta.Age = beta.Age,
                      se.Age = se.Age,
                      p.value.Age = p.value.Age,
                      beta.gene = beta.gene,
                      se.gene = se.gene,
                      p.value.gene = p.value.gene,
                      beta.gene_Age = beta.gene_Age,
                      se.gene_Age = se.gene_Age,
                      p.value.gene_Age = p.value.gene_Age,
                      anova_p = anova_p,
                      f2 = f2,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=cohort)),
          silent = T)
    }
  }
  
  
  y_x<-sapply( 
    as.data.frame(pheno_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(shortbred_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 19, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(shortbred_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", 
                        "Beta.Age","SE.Age", "P.Age",
                        "Beta.gene","SE.gene", "p.gene",
                        "Beta.gene_Age","SE.gene_Age", "p.gene_Age",
                        "anova_p","f2",
                        "N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}

## Mediation analysis for linear model
my_lm_mediation<-function(input.inv, input.dv, input.med, covDf){
  #input.inv<-pheno_cohort[sample_name,"Age"] !!!Exposure!!!
  #input.dv <-pheno_cohort[sample_name,"Weight"] !!!Outcome!!!
  #input.med<-shortbred_all[sample_name,"AL1_12460"]
  #covDf<-pheno_cov[sample_name,c("DNA.Concentration", "Sex", "log10_counts")]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  
  if (length(unique(input.df.rmna$input.dv)) == 2){input.df.rmna$input.dv=factor(input.df.rmna$input.dv,levels=c("0","1"))}
  #if (length(unique(input.df.rmna$input.inv)) == 2){input.df.rmna$input.inv=factor(input.df.rmna$input.inv,levels=c("0","1"))}
  if (length(unique(input.df.rmna$input.dv)) > 2){
    input.df.rmna$input.dv=as.numeric(input.df.rmna$input.dv)
    for (a in c("input.dv")){input.df.rmna[,a]=qtrans(input.df.rmna[,a])}
  }
  #input.df.rmna$input.med[input.df.rmna$input.med==0]=min(input.df.rmna$input.med[input.df.rmna$input.med>0])/10
  #if (length(unique(input.df.rmna$input.med)) > 2){for (b in c("input.med")){input.df.rmna[,b]=log2(input.df.rmna[,b])}}
  if (length(unique(input.df.rmna$input.med)) > 2){for (b in c("input.med")){input.df.rmna[,b]=qtrans(input.df.rmna[,b])}}
  if (c("Sex")%in%colnames(input.df.rmna)){input.df.rmna$Sex=as.factor(input.df.rmna$Sex)}
  if (c("ETHNICITY")%in%colnames(input.df.rmna)){input.df.rmna$ETHNICITY=as.numeric(input.df.rmna$ETHNICITY)}
  if (c("MSM")%in%colnames(input.df.rmna)){input.df.rmna$MSM=as.numeric(input.df.rmna$MSM)}
  if (length(unique(input.df.rmna$Sex)) < 2){input.df.rmna$Sex=NULL}
  if ("DNA.Concentration" %in% names(input.df.rmna)){
    for (d in c("DNA.Concentration")){input.df.rmna[,d]=qtrans(input.df.rmna[,d])}}
  
  # Change here for binary_exp_vSV_continuous_out and continuous_exp_vSV_continuous_out !!!
  #input.df.rmna$input.inv=as.factor(input.df.rmna$input.inv)
  #input.df.rmna$Sex=as.factor(input.df.rmna$input.inv)
  # input.df.rmna$input.dv=as.factor(input.df.rmna$input.dv)
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      # change for continuous_exp_vSV_binary_out !!! 3 glm for input.dv
      # fit.totaleffect=glm(input.dv~.,input.df.rmna[,-3],family = "binomial")
      fit.totaleffect=lm(input.dv~.,data=input.df.rmna[,-grep("input.med",colnames(input.df.rmna))])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-grep("input.dv",colnames(input.df.rmna))])
      fit.mediator.res<-summary(fit.mediator)
      
      # fit.dv=glm(input.dv~.,input.df.rmna,family = "binomial") 
      fit.dv=lm(input.dv~.,input.df.rmna) 
      fit.dv.res<-summary(fit.dv)
      
      
      print("problem")
      print(covDf$ETHNICITY)
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
    }
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
  }
  
  return(res_list)
}

## Bidirectional mediation analysis for linear model
lm_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-running_candidate[1,]
  #indvDf<-pheno_cohort[sample_name,"Age",drop=F]
  #dvDf1<-shortbred_all[sample_name,]
  #dvDf2<-pheno_cohort[sample_name,]
  #covDf<-pheno_cov_abun[sample_name,]
  #covar<-c('Sex','log10_counts')
  # test block
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lm_mediation(indv, dv1, dv2, covDf[,covar])# age > gene > phenotype
  dir2_res <- my_lm_mediation(indv, dv2, dv1, covDf[,covar])# age > phenotype > gene
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "species_metabolites_hypertension"}#age_gene_phenotype
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "species_hypertension_metabolites"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}

## Mediation analysis for logistic regresion model
my_glm_mediation<-function(input.inv, input.dv, input.med, covDf){
  #input.inv<-pheno_cohort[sample_name,"Age"] !!!Exposure!!!
  #input.dv <-pheno_cohort[sample_name,"Weight"] !!!Outcome!!!
  #input.med<-shortbred_all[sample_name,"AL1_12460"]
  #covDf<-pheno_cov[sample_name,c("DNA.Concentration", "Sex", "log10_counts")]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  
  if (length(unique(input.df.rmna$input.med)) > 2){for (b in c("input.med")){input.df.rmna[,b]=qtrans(input.df.rmna[,b])}}
  
  #input.df.rmna$input.med[input.df.rmna$input.med==0]=min(input.df.rmna$input.med[input.df.rmna$input.med>0])/10
  #if (length(unique(input.df.rmna$input.med)) > 2){for (b in c("input.med")){input.df.rmna[,b]=log2(input.df.rmna[,b])}}
  if (c("Sex")%in%colnames(input.df.rmna)){input.df.rmna$Sex=as.factor(input.df.rmna$Sex)}
  if (length(unique(input.df.rmna$Sex)) < 2){input.df.rmna$Sex=NULL}
  if ("DNA.Concentration" %in% names(input.df.rmna)){
    for (d in c("DNA.Concentration")){input.df.rmna[,d]=qtrans(input.df.rmna[,d])}}
  if (c("ETHNICITY")%in%colnames(input.df.rmna)){input.df.rmna$ETHNICITY=as.numeric(input.df.rmna$ETHNICITY)}
  if (c("MSM")%in%colnames(input.df.rmna)){input.df.rmna$MSM=as.numeric(input.df.rmna$MSM)}
  
  # Change here for binary_exp_vSV_continuous_out and continuous_exp_vSV_continuous_out !!!
  input.df.rmna$input.inv=as.factor(input.df.rmna$input.inv)
  #input.df.rmna$Sex=as.factor(input.df.rmna$input.inv)
  # input.df.rmna$input.dv=as.factor(input.df.rmna$input.dv)
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=glm(input.dv~.,data=input.df.rmna[,-grep("input.med",colnames(input.df.rmna))],family = "binomial")
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-grep("input.dv",colnames(input.df.rmna))])
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=glm(input.dv~.,input.df.rmna,family = "binomial") 
      fit.dv.res<-summary(fit.dv)
      
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
    }
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
  }
  
  return(res_list)
}

## Bidirectional mediation analysis for logistic regresion model
glm_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-running_candidate[1,]
  #indvDf<-pheno_cohort[sample_name,"Age",drop=F]
  #dvDf1<-shortbred_all[sample_name,]
  #dvDf2<-pheno_cohort[sample_name,]
  #covDf<-pheno_cov_abun[sample_name,]
  #covar<-c('Sex','log10_counts')
  # test block
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lm_mediation(indv, dv1, dv2, covDf[,covar])# age > gene > phenotype
  dir2_res <- my_glm_mediation(indv, dv2, dv1, covDf[,covar])# age > phenotype > gene
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "age_gene_phenotype"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "age_phenotype_gene"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}












