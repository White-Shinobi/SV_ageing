library(tidyverse)
library(doParallel)
library(foreach)
.libPaths(c("/home3/p303998/R/x86_64-pc-linux-gnu-library/4.3", .libPaths()))
library(ppcor)
library(scales)
#library(VennDetail)



removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
removeRowsAllNa <- function(x) {x[which(apply(x, 1, function(y) any(!is.na(y)))), ]}
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
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

SV_lm_glm<-function(pheno_mat, SV_mat, cov_mat_abun, covar, info,cohort){
  # pheno_mat<-pheno_3.1.1_Continuous[sample_name, setdiff(names(pheno_3.1.1_Continuous), covar)]
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov[sample_name,]
  # covar<-covar
  # abun<-sp_re_lld_baseline_clr[sample_name,]
  # info<-info_3.1.1
  
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
        y_x_i <- SV_lm_glm_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        # Y is from SV_mat_i and X is from pheno_mat
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- SV_lm_glm_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

SV_lm_glm_2<-function(pheno_mat,SV_mat,cov_mat_abun,covar,cohort){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #SV_mat<- SV_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar
  #covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from SV_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-SV_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_SV = x, X_phenotype = y,cov_mat_abun[,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    # lm_input<-data.frame(Y = vsgv_full[sample_name,c("Paraprevotella clara YIT 11840:3240_3241")], X = dag3_continuous[sample_name,c("hsCRP")],cov_mat_abun[sample_name,c("DNA.Concentration","Age","Sex","log10_counts","Paraprevotella clara YIT 11840")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input$X_phenotype)
    print(colnames(lm_input))
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_SV"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_SV"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
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
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))>2){
      try(lm_res <- summary(lm(Y_SV~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
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
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))==2){
      try(lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_SV))<2){
      
      try(beta    <- c("one level"),silent = T)
      try(se      <- c("one level"), silent = T)
      try(p.value <- NA,silent = T)
      
      
      try(return(list(beta = beta,
                      se = se,
                      p.value = p.value,
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
  y_x.unlist <- matrix(unlist(y_x), ncol = 11, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(SV_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

SV_lm_glm_age2<-function(pheno_mat, SV_mat, cov_mat_abun, covar, info,cohort){
  # pheno_mat<-pheno_3.1.1_Continuous[sample_name, setdiff(names(pheno_3.1.1_Continuous), covar)]
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov[sample_name,]
  # covar<-covar
  # abun<-sp_re_lld_baseline_clr[sample_name,]
  # info<-info_3.1.1
  
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
        y_x_i <- SV_lm_glm_age2_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        # Y is from SV_mat_i and X is from pheno_mat
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- SV_lm_glm_age2_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

SV_lm_glm_age2_2<-function(pheno_mat_use,SV_mat_use,cov_mat_abun,covar,cohort){
  require(reshape2)
  require(R.utils)
  
  my_lm<-function(y,x){
    # Note! Y is from SV_mat_i and X is from pheno_mat_use. Cause our model is treating SV as Y!
    # y<-pheno_mat_use[,1]
    #  x<-SV_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    p_anova         <-NA
    
    lm_input<-data.frame(Y_SV = x, X_phenotype = y,cov_mat_abun[,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    if (c("Age")%in%colnames(lm_input)){lm_input$Age_2=(lm_input$Age)^2}
    # print(lm_input$X_phenotype[1])
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_SV"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_SV"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
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
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))>2){
      try(lm_res <- summary(lm(Y_SV~.,data = lm_input)), silent = T)
      try({model1<-lm(Y_SV~.,data = lm_input[, !names(lm_input) %in% "Age_2"])}, silent = T)
      try({model2<-lm(Y_SV~.,data = lm_input)}, silent = T)
      try(anova_age<-anova(model1,model2)$`Pr(>F)`[2], silent = T)
      
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        if (exists("anova_age")){try(p_anova <- anova_age,silent = T)}else{try(p_anova <- NA,silent = T)}
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=Cohort,
                        p_anova=p_anova)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        try(p_anova <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=Cohort,
                        p_anova=p_anova)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))==2){
      try(lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial")), silent = T)
      try({model1 <- glm(Y_SV ~ ., data = lm_input[, !names(lm_input) %in% "Age_2"], family = "binomial")}, silent = TRUE)
      try({model2<-glm(Y_SV~.,data = lm_input,family = "binomial")}, silent = T)
      try(anova_age<-anova(model1,model2, test = "Chisq")$`Pr(>Chi)`[2], silent = T)
      
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        if (exists("anova_age")){try(p_anova <- anova_age,silent = T)}else{try(p_anova <- NA,silent = T)}
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=Cohort,
                        p_anova=p_anova)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        try(p_anova <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=Cohort,
                        p_anova=p_anova)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_SV))<2){
      
      try(beta    <- c("one level"),silent = T)
      try(se      <- c("one level"), silent = T)
      try(p.value <- NA,silent = T)
      try(p_anova <- NA,silent = T)
      
      
      try(return(list(beta = beta,
                      se = se,
                      p.value = p.value,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=Cohort,
                      p_anova=p_anova)),
          silent = T)
    }
    if (length(lm_input)<2){
      
      try(beta    <- c("lm_input_all_NA"),silent = T)
      try(se      <- c("lm_input_all_NA"), silent = T)
      try(p.value <- NA,silent = T)
      try(p_anova <- NA,silent = T)
      
      
      try(return(list(beta = beta,
                      se = se,
                      p.value = p.value,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Cohort=Cohort,
                      p_anova=p_anova)),
          silent = T)
    }
  }
  
  y_x<-sapply( 
    as.data.frame(pheno_mat_use),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(SV_mat_use)
    )
  )
  
  # saveRDS(y_x, paste0("~/Downloads/2023_09_02/", "test", ".rds")) This is the debug method: table(names(unlist(test))); test3=data.frame(num=seq(1:length(unlist(test))),name=names(unlist(test)),content=unlist(test)); test3$check=test3$num/12; test3=test3[test3$num %% 12 != 0,]; test1[[70]][[92]]
  
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat_use), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat_use)
  rownames(y_x.beta)<-colnames(SV_mat_use)
  y_x_edge_beta  <- melt(y_x.beta)
  
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  colnames(y_x_edge)<-c("X", "Y", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort","p_anova_for_Age2","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}




shortbred_correct<-function(shortbred_data,pheno_mat, SV_mat, cov_mat_abun, covar, proteins){
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
    spe_name<-unique(proteins$organism[proteins$SV_Name==i])
    # spe_name<-str_replace_all(spe_name,"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    SV_mat_i<-SV_mat[,i,drop=F]
    # shortbred_data=shortbred_all
    if(sum(colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i])>0){
      name=colnames(shortbred_data)[colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i]]
      shortbred_data_i=shortbred_data[,name,drop=F]
      covar_i<-c(covar,spe_name)
      y_x_i <- shortbred_lm_glm(pheno_mat, shortbred_data_i, SV_mat_i, cov_mat_abun, covar_i)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

shortbred_lm_glm<-function(pheno_mat,shortbred_mat,SV_mat_use,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  my_lm<-function(y,x){
    # Note! Y is from shortbred_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-shortbred_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_shortbred = x, X_phenotype = y,SV = SV_mat_use,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    colnames(lm_input)[3]=c("SV")
    print(lm_input)
    ####!!!! Change when you run the binary data!!!!
    lm_input<-lm_input[lm_input$Y_shortbred>0,]
    # lm_input<-data.frame(Y = vsgv_full[sample_name,c("Paraprevotella clara YIT 11840:3240_3241")], X = dag3_continuous[sample_name,c("hsCRP")],cov_mat_abun[sample_name,c("DNA.Concentration","Age","Sex","log10_counts","Paraprevotella clara YIT 11840")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_shortbred"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_shortbred"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
    if (length(unique(lm_input$X_phenotype)) > 2){
      lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
      for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
    }
    
    if (length(unique(lm_input$Y_shortbred)) > 2){for (b in c("Y_shortbred")){lm_input[,b]=log2(lm_input[,b])}}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
    if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
    if (c("SV")%in%colnames(lm_input)){for (c in c("SV")){lm_input[,c]=qtrans(lm_input[,c])}}
    
    if ("DNA.Concentration" %in% names(lm_input)){
      for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_shortbred))>2){
      try(lm_res <- summary(lm(Y_shortbred~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
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
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_shortbred))==2){
      try(lm_res <- summary(glm(Y_shortbred~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_shortbred))<2){
      
      try(beta    <- c("one level"),silent = T)
      try(se      <- c("one level"), silent = T)
      try(p.value <- NA,silent = T)
      
      
      try(return(list(beta = beta,
                      se = se,
                      p.value = p.value,
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
  y_x.unlist <- matrix(unlist(y_x), ncol = 11, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(shortbred_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

genes_correct_noSV<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run){
  # pheno_mat<-pheno_cohort
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      gene_data_i=gene_data[,name,drop=F]
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- genes_lm_glm_noSV(pheno_mat, gene_data_i, cov_mat_abun, covar)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

genes_lm_glm_noSV<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    if (length(covar_use)==1){lm_input<-data.frame(Y_genes = x, X_phenotype = y) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}else{lm_input<-data.frame(Y_genes = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(lm_input[1:10,])
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
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
    
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_genes))<2){
      
      try(beta_IDage    <- c("no lm_res"),silent = T)
      try(se_IDage      <- c("no lm_res"), silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("no lm_res"),silent = T)
      try(se_phenotype      <- c("no lm_res"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
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
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}


shortbred_ppcor1<-function(shortbred_data,pheno_mat, SV_mat, cov_mat_abun, covar, proteins){
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
    spe_name<-unique(proteins$organism[proteins$SV_Name==i])
    # spe_name<-str_replace_all(spe_name,"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    SV_mat_i<-SV_mat[,i,drop=F]
    # shortbred_data=shortbred_all
    if(sum(colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i])>0){
      name=colnames(shortbred_data)[colnames(shortbred_data)%in%proteins$GeneName[proteins$SV_Name==i]]
      shortbred_data_i=shortbred_data[,name,drop=F]
      covar_i<-c(covar)
      y_x_i <- shortbred_ppcor2(pheno_mat, shortbred_data_i, SV_mat_i, cov_mat_abun, covar_i)
      y_x_i$Adjust_Abundance<-"No"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

shortbred_ppcor2<-function(pheno_mat,shortbred_mat,SV_mat_use,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #shortbred_mat<- shortbred_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from shortbred_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-shortbred_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_shortbred = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_shortbred>0,]
    #### change all genes == 0 to smallest non-zero value/10
    lm_input$Y_shortbred[lm_input$Y_shortbred==0]=min(lm_input$Y_shortbred[lm_input$Y_shortbred>0])/10
    # lm_input<-data.frame(Y_shortbred = shortbred_all[sample_name,"BACUNI_00888"], X_phenotype = pheno_cohort[sample_name, "BMI"],pheno_cov_abundance[sample_name,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_shortbred"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_shortbred"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$X_phenotype)) > 2){
      lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
    }
    
    for (b in c("Y_shortbred")){lm_input[,b]=as.numeric(lm_input[,b])}
    for (b in c("Y_shortbred")){lm_input[,b]=log2(lm_input[,b])}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.numeric(lm_input$Sex)}
    # if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_shortbred))>1){
      #partial correlation
      try(lm_res <- pcor.test(x = lm_input$X_phenotype, y = lm_input$Y_shortbred, z = lm_input[, !colnames(lm_input)%in%c("X_phenotype","Y_shortbred")], method = 'spearman'), silent = T)
      if (exists("lm_res")){
        try(beta    <- lm_res$estimate,silent = T)
        try(se      <- lm_res$statistic, silent = T)
        try(p.value <- lm_res$p.value,silent = T)
        
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta    <- c("no lm_res"),silent = T)
        try(se      <- c("no lm_res"), silent = T)
        try(p.value <- NA,silent = T)
        
        try(return(list(beta = beta,
                        se = se,
                        p.value = p.value,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_shortbred))<2){
      
      try(beta    <- c("one level"),silent = T)
      try(se      <- c("one level"), silent = T)
      try(p.value <- NA,silent = T)
      
      
      try(return(list(beta = beta,
                      se = se,
                      p.value = p.value,
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
  y_x.unlist <- matrix(unlist(y_x), ncol = 11, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(shortbred_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Rho_estimate","statistic", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, se_col) {
  require(meta)
  require(metafor)
  
  #inVec<-test[1,]
  #study_name<-c("DMP", "300TZFG")
  #beta_col<-c(1,3)
  #se_col<-c(2,4)
  
  study_beta <- inVec[beta_col] %>% as.numeric
  study_se <- inVec[se_col] %>% as.numeric
  
  #study_beta<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  
  inDf <- data.frame(study_name, study_beta, study_se)
  
  m.hksj <- metagen(
    study_beta,
    study_se,
    data = inDf,
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.beta = m.hksj.res$random$TE,
      Meta.se = m.hksj.res$random$seTE,
      Meta.p = m.hksj.res$random$p,
      Meta.I2 = m.hksj.res$I2,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}

## Batch meta-analysis
my_batch_meta_lm <- function(inDf,study_name,beta_col,se_col,row_var_col = 1,col_var_col = 2) {
  #inDf<-cbind_vsv_ba_lm_edge[c(1:10),]
  #study_name<-c("LLD", "300OB")
  #beta_col<-c(3,10)
  #se_col<-c(4,11)
  #row_var_col<-1
  #col_var_col<-2
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_lm,
    study_name = study_name,
    beta_col = beta_col,
    se_col = se_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 5, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.beta', 'Meta.SE', "Meta.p", "Meta.I2", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  rowName <- unique(inDf[, row_var_col])
  colName <- unique(inDf[, col_var_col])
  
  N_row <- length(rowName)
  N_col <- length(colName)
  
  batch_res.beta <-
    matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
  colnames(batch_res.beta) <- colName
  rownames(batch_res.beta) <- rowName
  batch_res.beta[is.na(batch_res.beta)] <- 0
  
  batch_res.p <- matrix(batch_res.unlist[, 3], ncol = N_col, byrow = F)
  colnames(batch_res.p) <- colName
  rownames(batch_res.p) <- rowName
  batch_res.p[is.na(batch_res.p)] <- 0
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  # fdr matrix
  batch_res.fdr <- matrix(
    data  = batch_res_edge$Meta.fdr.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = F
  )
  colnames(batch_res.fdr) <- colName
  rownames(batch_res.fdr) <- rowName
  batch_res.fdr[is.na(batch_res.fdr)] <- 1
  
  # bonferroni matrix
  batch_res.bon <- matrix(
    data  = batch_res_edge$Meta.bonferroni.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = F
  )
  
  colnames(batch_res.bon) <- colName
  rownames(batch_res.bon) <- rowName
  batch_res.bon[is.na(batch_res.bon)] <- 1
  
  return(
    list(
      table      = batch_res_edge,
      beta       = batch_res.beta,
      p          = batch_res.p,
      fdr        = batch_res.fdr,
      bonferroni = batch_res.bon
    )
  )
}

SV_lm_glm_abundance <- function(pheno_mat, SV_mat, cov_mat_abun, covar, info,cohort){
  # pheno_mat<-pheno_3.1.1_Continuous[sample_name, se_phenotypetdiff(names(pheno_3.1.1_Continuous), covar)]
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov[sample_name,]
  # covar<-covar
  # abun<-sp_re_lld_base_phenotypeline_clr[sample_name,]
  # info<-info_3.1.1
  
  y_x.edge<-NULL
  for (i in unique(info$V3)){
    # i<-unique(running_list_hiv$V3)[1]
    cat(paste(i,"\n"))
    spe_name<-str_replace_all(i,"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    SV_mat_i<-SV_mat[,grep(spe_name,colnames(SV_mat)),drop=F]
    pheno_mat_i<-pheno_mat[,colnames(pheno_mat)%in%info$V2[grep(spe_name,info$V3)],drop=F]
    #print(SV_mat_i[1:5,,drop=F])
    print(pheno_mat_i[1:5,,drop=F])
    if(dim(SV_mat_i)[2]>0&dim(pheno_mat_i)[2]>0){
      if(i%in%colnames(cov_mat_abun)){
        covar_i<-c(covar,i)
        y_x_i <- SV_lm_glm_abundance_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        # Y is from SV_mat_i and X is from pheno_mat
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- SV_lm_glm_abundance_2(pheno_mat_i, SV_mat_i, cov_mat_abun, covar_i,cohort)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p.value_phenotype,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p.value_phenotype,method = 'bonferroni')
  
  return(y_x.edge)
}

SV_lm_glm_abundance_2 <- function(pheno_mat,SV_mat,cov_mat_abun,covar,cohort){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #SV_mat<- SV_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar
  #covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from SV_mat_i and X is from pheno_mat. Cause_phenotype our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-SV_mat[,1]
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    beta_abundance    <- NA
    se_abundance      <- NA
    p.value_abundance <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_SV = x, X_phenotype = y,cov_mat_abun[,covar]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    # lm_input<-data.frame(Y = vsgv_full[sample_name,c("Paraprevotella clara YIT 11840:3240_3241")], X = dag3_continuous[sample_name,c("hsCRP")],cov_mat_abun[sample_name,c("DNA.Concentration","Age","sex","log10_counts","Paraprevotella clara YIT 11840")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    print(lm_input$X_phenotype)
    print(colnames(lm_input))
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_SV"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_SV"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
    if (length(unique(lm_input$X_phenotype)) > 2){
      lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
      for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
    }
    
    if (length(unique(lm_input$Y_SV)) > 2){for (b in c("Y_SV")){lm_input[,b]=qtrans(lm_input[,b])}}
    if (c("sex")%in%colnames(lm_input)){lm_input$sex=as.factor(lm_input$sex)}
    if (length(unique(lm_input$sex)) < 2){lm_input$sex=NULL}
    if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
    
    if ("DNA.Concentration" %in% names(lm_input)){
      for (d in c("DNA.Concentration")){lm_input[,d]=qtrans(lm_input[,d])}}
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))>2){
      try(lm_res <- summary(lm(Y_SV~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        indv_abundance<-colnames(lm_input)[ncol(lm_input)]
        print(indv_abundance)
        
        try(beta_phenotype    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        try(beta_abundance    <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),1],silent = T)
        try(se_abundance      <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_abundance <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        beta_abundance = beta_abundance,
                        se_abundance = se_abundance,
                        p.value_abundance = p.value_abundance,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(beta_abundance    <- c("no lm_res"),silent = T)
        try(se_abundance      <- c("no lm_res"), silent = T)
        try(p.value_abundance <- NA,silent = T)
        
        try(return(list(beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        beta_abundance = beta_abundance,
                        se_abundance = se_abundance,
                        p.value_abundance = p.value_abundance,
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
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_SV))==2){
      try(lm_res <- summary(glm(Y_SV~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        indv_abundance<-colnames(lm_input)[ncol(lm_input)]
        print(indv_abundance)
        
        try(beta_phenotype    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
        
        try(beta_abundance    <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),1],silent = T)
        try(se_abundance      <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_abundance <- lm_res$coefficients[match(indv_abundance,rownames(lm_res$coefficients)),4],silent = T)
        
        try(return(list(beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        beta_abundance = beta_abundance,
                        se_abundance = se_abundance,
                        p.value_abundance = p.value_abundance,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(beta_abundance    <- c("no lm_res"),silent = T)
        try(se_abundance      <- c("no lm_res"), silent = T)
        try(p.value_abundance <- NA,silent = T)
        
        try(return(list(beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        beta_abundance = beta_abundance,
                        se_abundance = se_abundance,
                        p.value_abundance = p.value_abundance,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_SV))<2){
      
      try(beta_phenotype    <- c("one level"),silent = T)
      try(se_phenotype      <- c("one level"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(beta_abundance    <- c("no lm_res"),silent = T)
      try(se_abundance      <- c("no lm_res"), silent = T)
      try(p.value_abundance <- NA,silent = T)
      
      try(return(list(beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
                      beta_abundance = beta_abundance,
                      se_abundance = se_abundance,
                      p.value_abundance = p.value_abundance,
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
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(SV_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "beta_phenotype","se_phenotype", "p.value_phenotype",
                        "beta_abundance","se_abundance", "p.value_abundance",
                        "N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}


genes_species_correlation<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run){
  # pheno_mat<-pheno_cohort
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      gene_data_i=gene_data[,name,drop=F]
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- genes_species_correlation_2(pheno_mat, gene_data_i, cov_mat_abun, covar)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

genes_species_correlation_2<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_genes = x, X_phenotype = y) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(lm_input[1:5,])
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    lm_input$Y_genes[lm_input$Y_genes==0]=min(lm_input$Y_genes[lm_input$Y_genes>0])/2
    if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){lm_input[,b]=log2(lm_input[,b])}}
    
    if (length(unique(lm_input$X_phenotype))>1&length(unique(lm_input$Y_genes))>1){
      try(cor_test <- cor.test(lm_input[,1], lm_input[,2], method = "spearman"), silent = T)
      if (exists("cor_test")){
        #indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- NA,silent = T)
        try(se_IDage      <- NA, silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- cor_test$estimate,silent = T)
        try(se_phenotype      <- NA, silent = T)
        try(p.value_phenotype <- cor_test$p.value,silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$Y_genes))==1){
      
      try(beta_IDage    <- NA,silent = T)
      try(se_IDage      <- NA, silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("genes all zero"),silent = T)
      try(se_phenotype      <- c("genes all zero"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
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
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}

genes_correct_abundance<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run,gene_species){
  # pheno_mat<-pheno_Variables
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      cat("name",paste(name,"\n"))
      gene_data_i=gene_data[,name,drop=F]
      if (name%in%gene_species$Uniref90){
        #species=gene_species$X[which(gene_species$Y==name)]%>%unique(.) #for file Species_gene_smallest_P
        species=gene_species$organism[which(gene_species$Uniref90==name)]%>%unique(.)#for file Uniref90_species
        cat("species",paste(species,"\n"))
        covar_2<-c(covar,species)
        cat("covar_2",paste(covar_2,"\n"))
      }else{covar_2<-c(covar)}
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- genes_lm_glm_abundance(pheno_mat, gene_data_i, cov_mat_abun, covar_2)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

genes_lm_glm_abundance<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Variables          <-NA
    
    print(covar_use)
    if (length(covar_use)==1){lm_input<-data.frame(Y_genes = x, X_phenotype = y) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}else{lm_input<-data.frame(Y_genes = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(colnames(lm_input))
    Variables<-length(setdiff(colnames(lm_input),c("Y_genes","X_phenotype",covar)))
    print(Variables)
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    
    
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
    
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_genes))<2){
      
      try(beta_IDage    <- c("no lm_res"),silent = T)
      try(se_IDage      <- c("no lm_res"), silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("no lm_res"),silent = T)
      try(se_phenotype      <- c("no lm_res"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Variables=Variables)),
          silent = T)
    }
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(pheno_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Variables")
  
  return(y_x_edge)
}

shortbred_genes_correct_abundance<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run,gene_species){
  # pheno_mat<-pheno_Variables
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      cat("name",paste(name,"\n"))
      gene_data_i=gene_data[,name,drop=F]
      if (name%in%gene_species$GeneName){
        #species=gene_species$X[which(gene_species$Y==name)]%>%unique(.) #for file Species_gene_smallest_P
        species=gene_species$organism[which(gene_species$GeneName==name)]%>%unique(.)#for file Uniref90_species
        cat("species",paste(species,"\n"))
        covar_2<-c(covar,species)
        cat("covar_2",paste(covar_2,"\n"))
      }else{covar_2<-c(covar)}
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- shortbred_genes_lm_glm_abundance(pheno_mat, gene_data_i, cov_mat_abun, covar_2)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

shortbred_genes_lm_glm_abundance<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Variables          <-NA
    
    print(covar_use)
    if (length(covar_use)==1){lm_input<-data.frame(Y_genes = x, X_phenotype = y) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}else{lm_input<-data.frame(Y_genes = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(colnames(lm_input))
    Variables<-length(setdiff(colnames(lm_input),c("Y_genes","X_phenotype",covar)))
    print(Variables)
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    
    
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
    
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Variables=Variables)),
            silent = T)
      }
    }
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_genes))<2){
      
      try(beta_IDage    <- c("no lm_res"),silent = T)
      try(se_IDage      <- c("no lm_res"), silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("no lm_res"),silent = T)
      try(se_phenotype      <- c("no lm_res"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
                      N = N,
                      y_uniq_N     = y_uniq_N,
                      x_uniq_N     = x_uniq_N,
                      y_non_zero_N = y_non_zero_N,
                      x_non_zero_N = x_non_zero_N,
                      y_non_zero_rate = y_non_zero_rate,
                      x_non_zero_rate= x_non_zero_rate,
                      Variables=Variables)),
          silent = T)
    }
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(pheno_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Variables")
  
  return(y_x_edge)
}

my_cooccurrence_original <- function(x1, x2){
  #  x1 <- c(0, 1, 0, 1, 0, 1, 0, 1)
  #  x2 <- c(1, 1, 1, 1, 1, 1, 1, 1)
  return(sum(x1 == x2)/length(x1))
}

my_cooccurrence <- function(x1, x2) {
  valid_idx <- which(!is.na(x1) & !is.na(x2))  # Find indices where both x1 and x2 are NOT NA
  if (length(valid_idx) == 0) {  # If no valid data points, return NA
    return(NA)
  }
  return(sum(x1[valid_idx] == x2[valid_idx]) / length(valid_idx))  # Compute similarity only on valid indices
}


my_coocurrence_mat <- function(m1){
  #  m1<- t(dsv_i)
  
  y_x<-sapply( 
    as.data.frame(m1),
    function(x) Map(function(a,b) my_cooccurrence(a,b),
                    list(x),
                    as.data.frame(m1)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = ncol(m1), byrow = T)
  
  rownames(y_x.unlist) <- colnames(m1)
  colnames(y_x.unlist) <- colnames(m1)
  
  return(y_x.unlist)
}

SV_correct_nospecies<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run){
  # pheno_mat<-pheno_cohort
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      gene_data_i=gene_data[,name,drop=F]
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- SV_correct_nospecies_2(pheno_mat, gene_data_i, cov_mat_abun, covar)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

SV_correct_nospecies_2<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    lm_input<-data.frame(Y_genes = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(lm_input[1:10,])
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
    if (length(unique(lm_input$X_phenotype)) == 2){lm_input$X_phenotype=factor(lm_input$X_phenotype,levels=c("0","1"))}
    if (length(unique(lm_input$X_phenotype)) > 2){
      lm_input$X_phenotype=as.numeric(lm_input$X_phenotype)
      for (a in c("X_phenotype")){lm_input[,a]=qtrans(lm_input[,a])}
    }
    if (length(unique(lm_input$Y_genes)) > 2){for (b in c("Y_genes")){
      lm_input$Y_genes=as.numeric(lm_input$Y_genes)
      lm_input[,b]=qtrans(lm_input[,b])}}
    if (c("Sex")%in%colnames(lm_input)){lm_input$Sex=as.factor(lm_input$Sex)}
    if (length(unique(lm_input$Sex)) < 2){lm_input$Sex=NULL}
    if (c("Age")%in%colnames(lm_input)){for (c in c("Age")){lm_input[,c]=qtrans(lm_input[,c])}}
    
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_genes))<2){
      
      try(beta_IDage    <- c("no lm_res"),silent = T)
      try(se_IDage      <- c("no lm_res"), silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("no lm_res"),silent = T)
      try(se_phenotype      <- c("no lm_res"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
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
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}

shortbred_correct_noSV<-function(gene_data,pheno_mat, SV_mat, cov_mat_abun, covar, genes_to_run){
  # pheno_mat<-pheno_cohort
  # SV_mat<-vsgv_full[sample_name,]
  # cov_mat<-pheno_cov_abundance[sample_name,]
  # covar<-covar
  # abun<-pheno_cov_abundance
  # info<- 
  
  y_x.edge<-NULL
  for (i in 1:length(genes_to_run)){
    # proteins=SVs109__protein
    # i=c("Bacteroides uniformis ATCC 8492:1690_1707")
    # i<-unique(SVs109__protein$SV_Name)[1]
    cat(paste(i,"\n"))
    cat(paste(genes_to_run[i],"\n"))
    if(sum(colnames(gene_data)%in%genes_to_run[i])>0){
      name=colnames(gene_data)[colnames(gene_data)%in%genes_to_run[i]]
      gene_data_i=gene_data[,name,drop=F]
      #print(gene_data_i[1:10,,drop=F])
      y_x_i <- shortbred_correct_noSV_2(pheno_mat, gene_data_i, cov_mat_abun, covar)
      y_x_i$Adjust_Abundance<-"Yes"
      y_x.edge<-rbind(y_x.edge,y_x_i)
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  
  return(y_x.edge)
}

shortbred_correct_noSV_2<-function(pheno_mat,genes_mat,cov_mat_abun,covar_use){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #pheno_mat<- pheno_mat_i #lld_intri[,c(1:3)]
  #genes_mat<- genes_mat #lld_vsv[,c(1:3)]
  #cov_mat_abun<- cov_mat_abun# lld_covar_use
  #covar_use<- covar_use_i #covar_use
  ## test block
  
  my_lm<-function(y,x){
    # Note! Y is from genes_mat_i and X is from pheno_mat. Cause our model is treating SV as Y!
    # y<-pheno_mat[,1]
    #  x<-genes_mat[,1]
    
    beta_IDage    <- NA
    se_IDage      <- NA
    p_IDvalue <- NA
    
    beta_phenotype    <- NA
    se_phenotype      <- NA
    p.value_phenotype <- NA
    
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    Cohort          <-NA
    
    if (length(covar_use)==1){lm_input<-data.frame(Y_genes = x, X_phenotype = y) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}else{lm_input<-data.frame(Y_genes = x, X_phenotype = y,cov_mat_abun[,covar_use]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame}
    #lm_input$Y_genes=lm_input$Y_genes/lm_input$V2
    #lm_input$V2=NULL
    print(lm_input[1:10,])
    # lm_input<-data.frame(Y_genes = genes_LLD[sample_name,c("HMPREF9436_01361")], X_phenotype = bio_aging[sample_name,c("Methyl.Age.Horvath.")],pheno_cov[sample_name,c("Age","Sex","log10_counts")]) %>% sapply(as.numeric) %>% na.omit %>% as.data.frame
    #print(lm_input)
    #### exclude all genes == 0
    #lm_input<-lm_input[lm_input$Y_genes>0,]
    #### change all genes == 0 to smallest non-zero value/10
    #print(lm_input$X_phenotype)
    
    
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y_genes"]))
    x_uniq_N   <- length(unique(lm_input[,"X_phenotype"]))
    
    y_non_zero_N <- sum(!lm_input[,"Y_genes"]==0)
    x_non_zero_N <- sum(!lm_input[,"X_phenotype"]==0)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    Cohort<-cohort
    
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
    
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))>2){
      try(lm_res <- summary(lm(Y_genes~.,data = lm_input)), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        #print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))>2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))==2&length(unique(lm_input$Y_genes))==2){
      try(lm_res <- summary(glm(Y_genes~.,data = lm_input,family = "binomial")), silent = T)
      if (exists("lm_res")){
        indv<-rownames(lm_res$coefficients)[2]
        print(indv)
        
        try(beta_IDage    <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_IDage      <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_IDage <- lm_res$coefficients[match(c("Age"),rownames(lm_res$coefficients)),4],silent = T)
        try(beta_phenotype    <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),1],silent = T)
        try(se_phenotype      <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),2], silent = T)
        try(p.value_phenotype <- lm_res$coefficients[match(c("X_phenotype1"),rownames(lm_res$coefficients)),4],silent = T)
        
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
                        N = N,
                        y_uniq_N     = y_uniq_N,
                        x_uniq_N     = x_uniq_N,
                        y_non_zero_N = y_non_zero_N,
                        x_non_zero_N = x_non_zero_N,
                        y_non_zero_rate = y_non_zero_rate,
                        x_non_zero_rate= x_non_zero_rate,
                        Cohort=cohort)),
            silent = T)
      }else{
        try(beta_IDage    <- c("no lm_res"),silent = T)
        try(se_IDage      <- c("no lm_res"), silent = T)
        try(p.value_IDage <- NA,silent = T)
        try(beta_phenotype    <- c("no lm_res"),silent = T)
        try(se_phenotype      <- c("no lm_res"), silent = T)
        try(p.value_phenotype <- NA,silent = T)
        
        try(return(list(beta_IDage = beta_IDage,
                        se_IDage = se_IDage,
                        p.value_IDage = p.value_IDage,
                        beta_phenotype = beta_phenotype,
                        se_phenotype = se_phenotype,
                        p.value_phenotype = p.value_phenotype,
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
    if (length(unique(lm_input$X_phenotype))<2|length(unique(lm_input$Y_genes))<2){
      
      try(beta_IDage    <- c("no lm_res"),silent = T)
      try(se_IDage      <- c("no lm_res"), silent = T)
      try(p.value_IDage <- NA,silent = T)
      try(beta_phenotype    <- c("no lm_res"),silent = T)
      try(se_phenotype      <- c("no lm_res"), silent = T)
      try(p.value_phenotype <- NA,silent = T)
      
      try(return(list(beta_IDage = beta_IDage,
                      se_IDage = se_IDage,
                      p.value_IDage = p.value_IDage,
                      beta_phenotype = beta_phenotype,
                      se_phenotype = se_phenotype,
                      p.value_phenotype = p.value_phenotype,
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
                    as.data.frame(genes_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 14, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(pheno_mat), byrow = F)
  colnames(y_x.beta)<-colnames(pheno_mat)
  rownames(y_x.beta)<-colnames(genes_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist))
  
  
  colnames(y_x_edge)<-c("X", "Y", "Beta_IDage","SE_IDage", "p_IDage","Beta_phenotype","SE_phenotype", "p_phenotype","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","Cohort")
  
  return(y_x_edge)
}

