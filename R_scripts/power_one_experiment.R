plot_one <- function(method, rep, df_REML, df_HE, m_results, cutoffs, df_GWAS, out_dir){
  pdf(paste(out_dir,"/significant_rep", rep, "_", method, ".pdf", sep=''), width=5, height=5)
  par(mfrow=c(3,1))
  if(sum(df_REML$p_values == 0) > 0){
    df_REML$p_values[which(df_REML$p_values == 0)] <- .Machine$double.xmin
  }
  plot(df_REML$start, -log10(df_REML$p_values), main=paste("REML rep", rep), xlab="position", ylab="-log10(p)")
  abline(v=m_results$pos_causal[rep], col="red") 
  abline(h=cutoffs$cutoff_p_REML, col="gray")
  
  if(sum(df_HE$p_values_HESD_Jackknife == 0) > 0){
    df_HE$p_values_HESD_Jackknife[which(df_HE$p_values_HESD_Jackknife == 0)] <- .Machine$double.xmin
  }
  plot(df_HE$start, -log10(df_HE$p_values_HESD_Jackknife), main=paste("HE SD rep", rep), xlab="position", ylab="-log10(p)")
  abline(v=m_results$pos_causal[rep], col="red") 
  abline(h=cutoffs$cutoff_p_HE_SD, col="gray")
  
  plot(df_GWAS$start, -log10(df_GWAS$p_value), main=paste("GWAS rep", rep), xlab="position", ylab="-log10(p)")
  abline(v=m_results$pos_causal[rep], col="red")
  abline(h=cutoffs$cutoff_p_GWAS, col="gray")
  dev.off()
}

power_one_experiment <- function(hsquared, REPS, folder, tree_type){
  #which power do I have with heritability 0.005, 1000 diploids
  expected_power_GWAS <- pwr.norm.test(d=hsquared^2, n=1000, sig.level = 0.05, power=NULL)$power
  
  out_dir = paste(folder,"/h", hsquared, sep='')
  
  position_interval <- c(49350000, 49650000)
  cutoffs <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/eGRM/", tree_type, "/p_value_cutoffs.csv", sep=''))
  
  reps <- REPS
  
  m_results <- matrix(nrow=reps, ncol=12)
  colnames(m_results) <- c("REML", "HE_SD", "HE_CP", "GWAS", "REML_equal0.5", "pos_min_REML", "pos_min_HESD", "pos_min_HECP", "pos_min_GWAS", "pos_causal", "af_causal", "h2_empiric")
  m_results <- data.frame(m_results)
  

  
  for(rep in 1:reps){
    # read REML results
    df_REML <- read.csv(paste(out_dir,"/rep", rep, "/power_sims_", rep, "_trees_REML_results.csv", sep=''))
    df_REML <- df_REML[-1,]
    df_REML <- df_REML[which(df_REML$start > position_interval[1] & df_REML$start < position_interval[2]),]
    # df_REML <- df_REML[-nrow(df_REML),]
    df_REML$p_values[df_REML$p_values == 0] <- .Machine$double.xmin
    
    # read HE results
    df_HE <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_trees_HE_results.csv", sep=''))
    df_HE <- df_HE[-1,]
    df_HE <- df_HE[which(df_HE$start > position_interval[1] & df_HE$start < position_interval[2]),]
    # df_HE <- df_HE[-nrow(df_HE),]
    if(sum(df_HE$p_values_HESD_Jackknife == 0) > 0){
      df_HE$p_values_HESD_Jackknife[df_HE$p_values_HESD_Jackknife == 0] <- .Machine$double.xmin
    }
    if(sum(df_HE$p_values_HECP_Jackknife == 0) > 0){
      df_HE$p_values_HECP_Jackknife[df_HE$p_values_HECP_Jackknife == 0] <- .Machine$double.xmin
    }
    
    #read GWAS results
    df_GWAS <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_variants_results.csv", sep=''))
    df_GWAS <- df_GWAS[which(df_GWAS$start > position_interval[1] & df_GWAS$start < position_interval[2]),]
    
    #read phenotype spects
    df_pheno <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_pheno_causal_vars.csv", sep=''))
    m_results$pos_causal[rep] <- df_pheno$start[df_pheno$causal == TRUE]
    m_results$af_causal[rep] <- df_pheno$allele_freq[df_pheno$causal == TRUE]
    m_results$h2_empiric[rep] <- df_pheno$var_genotypic_empiric[1] / df_pheno$var_phenotypic[1]
    
    m_results$REML[rep] <- -log10(min(df_REML$p_values))
    pos_with_min_p <- df_REML$start[which(df_REML$p_values == min(df_REML$p_values))]
    distances <- abs(pos_with_min_p - m_results$pos_causal[rep])
    m_results$pos_min_REML[rep] <- pos_with_min_p[distances == min(distances)]
    
    m_results$HE_SD[rep] <- -log10(min(df_HE$p_values_HESD_Jackknife))
    pos_with_min_p <- df_HE$start[which(df_HE$p_values_HESD_Jackknife == min(df_HE$p_values_HESD_Jackknife))]
    distances <- abs(pos_with_min_p - m_results$pos_causal[rep])
    m_results$pos_min_HESD[rep] <- pos_with_min_p[distances == min(distances)]
    
    m_results$HE_CP[rep] <- -log10(min(df_HE$p_values_HECP_Jackknife))
    pos_with_min_p <- df_HE$start[which(df_HE$p_values_HECP_Jackknife == min(df_HE$p_values_HECP_Jackknife))]
    distances <- abs(pos_with_min_p - m_results$pos_causal[rep])
    m_results$pos_min_HECP[rep] <- pos_with_min_p[distances == min(distances)]
    
    m_results$GWAS[rep] <- -log10(min(df_GWAS$p_value))
    pos_with_min_p <- df_GWAS$start[which(df_GWAS$p_value == min(df_GWAS$p_value))]
    distances <- abs(pos_with_min_p - m_results$pos_causal[rep])
    m_results$pos_min_GWAS[rep] <- pos_with_min_p[distances == min(distances)]
    
    if(m_results$REML[rep] > cutoffs$cutoff_p_REML){
      plot_one("REML", rep, df_REML, df_HE, m_results, cutoffs, df_GWAS, out_dir)
    }
    
    if(m_results$GWAS[rep] > cutoffs$cutoff_p_GWAS){
      plot_one("GWAS", rep, df_REML, df_HE, m_results, cutoffs, df_GWAS, out_dir)
    }
    
    if(m_results$HE_SD[rep] > cutoffs$cutoff_p_HE_SD){
      plot_one("HE_SD", rep, df_REML, df_HE, m_results, cutoffs, df_GWAS, out_dir)
    }
    
    if(m_results$REML[rep] > cutoffs$cutoff_p_REML & abs(m_results$pos_min_REML[rep] - m_results$pos_causal[rep]) < 5000){
      # print(paste("within region:", abs(m_results$pos_min_REML[rep] - m_results$pos_causal[rep])))
      plot_one("REML_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
    }
    
    if(m_results$GWAS[rep] > cutoffs$cutoff_p_GWAS & abs(m_results$pos_min_GWAS[rep] - m_results$pos_causal[rep]) < 5000){
      plot_one("GWAS_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
    }
    
    if(m_results$HE_SD[rep] > cutoffs$cutoff_p_HE_SD & abs(m_results$pos_min_HESD[rep] - m_results$pos_causal[rep]) < 5000){
      plot_one("HE_SD_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
    }
  }
  
  #--------------------------
  # power
  #--------------------------
  power_REML <- sum(m_results[,"REML"] > cutoffs$cutoff_p_REML) / reps
  power_GWAS <- sum(m_results[,"GWAS"] > cutoffs$cutoff_p_GWAS) / reps
  power_HE_SD <- sum(m_results[,"HE_SD"] > cutoffs$cutoff_p_HE_SD) / reps
  power_HE_CP <- sum(m_results[,"HE_CP"] > cutoffs$cutoff_p_HE_CP) / reps
  
  power_REML_region <- sum(m_results[,"REML"] > cutoffs$cutoff_p_REML & abs(m_results[,"pos_min_REML"] - m_results[, "pos_causal"]) < 5000) / reps
  power_GWAS_region <- sum(m_results[,"GWAS"] > cutoffs$cutoff_p_GWAS & abs(m_results[,"pos_min_GWAS"] - m_results[, "pos_causal"]) < 5000) / reps
  power_HE_SD_region <- sum(m_results[,"HE_SD"] > cutoffs$cutoff_p_HE_SD & abs(m_results[,"pos_min_HESD"] - m_results[, "pos_causal"]) < 5000) / reps
  power_HE_CP_region <- sum(m_results[,"HE_CP"] > cutoffs$cutoff_p_HE_CP & abs(m_results[,"pos_min_HECP"] - m_results[, "pos_causal"]) < 5000) / reps
  
  out_t <- cbind(power_REML, power_GWAS, power_HE_SD, power_HE_CP, expected_power_GWAS, power_REML_region, power_GWAS_region, power_HE_SD_region, power_HE_CP_region)
  write.table(out_t, file=paste(out_dir, "/power_results.txt", sep=''), quote=FALSE, row.names=FALSE)
  
  #--------------------------
  # location of smallest p-value
  #--------------------------
  plot(m_results$pos_min_REML, m_results$REML, xlab="location smallest p-value", ylab="smallest p-value")
  
  #--------------------------
  # emiric heritability
  #--------------------------
  pdf(paste(out_dir,"/hist_heritability.pdf", sep=''), width=4, height=4)
  hist(m_results$h2_empiric, xlab="empiric heritability")
  dev.off()
  
  #--------------------------
  # causal allele freq
  #--------------------------
  pdf(paste(out_dir,"/hist_causal_allele_freq.pdf", sep=''), width=4, height=4)
  hist(m_results$af_causal, xlab="allele freq")
  dev.off()
  
  #--------------------------
  # hist of indeces of minimums
  #--------------------------
  pdf(paste(out_dir,"/hist_lowest_pvalue_index.pdf", sep=''), width=8, height=8)
  par(mfrow=c(2,2))
  hist(m_results$pos_min_REML, xlab ="index of min p-value in ARG REML", breaks=20, main="")
  hist(m_results$pos_min_GWAS, xlab ="index of min p-value of all variants GWAS", breaks=20, main="")
  hist(m_results$pos_min_HESD, xlab ="index of min p-value in ARG HE SD", breaks=20, main="")
  hist(m_results$pos_min_HECP, xlab ="index of min p-value in ARG HE CP", breaks=20, main="")
  dev.off()
  
  #--------------------------
  # cutoff plot
  #--------------------------
  
  # pdf(paste(out_dir,"/p_values_cutoff_separate.pdf", sep=''), width=10, height=5)
  # 
  # par(mfrow=c(1,3))
  # par(pty="s")
  # plot(sort(m_results[,"REML"], decreasing=T), col="dodgerblue", yaxt='n', ylim=c(1,6.5), xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
  # lines(sort(m_results[,"GWAS"], decreasing=T), col="black")
  # axis(side=2, las=2)
  # legend(legend=c("GWAS", "REML"), x="topright", ncol=1, col=c("black", "dodgerblue"), lty=1, bty='n')
  # 
  # plot(sort(m_results[,"HE_SD"], decreasing=T), col="red", yaxt='n', xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
  # lines(sort(m_results[,"HE_CP"], decreasing=T), col="orange2")
  # axis(side=2, las=2)
  # legend(legend=c("HE_SD","HE_CP"), x="topright", ncol=1, col=c("red", "orange2"), lty=1, bty='n')
  # dev.off()
  # 
  # pdf(paste(out_dir,"/p_values_cutoff_all.pdf", sep=''), width=5, height=5)
  # 
  # par(pty="s")
  # plot(sort(m_results[,"REML"], decreasing=T), col="dodgerblue", yaxt='n', ylim=c(0.5,6), xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
  # lines(sort(m_results[,"GWAS"], decreasing=T), col="black")
  # lines(sort(m_results[,"HE_SD"], decreasing=T), col="red")
  # lines(sort(m_results[,"HE_CP"], decreasing=T), col="orange2")
  # axis(side=2, las=2)
  # legend(legend=c("GWAS", "REML","HE_SD","HE_CP"), x="topright", ncol=1, col=c("black", "dodgerblue", "red", "orange2"), lty=1, bty='n')
  # 
  # dev.off()
}
