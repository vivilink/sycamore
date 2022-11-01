plot_one <- function(method, rep, df_REML, df_HE, cutoff_REML, cutoff_HE, cutoff_GWAS, df_GWAS, out_dir, df_pheno){
  causal_pos_neg <- df_pheno$start[which(df_pheno$causal == TRUE & df_pheno$betas < 0)]
  causal_pos_pos <- df_pheno$start[which(df_pheno$causal == TRUE & df_pheno$betas > 0)]
  
  pdf(paste(out_dir,"/significant_rep", rep, "_", method, ".pdf", sep=''), width=6, height=5)
  par(mfrow=c(3,1))
  
  #plot REML
  if(sum(df_REML$p_values == 0, na.rm = TRUE) > 0){
    df_REML$p_values[which(df_REML$p_values == 0)] <- .Machine$double.xmin
  }
  plot(df_REML$start, -log10(df_REML$p_values), main=paste("REML rep", rep), xlab="position", ylab="-log10(p)", ylim=c(0,10))
  abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.2))
  abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.2))
  abline(h=cutoff_REML, col="gray")
  
  #plot HE
  if(sum(df_HE$p_values_HESD_Jackknife == 0, na.rm = TRUE) > 0){
    df_HE$p_values_HESD_Jackknife[which(df_HE$p_values_HESD_Jackknife == 0)] <- .Machine$double.xmin
  }
  plot(df_HE$start, -log10(df_HE$p_values_HESD_Jackknife), main=paste("HE SD rep", rep), xlab="position", ylab="-log10(p)")
  abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.2))
  abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.2))
  abline(h=cutoff_HE, col="gray")

  #plot GWAS
  plot(df_GWAS$start, -log10(df_GWAS$p_value), main=paste("GWAS rep", rep), xlab="position", ylab="-log10(p)", ylim=c(0,10))
  abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.2))
  abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.2))
  abline(h=cutoff_GWAS, col="gray")
  dev.off()
}

power_one_experiment <- function(hsquared, REPS, folder, tree_type, region_type){
  #which power do I have with heritability 0.005, 1000 diploids
  expected_power_GWAS <- pwr.norm.test(d=hsquared^2, n=1000, sig.level = 0.05, power=NULL)$power
  
  out_dir = paste(folder,"/h", hsquared, sep='')
  
  position_interval <- c(49000000, 49988000)
  # position_interval <- c(49350000, 49650000)
  

  reps <- REPS
  
  #prepare result data frames
  m_results_REML <- matrix(nrow=reps, ncol=10)
  colnames(m_results_REML) <- c("REML", "HE_SD", "HE_CP", "REML_equal0.5", "pos_min_REML", "pos_min_HESD", "pos_min_HECP", "pos_causal", "af_causal", "h2_empiric")
  m_results_REML <- data.frame(m_results_REML)

  m_results_GWAS <- matrix(nrow=reps, ncol=2)
  colnames(m_results_GWAS) <- c("GWAS", "pos_min_GWAS")
  m_results_GWAS <- data.frame(m_results_GWAS)
  cutoff_GWAS <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/", region_type, "/p_value_cutoffs_", "GWAS", ".csv", sep=''))$x
  
  #covariance_types 
  covariance_types <- c("eGRM", "GRM")
  result_matrices <- list()
  cutoff_files <- list()
  for(i in 1:length(covariance_types)){
    result_matrices[[covariance_types[i]]] <- m_results_REML
    
    cutoff_files[[covariance_types[i]]] <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type,  "/", region_type, "/p_value_cutoffs_", covariance_types[i], ".csv", sep=''))

  }
  result_matrices[["GWAS"]] <- m_results_GWAS
  
  for(rep in 1:reps){
    print(rep)
    #read phenotype specs
    df_pheno <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_pheno_causal_vars.csv", sep=''))
    
    #read GWAS results
    df_GWAS <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_GWAS_variants_results.csv", sep=''))
    df_GWAS <- df_GWAS[which(df_GWAS$start > position_interval[1] & df_GWAS$start < position_interval[2]),]
    
    result_matrices[["GWAS"]]$GWAS[rep] <- -log10(min(df_GWAS$p_value))
    pos_with_min_p <- df_GWAS$start[which(df_GWAS$p_value == min(df_GWAS$p_value))]
    distances <- abs(pos_with_min_p - result_matrices[[i]]$pos_causal[rep])
    # result_matrices[[i]]$pos_min_GWAS[rep] <- pos_with_min_p[distances == min(distances)]
    
    
    for(i in 1:length(covariance_types)){
      # read REML results
      df_REML <- read.csv(paste(out_dir,"/rep", rep, "/power_sims_", rep, "_", covariance_types[[i]], "_trees_REML_results.csv", sep=''))
      df_REML <- df_REML[-1,]
      df_REML <- df_REML[-nrow(df_REML),]
      df_REML <- df_REML[which(df_REML$start > position_interval[1] & df_REML$start < position_interval[2]),]
      df_REML$p_values[df_REML$p_values == 0] <- .Machine$double.xmin
      
      result_matrices[[i]]$REML[rep] <- -log10(min(df_REML$p_values, na.rm = TRUE))
      pos_with_min_p <- df_REML$start[which(df_REML$p_values == min(df_REML$p_values, na.rm = TRUE))]
      distances <- abs(pos_with_min_p - result_matrices[[i]]$pos_causal[rep])
      # result_matrices[[i]]$pos_min_REML[rep] <- pos_with_min_p[distances == min(distances)]
      
      # read HE results
      df_HE <- read.csv(paste(out_dir,"/rep", rep,"/power_sims_", rep, "_", covariance_types[[i]], "_trees_HE_results.csv", sep=''))
      df_HE <- df_HE[-1,]
      df_HE <- df_HE[which(df_HE$start > position_interval[1] & df_HE$start < position_interval[2]),]
      # df_HE <- df_HE[-nrow(df_HE),]
      if(sum(df_HE$p_values_HESD_Jackknife == 0, na.rm = TRUE) > 0){
        df_HE$p_values_HESD_Jackknife[df_HE$p_values_HESD_Jackknife == 0] <- .Machine$double.xmin
      }
      if(sum(df_HE$p_values_HECP_Jackknife == 0, na.rm = TRUE) > 0){
        df_HE$p_values_HECP_Jackknife[df_HE$p_values_HECP_Jackknife == 0] <- .Machine$double.xmin
      }

      result_matrices[[i]]$HE_SD[rep] <- -log10(min(df_HE$p_values_HESD_Jackknife, na.rm = TRUE))
      pos_with_min_p <- df_HE$start[which(df_HE$p_values_HESD_Jackknife == min(df_HE$p_values_HESD_Jackknife))]
      distances <- abs(pos_with_min_p - result_matrices[[i]]$pos_causal[rep])
      # result_matrices[[m]]$pos_min_HESD[rep] <- pos_with_min_p[distances == min(distances)]

      result_matrices[[i]]$HE_CP[rep] <- -log10(min(df_HE$p_values_HECP_Jackknife, na.rm = TRUE))
      pos_with_min_p <- df_HE$start[which(df_HE$p_values_HECP_Jackknife == min(df_HE$p_values_HECP_Jackknife))]
      distances <- abs(pos_with_min_p - result_matrices[[i]]$pos_causal[rep])
      # result_matrices[[m]]$pos_min_HECP[rep] <- pos_with_min_p[distances == min(distances)]
      
      
      #plot significant replicates
      if(result_matrices[[i]]$REML[rep] > cutoff_files[[covariance_types[i]]]$cutoff_p_REML){
        plot_one(method=paste("REML", covariance_types[i], sep='_'), rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
      }

      if(result_matrices[[i]]$HE_SD[rep] > cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD){
        plot_one(method=paste("HE_SD", covariance_types[i], sep='_'), rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
      }
      
      if(result_matrices[["GWAS"]]$GWAS[rep] > cutoff_GWAS){
        plot_one(method="GWAS", rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
      }
      
      if(rep %in% c(23)){
        plot_one(method=paste("REML", covariance_types[i], sep='_'), rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
        plot_one(method=paste("HE_SD", covariance_types[i], sep='_'), rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
        plot_one(method="GWAS", rep=rep, df_REML=df_REML, df_HE=df_HE, cutoff_REML=cutoff_files[[covariance_types[i]]]$cutoff_p_REML, cutoff_HE=cutoff_files[[covariance_types[i]]]$cutoff_p_HE_SD, cutoff_GWAS=cutoff_GWAS, df_GWAS=df_GWAS, out_dir=out_dir, df_pheno=df_pheno)
      }
    }
  }
  

  
  # if(m_results$REML[rep] > cutoffs$cutoff_p_REML & abs(m_results$pos_min_REML[rep] - m_results$pos_causal[rep]) < 5000){
  #   # print(paste("within region:", abs(m_results$pos_min_REML[rep] - m_results$pos_causal[rep])))
  #   plot_one("REML_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
  # }
  # 
  # if(m_results$GWAS[rep] > cutoffs$cutoff_p_GWAS & abs(m_results$pos_min_GWAS[rep] - m_results$pos_causal[rep]) < 5000){
  #   plot_one("GWAS_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
  # }
  # 
  # if(m_results$HE_SD[rep] > cutoffs$cutoff_p_HE_SD & abs(m_results$pos_min_HESD[rep] - m_results$pos_causal[rep]) < 5000){
  #   plot_one("HE_SD_region", rep=rep, df_REML=df_REML, m_results=m_results, cutoffs=cutoffs, df_HE=df_HE, df_GWAS=df_GWAS, out_dir=out_dir)
  # }
  

    
  
  #--------------------------
  # power
  #--------------------------
  power_REML_eGRM <- sum(result_matrices[["eGRM"]][,"REML"] > cutoff_files[["eGRM"]]$cutoff_p_REML) / reps
  power_REML_GRM <- sum(result_matrices[["GRM"]][,"REML"] > cutoff_files[["GRM"]]$cutoff_p_REML) / reps
  power_GWAS <- sum(result_matrices[["GWAS"]][,"GWAS"] > cutoff_GWAS) / reps

  
  # power_REML_region <- sum(m_results[,"REML"] > cutoffs$cutoff_p_REML & abs(m_results[,"pos_min_REML"] - m_results[, "pos_causal"]) < 5000) / reps
  # power_GWAS_region <- sum(m_results[,"GWAS"] > cutoffs$cutoff_p_GWAS & abs(m_results[,"pos_min_GWAS"] - m_results[, "pos_causal"]) < 5000) / reps
  # power_HE_SD_region <- sum(m_results[,"HE_SD"] > cutoffs$cutoff_p_HE_SD & abs(m_results[,"pos_min_HESD"] - m_results[, "pos_causal"]) < 5000) / reps
  # power_HE_CP_region <- sum(m_results[,"HE_CP"] > cutoffs$cutoff_p_HE_CP & abs(m_results[,"pos_min_HECP"] - m_results[, "pos_causal"]) < 5000) / reps
  
  out_t <- cbind(power_REML_eGRM, power_REML_GRM, power_GWAS, expected_power_GWAS)#, power_REML_region, power_GWAS_region, power_HE_SD_region, power_HE_CP_region)
  write.table(out_t, file=paste(out_dir, "/power_results.txt", sep=''), quote=FALSE, row.names=FALSE)
  
}
