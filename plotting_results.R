m <- matrix(nrow=4, ncol=2)
colnames(m) <- c("freq", "index")
m[1,] <- c(0.4, 1216)
m[2,] <- c(0.2, 1423)
m[3,] <- c(0.102, 1219)
m[4,] <- c(0.05, 1217)

for(noise in c("_withNoise", "")){
  for(beta in c(1,0.01)){
    pdf(paste("/data/ARGWAS/experiments/test_2_results_beta", beta, noise, ".pdf", sep=''), width=20, height=10)
    m_layout <- matrix(nrow=3, ncol=4)
    m_layout[,1] <- c(1:3)
    m_layout[,2] <- c(4:6)
    m_layout[,3] <- c(7:9)
    m_layout[,4] <- c(10:12)
    
    layout(m_layout)
    for(i in 1:4){
      ARGWAS <- read.csv(paste("/data/ARGWAS/experiments/freq", m[i,1], "_index", m[i,2] - 1, "_beta", beta, noise, "_trees_results.csv", sep=''))
      ARGWAS <- ARGWAS[-1,]
      GWAS <- read.csv(paste("/data/ARGWAS/experiments/freq", m[i,1], "_index", m[i,2] - 1, "_beta", beta, noise, "_variants.csv", sep=''))
      
      ARGWAS_p <- ARGWAS$p_values_HECP_Jackknife
      ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
      GWAS_p <- GWAS$p_value
      GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
      max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
      
      plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("ARGWAS CP - Allele freq:", m[i,1], "Beta:", beta), xlab="position", ylab="-log10(p)")
      abline(v=GWAS$start[which(GWAS$causal == TRUE)], col="red")
      
      ARGWAS_p <- ARGWAS$p_values_HESD_Jackknife
      ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
      GWAS_p <- GWAS$p_value
      GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
      max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
      
      plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("ARGWAS SD - Allele freq:", m[i,1], "Beta:", beta), xlab="position", ylab="-log10(p)")
      abline(v=GWAS$start[which(GWAS$causal == TRUE)], col="red")
      
      plot(GWAS$start, -log10(GWAS_p), main=paste("GWAS - Allele freq:", m[i,1], "Beta:", beta), xlab="position", ylab="-log10(p)")
      abline(v=GWAS$start[which(GWAS$causal == TRUE)], col="red")
      # if(ARGWAS$start[which(ARGWAS$causal == TRUE)] != GWAS$start[which(GWAS$causal == TRUE)]){
      #   stop("causal positions of ARGWAS and GWAS don't match")
      # }
    }
    dev.off()
  }
}

for(af in c(0.4)){#c(0.4, 0.2, 0.1, 0.05)
  for(noise in c("", "_withNoise")){
    for(beta in c(1) ){ #c(1,0.01)
      pdf(paste("/data/ARGWAS/experiments/test_2_results_beta", beta, "_freq", af, noise, ".pdf", sep=''), width=20, height=10)
      
      m_layout <- matrix(nrow=3, ncol=8)
      m_layout[,1] <- c(1:3)
      m_layout[,2] <- c(4:6)
      m_layout[,3] <- c(7:9)
      m_layout[,4] <- c(10:12)
      m_layout[,5] <- c(13:15)
      m_layout[,6] <- c(16:18)
      m_layout[,7] <- c(19:21)
      m_layout[,8] <- c(22:24)
      
      layout(m_layout)
      /data/ARGWAS/experiments/freq0.4_indexUntyped_beta1_propTyped0.7_trees_results.csv
      for(prop in c(0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){ #1, 0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
        ARGWAS <- read.csv(paste("/data/ARGWAS/experiments/freq", af, "_indexUntyped_beta", beta, "_propTyped", prop , noise, "_trees_results.csv", sep=''))
        ARGWAS <- ARGWAS[-1,]
        GWAS <- read.csv(paste("/data/ARGWAS/experiments/freq", af, "_indexUntyped_beta", beta, "_propTyped", prop, noise, "_variants_results.csv", sep=''))
        pt_info <- read.csv(paste("/data/ARGWAS/experiments/freq0.4_indexUntyped_beta1_propTyped", prop, noise, "_pheno_causal_vars.csv", sep = ''))
        causal_pos <- pt_info$start[which(pt_info$causal == TRUE)]
        
        #get actual allele freq simulated
        # af <- read.csv()
        
        ARGWAS_p <- ARGWAS$p_values_HECP_Jackknife
        ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
        GWAS_p <- GWAS$p_value
        GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
        max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
        
        plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("ARGWAS CP - AFreq:", af, "Beta:", beta, "PrTyped:", prop), xlab="position", ylab="-log10(p)")
        abline(v=causal_pos, col="red")
        
        ARGWAS_p <- ARGWAS$p_values_HESD_Jackknife
        ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
        GWAS_p <- GWAS$p_value
        GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
        max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
        
        plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("ARGWAS SD - AFreq:", af, "Beta:", beta, "PrTyped:", prop), xlab="position", ylab="-log10(p)")
        abline(v=causal_pos, col="red")
        
        plot(GWAS$start, -log10(GWAS_p), main=paste("GWAS - AFreq:", af, "Beta:", beta, "PrTyped:", prop), xlab="position", ylab="-log10(p)")
        abline(v=GWAS$start[which(GWAS$causal == TRUE)], col="red")
        # if(ARGWAS$start[which(ARGWAS$causal == TRUE)] != GWAS$start[which(GWAS$causal == TRUE)]){
        #   stop("causal positions of ARGWAS and GWAS don't match")
      }
    }
    dev.off()
  }
}

