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
      
      plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE CP - Allele freq:", m[i,1], "Beta:", beta), xlab="position", ylab="-log10(p)")
      abline(v=GWAS$start[which(GWAS$causal == TRUE)], col="red")
      
      ARGWAS_p <- ARGWAS$p_values_HESD_Jackknife
      ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
      GWAS_p <- GWAS$p_value
      GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
      max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
      
      plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE SD - Allele freq:", m[i,1], "Beta:", beta), xlab="position", ylab="-log10(p)")
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

#------------------------------------
# allelic heterogeneity
#------------------------------------


dir <- "/data/ARGWAS/experiments_N500/aH/oneTree/HE/"

for(af in c("")){#c(0.4, 0.2, 0.1, 0.05)
  for(noise in c("", "_withNoise", "_withLotsNoise")){
    for(beta in c("") ){ #c(1,0.01)
      for(tree in c("")){
        pdf(paste(dir,"test_2_results_beta", beta, "_freq", af, tree, noise, ".pdf", sep=''), width=10, height=13)
        
        par(mfrow=c(7,3))
        
        # m_layout <- matrix(nrow=3, ncol=8)
        # m_layout[,1] <- c(1:3)
        # m_layout[,2] <- c(4:6)
        # m_layout[,3] <- c(7:9)
        # m_layout[,4] <- c(10:12)
        # m_layout[,5] <- c(13:15)
        # m_layout[,6] <- c(16:18)
        # m_layout[,7] <- c(19:21)
        # m_layout[,8] <- c(22:24)
        # 
        # layout(m_layout)

        for(prop in c(1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){ #1, 0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
          ARGWAS <- read.csv(paste(dir,"freq", af, "_indexUntyped_aH_beta", beta, "_propTyped", prop , tree, noise, "_trees_results.csv", sep=''))
          ARGWAS <- ARGWAS[-1,]
          GWAS <- read.csv(paste(dir,"freq", af, "_indexUntyped_aH_beta", beta, "_propTyped", prop, tree, noise, "_variants_results.csv", sep=''))
          pt_info <- read.csv(paste(dir,"freq", af, "_indexUntyped_aH_beta", beta, "_propTyped", prop, tree, noise, "_pheno_causal_vars.csv", sep = ''))
          causal_pos_neg <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas < 0)]
          causal_pos_pos <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas > 0)]
          
          ARGWAS_p <- ARGWAS$p_values_HECP_Jackknife
          ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
          GWAS_p <- GWAS$p_value
          GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
          max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
          
          plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE CP - PrTyped:", prop), xlab="position", ylab="-log10(p)")
          abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
          abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
          
          ARGWAS_p <- ARGWAS$p_values_HESD_Jackknife
          ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
          GWAS_p <- GWAS$p_value
          GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
          max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
          
          plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE SD - PrTyped:", prop), xlab="position", ylab="-log10(p)")
          abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
          abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
          
          plot(GWAS$end, -log10(GWAS_p), main=paste("GWAS - PrTyped:", prop), xlab="position", ylab="-log10(p)")
          abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
          abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
          abline(h=8, col="gray")
          
          # if(ARGWAS$start[which(ARGWAS$causal == TRUE)] != GWAS$start[which(GWAS$causal == TRUE)]){
          #   stop("causal positions of ARGWAS and GWAS don't match")
        }
        dev.off()
      }
    }
  }
}


#---------------
# true tree
#---------------

for(noise in c("", "_withNoise", "_withLotsNoise")){
  
  pdf(paste(dir,"test_2_results", noise, "_trueTree.pdf", sep=''), width=6, height=8)
  par(mfrow=c(4,2))
  
  for(prop in c(1)){ #1, 0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
    ARGWAS <- read.csv(paste(dir,"freq_indexUntyped_aH_beta_propTyped", prop , "_trueTree", noise, "_trees_results.csv", sep=''))
    ARGWAS <- ARGWAS[-1,]
    
    ARGWAS_p <- ARGWAS$p_values_HECP_Jackknife
    ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
    max_q <- max(-log10(ARGWAS_p))
    
    pt_info <- read.csv(paste(dir,"freq_indexUntyped_aH_beta_propTyped", prop, "_trueTree", noise, "_pheno_causal_vars.csv", sep = ''))
    causal_pos_neg <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas < 0)]
    causal_pos_pos <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas > 0)]
    
    plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE CP - PrTyped:", prop), xlab="position", ylab="-log10(p)")
    abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
    abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
    
    ARGWAS_p <- ARGWAS$p_values_HESD_Jackknife
    ARGWAS_p[which(ARGWAS_p == 0)] <- .Machine$double.xmin
    max_q <- max(-log10(ARGWAS_p))
    
    plot(ARGWAS$start, -log10(ARGWAS_p), main=paste("HE SD - PrTyped:", prop), xlab="position", ylab="-log10(p)")
    abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
    abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
  }
  
  for(prop in c(1, 0.5, 0.2, 0.1, 0.05, 0.02)){ #1, 0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
    GWAS <- read.csv(paste(dir,"freq_indexUntyped_aH_beta_propTyped", prop, "_trueTree", noise, "_variants_results.csv", sep=''))
    pt_info <- read.csv(paste(dir,"freq_indexUntyped_aH_beta_propTyped", prop, "_trueTree", noise, "_pheno_causal_vars.csv", sep = ''))
    causal_pos_neg <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas < 0)]
    causal_pos_pos <- pt_info$start[which(pt_info$causal == TRUE & pt_info$betas > 0)]
    
    GWAS_p <- GWAS$p_value
    GWAS_p[which(GWAS_p == 0)] <- .Machine$double.xmin
    max_q <- max(-log10(ARGWAS_p), -log10(GWAS_p))
    
    plot(GWAS$end, -log10(GWAS_p), main=paste("GWAS - PrTyped:", prop), xlab="position", ylab="-log10(p)")
    abline(v=causal_pos_pos, col = adjustcolor("red", alpha = 0.05))
    abline(v=causal_pos_neg, col = adjustcolor("blue", alpha = 0.05))
    abline(h=8, col="gray")
    
  }
  
  dev.off()
}

# freq <- read.table("/home/vivian/postdoc_USC/AIM/experiments_N5k/simulation5k_simulated_sample_variants.csv", sep=',', header=TRUE)
# hist(freq$allele_freq)
# 
# freq_0.01 <- read.table("/home/vivian/postdoc_USC/AIM/experiments_N5k/relate_simulation5k_propTyped0.01/simulation5k_propTyped0.01_filtered_sample_variants.csv", sep=',', header=TRUE)
# hist(freq_0.01$allele_freq[freq_0.01$typed == "False" & freq_0.01$position >= 49461796 & freq_0.01$position <= 49476796], breaks=seq(0,0.5,0.001))
# 
#-----------------------------
# downsampled RELATE
#-----------------------------


dir <- "/home/vivian/postdoc_USC/AIM/experiments/"

for(af in c(0.4)){#c(0.4, 0.2, 0.1, 0.05)
  for(noise in c("", "_withNoise")){
    for(beta in c(1) ){ #c(1,0.01)
      pdf(paste(dir,"test_2_results_beta", beta, "_freq", af, noise, ".pdf", sep=''), width=20, height=10)
      
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
      # /data/ARGWAS/experiments/freq0.4_indexUntyped_beta1_propTyped0.7_trees_results.csv
      for(prop in c(0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01)){ #1, 0.7, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01
        ARGWAS <- read.csv(paste(dir,"freq", af, "_indexUntyped_beta", beta, "_propTyped", prop , noise, "_trees_results.csv", sep=''))
        ARGWAS <- ARGWAS[-1,]
        GWAS <- read.csv(paste(dir,"freq", af, "_indexUntyped_beta", beta, "_propTyped", prop, noise, "_variants_results.csv", sep=''))
        pt_info <- read.csv(paste(dir,"freq0.4_indexUntyped_beta1_propTyped", prop, noise, "_pheno_causal_vars.csv", sep = ''))
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
        
        plot(GWAS$end, -log10(GWAS_p), main=paste("GWAS - AFreq:", af, "Beta:", beta, "PrTyped:", prop), xlab="position", ylab="-log10(p)")
        abline(v=causal_pos, col="red")
        # if(ARGWAS$start[which(ARGWAS$causal == TRUE)] != GWAS$start[which(GWAS$causal == TRUE)]){
        #   stop("causal positions of ARGWAS and GWAS don't match")
      }
    }
    dev.off()
  }
}

t <- read.table("/home/vivian/postdoc_USC/AIM/experiments/freq0.2_indexUntyped_beta1_propTyped0.5_withNoise_pheno_causal_vars.csv", sep=',', header=TRUE)
sum(t$causal)
