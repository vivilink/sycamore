library(ACAT)

setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k")

reps <- 300
cutoff_rep <- 0.05 * reps
position_interval <- c(49000000, 50000000)
run_acat <- TRUE
# position_interval <- c(0, 1000000)
# result_file_prefix <- "two_pops"
result_file_prefix <- "cutoff_sims"

covariance_types <- c("eGRM", "GRM")
m_results_VC <- list()
m_results_GWAS <- matrix(nrow=reps, ncol=3)
colnames(m_results_GWAS) <- c( "GWAS", "index_min_GWAS", "p_value_10")
m_results_GWAS <- data.frame(m_results_GWAS)
if(run_acat){
  m_results_acat <- matrix(nrow=reps, ncol=3)
  colnames(m_results_acat) <- c( "acat", "index_min_acat", "p_value_10")
  m_results_acat <- data.frame(m_results_acat)
}
tmp <- matrix(nrow=reps, ncol=7)
colnames(tmp) <- c("REML", "HE_SD", "HE_CP","REML_equal0.5", "index_min_REML", "index_min_HESD", "index_min_HECP")

m_results_VC[[1]] <- data.frame(tmp)
m_results_VC[[2]] <- data.frame(tmp)


for(i in 1:length(covariance_types)){
  m_p_values_REML <- matrix(nrow=reps, ncol=4)
  colnames(m_p_values_REML) <- c("p_value_10", "p_value_100", "p_value_200", "p_value_400")
  m_p_values_REML <- data.frame(m_p_values_REML)
  
  m_p_values_HESD <- matrix(nrow=reps, ncol=4)
  colnames(m_p_values_HESD) <- c("p_value_10", "p_value_100", "p_value_200", "p_value_400")
  m_p_values_HESD <- data.frame(m_p_values_HESD)
  
  m_p_values_HECP <- matrix(nrow=reps, ncol=4)
  colnames(m_p_values_HECP) <- c("p_value_10", "p_value_100", "p_value_200", "p_value_400")
  m_p_values_HECP <- data.frame(m_p_values_HECP)

  for(rep in 1:reps){
    print(rep)
    # read REML results
    df_REML <- read.csv(paste(result_file_prefix, "_", rep, "_", covariance_types[i], "_trees_REML_results.csv", sep=''))
    # df_REML <- df_REML[-1,]
    df_REML <- df_REML[which(df_REML$start > position_interval[1] & df_REML$start < position_interval[2]),]
    # df_REML <- df_REML[-nrow(df_REML),]

    # read HE results
    df_HE <- read.csv(paste(result_file_prefix, "_", rep, "_", covariance_types[i], "_trees_HE_results.csv", sep=''))
    df_HE <- df_HE[-1,]
    df_HE <- df_HE[which(df_HE$start > position_interval[1] & df_HE$start < position_interval[2]),]
    # df_HE_eGRM <- df_HE_eGRM[-nrow(df_HE_eGRM),]
    df_HE$p_values_HESD_Jackknife[df_HE$p_values_HESD_Jackknife == 0] <- .Machine$double.xmin
    df_HE$p_values_HECP_Jackknife[df_HE$p_values_HECP_Jackknife == 0] <- .Machine$double.xmin

    m_results_VC[[i]]$REML[rep] <- -log10(min(df_REML$p_values, na.rm=TRUE))
    m_results_VC[[i]]$index_min_REML[rep] <- which(df_REML$p_values == min(df_REML$p_values, na.rm=TRUE))[1]
    m_results_VC[[i]]$HE_SD[rep] <- -log10(min(df_HE$p_values_HESD_Jackknife, na.rm=TRUE))
    m_results_VC[[i]]$index_min_HESD[rep] <- which(df_HE$p_values_HESD_Jackknife == min(df_HE$p_values_HESD_Jackknife, na.rm=TRUE))[1]
    m_results_VC[[i]]$HE_CP[rep] <- -log10(min(df_HE$p_values_HECP_Jackknife, na.rm=TRUE))
    m_results_VC[[i]]$index_min_HECP[rep] <- which(df_HE$p_values_HECP_Jackknife == min(df_HE$p_values_HECP_Jackknife, na.rm=TRUE))[1]

    m_results_VC[[i]]$REML_equal0.5[rep] <- sum(df_REML$p_values == 0.5)/nrow(df_REML)
    
    m_p_values_REML$p_value_10[rep] <- df_REML$p_values[10]
    m_p_values_REML$p_value_100[rep] <- df_REML$p_values[100]
    m_p_values_REML$p_value_200[rep] <- df_REML$p_values[200]
    m_p_values_REML$p_value_400[rep] <- df_REML$p_values[400]
    
    m_p_values_HESD$p_value_10[rep] <- df_HE$p_values_HESD_Jackknife[10]
    m_p_values_HESD$p_value_100[rep] <- df_HE$p_values_HESD_Jackknife[100]
    m_p_values_HESD$p_value_200[rep] <- df_HE$p_values_HESD_Jackknife[200]
    m_p_values_HESD$p_value_400[rep] <- df_HE$p_values_HESD_Jackknife[400]
    
    # m_p_values_HECP$p_value_10[rep] <- df_HE$p_values_HECP_Jackknife[10]
    # m_p_values_HECP$p_value_100[rep] <- df_HE$p_values_HECP_Jackknife[100]
    # m_p_values_HECP$p_value_200[rep] <- df_HE$p_values_HECP_Jackknife[200]
    # m_p_values_HECP$p_value_400[rep] <- df_HE$p_values_HECP_Jackknife[400]
  }
  
  cutoff_p_REML <- sort(m_results_VC[[i]][,"REML"], decreasing=T)[cutoff_rep]
  cutoff_p_HE_SD <- sort(m_results_VC[[i]][,"HE_SD"], decreasing=T)[cutoff_rep]
  cutoff_p_HE_CP <- sort(m_results_VC[[i]][,"HE_CP"], decreasing=T)[cutoff_rep]
  
  
  write.csv(cbind(cutoff_p_REML, cutoff_p_HE_SD, cutoff_p_HE_CP), file=paste("p_value_cutoffs_", covariance_types[i], sep='', ".csv"), row.names=FALSE, quote=FALSE)
  
  #--------------------------
  # hist of indeces of minimums
  #--------------------------
  pdf(paste("hist_lowest_pvalue_index_", covariance_types[i], ".pdf", sep=""), width=8, height=8)
  par(mfrow=c(2,2))
  hist(m_results_VC[[i]]$index_min_REML, xlab ="index of min p-value in ARG REML eGRM", breaks=20, main="")
  hist(m_results_VC[[i]]$index_min_HESD, xlab ="index of min p-value in ARG HE SD eGRM", breaks=20, main="")
  hist(m_results_VC[[i]]$index_min_HECP, xlab ="index of min p-value in ARG HE CP GRM", breaks=20, main="")
  dev.off()

  #--------------------------
  # qqplots
  #--------------------------
  
  pdf(paste("p_values_REML_", covariance_types[i], "_qqplot.pdf", sep=""), width=4, height=4)
  par(mfrow=c(1,1))
  
  uniform <- (runif(1000))
  qqplot(uniform, (m_p_values_REML$p_value_10), ylab="window index 10", las=2)
  abline(0,1)

  dev.off()
  
  pdf(paste("p_values_HESD_", covariance_types[i], "_qqplot.pdf", sep=""), width=15, height=5)
  par(mfrow=c(1,1))
  
  qqplot(uniform, (m_p_values_HESD$p_value_10), ylab="window index 10", las=2)
  abline(0,1)
  
  # qqplot(uniform, (m_p_values_HESD$p_value_100), ylab="tree index 100")
  # abline(0,1)
  
  dev.off()
  
  # pdf(paste("p_values_HECP_", covariance_types[i], "_qqplot.pdf", sep=""), width=15, height=5)
  # par(mfrow=c(1,1))
  # 
  # qqplot(uniform, (m_p_values_HECP$p_value_10), ylab="tree index 10")
  # abline(0,1)
  # 
  # # qqplot(uniform, (m_p_values_HECP$p_value_100), ylab="tree index 100")
  # # abline(0,1)
  # 
  # dev.off()
  
  #--------------------------
  # REML p-values
  #--------------------------
  # pdf(paste("REML_pvalue0.5_", covariance_types[i], ".pdf", sep=''), width=5, height=5)
  # hist(m_results_VC[[i]]$REML_equal0.5, xlab="Fraction of REML p-value that equal 0.5", main="")
  # dev.off()
  
}

for(rep in 1:reps){
  #read GWAS results
  df_GWAS <- read.csv(paste(result_file_prefix, "_", rep, "_GWAS_variants_results.csv", sep=''))
  df_GWAS <- df_GWAS[which(df_GWAS$start > position_interval[1] & df_GWAS$start < position_interval[2]),]
  m_results_GWAS$GWAS[rep] <- -log10(min(df_GWAS$p_value))
  m_results_GWAS$index_min_GWAS[rep] <- which(df_GWAS$p_value == min(df_GWAS$p_value))[1]
  m_results_GWAS$p_value_10[rep] <- df_GWAS$p_value[10]
  
  #ACAT
  if(run_acat){
    df_REML <- read.csv(paste(result_file_prefix, "_", rep, "_", covariance_types[i], "_trees_REML_results.csv", sep=''))
    df_REML <- df_REML[which(df_REML$start > position_interval[1] & df_REML$start < position_interval[2]),]
    
    ps_acat <- numeric(length=nrow(df_REML))
    for(r in 1:nrow(df_REML)){
      start_w <- df_REML$start[r]
      end_w <- df_REML$end[r]
      ps_acat[r] <- ACAT(df_GWAS$p_value[df_GWAS$start >= start_w & df_GWAS$start < end_w])
    }
    m_results_acat$acat[rep] <- -log10(min(ps_acat, na.rm=TRUE))
    m_results_acat$index_min_acat[rep] <- which(ps_acat == min(ps_acat, na.rm=TRUE))[1]
    m_results_acat$p_value_10[rep] <- ps_acat[10]
  }
}

pdf(paste("hist_lowest_pvalue_index_", "GWAS", ".pdf", sep=""), width=8, height=8)
hist(m_results_GWAS$index_min_GWAS, xlab ="index of min p-value of all variants GWAS", breaks=20, main="")
dev.off()

cutoff_p_GWAS <- sort(m_results_GWAS[,"GWAS"], decreasing=T)[cutoff_rep]
write.csv(cutoff_p_GWAS, file=paste("p_value_cutoffs_", "GWAS", sep='', ".csv"), row.names=FALSE, quote=FALSE)

pdf(paste("p_values_GWAS_qqplot.pdf", sep=""), width=4, height=4)
par(mfrow=c(1,1))

uniform <- (runif(1000))
qqplot(uniform, m_results_GWAS$p_value_10, ylab="window index 10", las=2)
abline(0,1)

dev.off()

#--------------------------
# ACAT
#--------------------------

if(run_acat){
  pdf(paste("hist_lowest_pvalue_index_", "acat", ".pdf", sep=""), width=8, height=8)
  hist(m_results_GWAS$index_min_GWAS, xlab ="index of min p-value of all variants acat", breaks=20, main="")
  dev.off()
  
  cutoff_p_acat <- sort(m_results_acat[,"acat"], decreasing=T)[cutoff_rep]
  write.csv(cutoff_p_acat, file=paste("p_value_cutoffs_", "acat", sep='', ".csv"), row.names=FALSE, quote=FALSE)
  
  pdf(paste("p_values_acat_qqplot.pdf", sep=""), width=4, height=4)
  par(mfrow=c(1,1))
  
  uniform <- (runif(1000))
  qqplot(uniform, m_results_GWAS$p_value_10, ylab="window index 10", las=2)
  abline(0,1)
  
  dev.off()
}


#--------------------------
# cutoff plot all
#--------------------------

# pdf("p_values_cutoff_separate.pdf", width=10, height=5)
# 
# par(mfrow=c(1,3))
# par(pty="s")
# plot(sort(m_results[,"REML"], decreasing=T), col="dodgerblue", yaxt='n', ylim=c(1,6.5), xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
# lines(sort(m_results[,"GWAS"], decreasing=T), col="black")
# axis(side=2, las=2)
# legend(legend=c("GWAS", "REML"), x="topright", ncol=1, col=c("black", "dodgerblue"), lty=1, bty='n')
# abline(v=0.05*reps, lty=2)
# 
# plot(sort(m_results[,"HE_SD"], decreasing=T), col="red", yaxt='n', xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
# lines(sort(m_results[,"HE_CP"], decreasing=T), col="orange2")
# axis(side=2, las=2)
# legend(legend=c("HE_SD","HE_CP"), x="topright", ncol=1, col=c("red", "orange2"), lty=1, bty='n')
# abline(v=0.05*reps, lty=2)
# dev.off()


pdf("p_values_cutoff_all.pdf", width=5, height=5)

par(pty="s")
plot(sort(m_results_VC[[1]][,"REML"], decreasing=T), col="dodgerblue", yaxt='n', ylim=c(0.5,6), xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l')
lines(sort(m_results_VC[[2]][,"REML"], decreasing=T), col="dodgerblue", lty=2)
lines(sort(m_results_VC[[1]][,"HE_SD"], decreasing=T), col="red")
lines(sort(m_results_VC[[2]][,"HE_SD"], decreasing=T), col="red", lty=2)
lines(sort(m_results_VC[[1]][,"HE_CP"], decreasing=T), col="orange2")
lines(sort(m_results_VC[[2]][,"HE_CP"], decreasing=T), col="orange2", lty=2)
lines(sort(m_results_GWAS[,"GWAS"], decreasing=T), col="black")

axis(side=2, las=2)
legend(legend=c(paste("GWAS:", cutoff_p_GWAS), 
                paste("REML eGRM:", cutoff_p_REML), 
                paste("REML GRM:", cutoff_p_REML), 
                paste("HE_SD eGRM:", cutoff_p_HE_SD), 
                paste("HE_SD GRM:", cutoff_p_HE_SD), 
                paste("HE_CP eGRM:", cutoff_p_HE_CP), 
                paste("HE_CP GRM:", cutoff_p_HE_CP)), 
       x="topright", ncol=1, col=c("black", "dodgerblue","dodgerblue", "red", "red","orange2", "orange2"), lty=c(1,rep(c(1,2), 3)), bty='n')
abline(v=cutoff_rep, lty=2)

dev.off()


#--------------------------
# cutoff plot paper
#--------------------------


pdf("p_values_cutoff_paper.pdf", width=5, height=5)

par(pty="s")
plot(sort(m_results_VC[[1]][,"REML"], decreasing=T), col="dodgerblue", yaxt='n', ylim=c(0.5,6), xlim=c(1,reps), ylab=expression("-log"[10]*"(p)"), xlab="ordered simulation number", type='l', bty='n')
lines(sort(m_results_VC[[2]][,"REML"], decreasing=T), col="maroon2", lty=1)
lines(sort(m_results_GWAS[,"GWAS"], decreasing=T), col="orange2")
lines(sort(m_results_acat[,"acat"], decreasing=T), col="black")

axis(side=2, las=2)
legend(legend=c("GWAS", 
                "local eGRM",  
                "local GRM",
                "ACAT-V"), 
                x="topright", ncol=1, col=c("orange2", "dodgerblue","maroon2", "black"), lty=c(1), bty='n')
abline(v=cutoff_rep, lty=2)

dev.off()
