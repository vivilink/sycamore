pdf("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/null_p_value_distributions.pdf", width=7, height=7)
par(mfrow=c(2,2))
par(pty="s")
par(mar=c(3.5,3.5,3.5,3.5))


plot_qq <- function(p_values, MAIN){
  unif <- runif(5000)
  qqplot(unif, p_values, xlim=c(0,1), ylim=c(0,1), main=MAIN, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', pch=20)
  axis(side=1, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1)) 
  axis(side=2, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1), las=2)
  title(ylab="p", line=2.8)
  title(xlab="Uniform(0,1)", line=2.2)
  abline(a=0, b=1)
}

setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/relate_trees/window_based/5k")
n_reps <- 300
p_values <- c()
for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}
plot_qq(p_values=p_values, MAIN="relate trees local eGRM")


setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/relate_trees/window_based/5k")
n_reps <- 300
p_values <- c()
for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_GRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}
plot_qq(p_values=p_values, MAIN="typed variants local GRM")


setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k")
n_reps <- 300
p_values <- c()
for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}
plot_qq(p_values=p_values, MAIN="true trees local eGRM")


setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k")
n_reps <- 300
p_values <- c()
for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_GRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}
plot_qq(p_values=p_values, MAIN="all variants local GRM")
dev.off()
