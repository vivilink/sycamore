pdf("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/null_p_value_distributions.pdf", width=7, height=7)
par(mfrow=c(2,2))
par(pty="s")


setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/relate_trees/window_based/5k")

n_reps <- 300
p_values <- c()

for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}


qqplot(runif(5000), p_values, xlim=c(0,1), ylim=c(0,1), main="Relate trees, local eGRM", xlab="U(0,1)", ylab="p-values", bty='n')
abline(a=0, b=1)



setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k")

n_reps <- 300
p_values <- c()

for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}


qqplot(runif(5000), p_values, xlim=c(0,1), ylim=c(0,1), main="True trees, local eGRM", xlab="U(0,1)", ylab="p-values", bty='n')
abline(a=0, b=1)


setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/relate_trees/window_based/5k")

n_reps <- 300
p_values <- c()

for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}


qqplot(runif(5000), p_values, xlim=c(0,1), ylim=c(0,1), main="Relate trees, local GRM", xlab="U(0,1)", ylab="p-values", bty='n')
abline(a=0, b=1)



setwd("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k")

n_reps <- 300
p_values <- c()

for(r in 1:n_reps){
  df <- read.csv(paste("cutoff_sims_",r,"_GRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}


qqplot(runif(5000), p_values, xlim=c(0,1), ylim=c(0,1), main="True trees, local GRM", xlab="U(0,1)", ylab="p-values", bty='n')
abline(a=0, b=1)


dev.off()
