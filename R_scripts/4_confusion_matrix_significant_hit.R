library(xtable)
base_dir <- "/data/ARGWAS/power_sims/stdpopsim/"

make_confusion_m <- function(dir, tree_type, window_size_causal, window_size_testing, propCausal, null_sims = FALSE){
  if(null_sims){
    reps <- 300
  } else {
    reps <- 200
  }
  
  if(null_sims == TRUE){
    p_GWAS <- read.table(paste(dir,"/p_values_replicates_GWAS.csv", sep=''), sep=',', header=TRUE)$GWAS_min_p_value
    p_acat <- read.table(paste(dir,"/p_values_replicates_acat.csv", sep=''), sep=',', header=TRUE)$ACAT_min_p_value
    p_eGRM <- read.table(paste(dir,"/p_values_replicates_eGRM.csv", sep=''), sep=',', header=TRUE)$REML_min_p_value
    p_GRM <- read.table(paste(dir,"/p_values_replicates_GRM.csv", sep=''), sep=',', header=TRUE)$REML_min_p_value
  } else{
    p_GWAS <- read.table(paste(dir,"/association_results_GWAS.txt", sep=''), header=TRUE)$GWAS
    p_acat <- read.table(paste(dir,"/association_results_acat.txt", sep=''), header=TRUE)$acat
    p_eGRM <- read.table(paste(dir,"/association_results_eGRM.txt", sep=''), header=TRUE)$REML
    p_GRM <- read.table(paste(dir,"/association_results_GRM.txt", sep=''), header=TRUE)$REML
  }
  
  cutoff_GWAS <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, "/p_value_cutoffs_GWAS.csv", sep=''))$x
  cutoff_acat <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, "/p_value_cutoffs_acat.csv", sep=''))$x
  cutoff_eGRM  <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, "/p_value_cutoffs_eGRM.csv", sep=''))$cutoff_p_REML
  cutoff_GRM <- read.csv(paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, "/p_value_cutoffs_GRM.csv", sep=''))$cutoff_p_REML
  
  
  sig_GWAS <- which(p_GWAS > cutoff_GWAS)
  sig_acat <- which(p_acat > cutoff_acat)
  sig_eGRM <- which(p_eGRM > cutoff_eGRM)
  sig_GRM <- which(p_GRM > cutoff_GRM)
  
  #enum
  GWAS <- 1
  acat <- 2
  eGRM <- 3
  GRM <- 4
  m <- matrix(ncol=4, nrow=4)
  colnames(m) <- c("GWAS", "ACAT", "local eGRM", "local GRM")
  rownames(m) <- c("GWAS", "ACAT", "local eGRM", "local GRM")
  
  m[GWAS,acat] <- length(intersect(sig_GWAS, sig_acat)) / reps
  m[GWAS,eGRM] <- length(intersect(sig_GWAS, sig_eGRM)) / reps
  m[GWAS,GRM] <- length(intersect(sig_GWAS, sig_GRM)) / reps
  m[acat,eGRM] <- length(intersect(sig_acat, sig_eGRM)) / reps
  m[acat,GRM] <- length(intersect(sig_acat, sig_GRM)) / reps
  m[eGRM,GRM] <- length(intersect(sig_eGRM, sig_GRM)) / reps
  
  m[GWAS,GWAS] <- length(sig_GWAS) / reps
  m[acat,acat] <- length(sig_acat) / reps
  m[eGRM,eGRM] <- length(sig_eGRM) / reps
  m[GRM,GRM] <- length(sig_GRM) / reps
  
  #write latex table
  print(xtable(m, type = "latex"), file = paste(base_dir,"confusion_matrix_significant_hits_", tree_type, "_windowTesting", window_size_testing, "_propCausal",propCausal,".tex", sep=''), floating = FALSE)
  
}

# AH

tree_type <- "relate_trees"
window_size_testing <- "5k"
window_size_causal <- "5k"
propCausal <- 0.1
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneRegion/eGRM_GRM/window_based/", window_size_causal, "/tested",window_size_testing,"/propCausal",propCausal,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)


tree_type <- "true_trees"
window_size_testing <- "5k"
window_size_causal <- "5k"
propCausal <- 0.1
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneRegion/eGRM_GRM/window_based/", window_size_causal, "/tested",window_size_testing,"/propCausal",propCausal,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)


# one causal variant

tree_type <- "relate_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "commonVariant"
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneVariant/commonVariant/eGRM_GRM/window_based/tested",window_size_testing,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)

tree_type <- "true_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "commonVariant"
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneVariant/commonVariant/eGRM_GRM/window_based/tested",window_size_testing,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)

tree_type <- "relate_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "commonVariant"
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneVariant/commonVariant/eGRM_GRM/window_based/tested",window_size_testing,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)

tree_type <- "true_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "commonVariant"
dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type,"/oneVariant/commonVariant/eGRM_GRM/window_based/tested",window_size_testing,"/h0.02", sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal)


# null simulations

tree_type <- "relate_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "NULLsims"
dir <- paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal, null_sims=TRUE)

tree_type <- "true_trees"
window_size_testing <- "5k"
window_size_causal <- NA
propCausal <- "NULLsims"
dir <- paste("/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/", tree_type, "/window_based/", window_size_testing, sep='')
make_confusion_m(dir=dir, tree_type=tree_type, window_size_causal=window_size_causal, window_size_testing=window_size_testing, propCausal=propCausal, null_sims=TRUE)