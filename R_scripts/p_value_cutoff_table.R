
base_dir <- "/data/ARGWAS/experiments_cutoff_N2K/"
setwd("/data/ARGWAS/power_sims/stdpopsim")

options(scipen = 100, digits = 4)
run_acat <- TRUE

t <- matrix(ncol=4)
colnames(t) <- c("tree type", "testing window size","method", "p-value cutoff")

for(tree_type in c("true_trees", "relate_trees")){ #, , true_trees
  for(region_type in c("window_based")){
    for(ws_testing in c("5k","10k")){ #, "20k", "50k" , "10k", ,"10k"
      power_results_aH <- data.frame()
      folder=paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneVariant/",variant_freq, "/eGRM_GRM/", region_type, "/", ws_causal, "/" ,sep="")
      print(paste("analyzing folder", folder))
      cutoff_GWAS <- read.csv(paste(base_dir, "diploid/GRM_eGRM/", tree_type, "/", region_type, "/", ws_testing, "/p_value_cutoffs_GWAS.csv", sep=''))$x
      cutoff_acat <- read.csv(paste(base_dir, "diploid/GRM_eGRM/", tree_type, "/", region_type, "/", ws_testing, "/p_value_cutoffs_acat.csv", sep=''))$x
      cutoff_eGRM <- read.csv(paste(base_dir, "diploid/GRM_eGRM/", tree_type,  "/", region_type, "/", ws_testing, "/p_value_cutoffs_eGRM.csv", sep=''))$cutoff_p_REML
      cutoff_GRM <- read.csv(paste(base_dir, "diploid/GRM_eGRM/", tree_type,  "/", region_type, "/", ws_testing, "/p_value_cutoffs_GRM.csv", sep=''))$cutoff_p_REML
      
      t <- rbind(t, c(tree_type, ws_testing, "GWAS", round(cutoff_GWAS,4)))
      t <- rbind(t, c(tree_type, ws_testing, "ACAT-V", round(cutoff_acat,4)))
      t <- rbind(t, c(tree_type, ws_testing, "local eGRM", round(cutoff_eGRM,4)))
      t <- rbind(t, c(tree_type, ws_testing, "local GRM", round(cutoff_GRM,4)))
      
    }
  }
}

t <- t[-1,]
t[which(t[,1] == "true_trees" & t[,3] != "local eGRM"),1] <- "all variants"
t[which(t[,1] == "true_trees" & t[,3] == "local eGRM"),1] <- "true trees"

t[which(t[,1] == "relate_trees" & t[,3] != "local eGRM"),1] <- "typed variants"
t[which(t[,1] == "relate_trees" & t[,3] == "local eGRM"),1] <- "Relate trees"

t[which(t == "5k")] <- "5kb"
t[which(t == "10k")] <- "10kb"

print(xtable(t, type = "latex"), file = paste(base_dir,"p_value_cutoffs.tex", sep=''))



