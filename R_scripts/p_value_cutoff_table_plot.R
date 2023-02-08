#---------------------------
# latex table
#---------------------------

base_dir <- "/data/ARGWAS/experiments_cutoff_N2K/"
setwd(base_dir)

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


#---------------------------
# p-value cutoff plot
#---------------------------

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

plot_lines <- function(folder, LTY){
  p_values <- read.csv(file=paste(folder, "p_values_replicates_eGRM.csv",sep=''))
  lines(sort(p_values$REML_min_p_value, decreasing=TRUE), col=blu, lty=LTY)
  
  p_values <- read.csv(file=paste(folder, "p_values_replicates_GRM.csv",sep=''))
  lines(sort(p_values$REML_min_p_value, decreasing=TRUE), col=pin, lty=LTY)
  
  p_values <- read.csv(file=paste(folder, "p_values_replicates_GWAS.csv",sep=''))
  lines(sort(p_values$GWAS_min_p_value, decreasing=TRUE), col=org, lty=LTY)
  
  p_values <- read.csv(file=paste(folder, "p_values_replicates_acat.csv",sep=''))
  lines(sort(p_values$ACAT_min_p_value, decreasing=TRUE), col="black", lty=LTY)
}



pdf("p_value_cutoffs_true_relate.pdf", width=5, height = 5)
reps <- 300
cutoff_rep <- 0.05 * reps

plot(0, type='n', yaxt='n', ylim=c(0.5,6), xlim=c(1,reps), ylab='', xlab="",  bty='n')
title(ylab=expression("-log"[10]*"(p)"), line=2)
title(xlab="ordered simulation number", line=2.2)
axis(side=2, las=2)



#true trees
folder <- "/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/true_trees/window_based/5k/"
plot_lines(folder=folder, LTY=2)

#relate trees
folder <- "/data/ARGWAS/experiments_cutoff_N2K/diploid/GRM_eGRM/relate_trees/window_based/5k/"
plot_lines(folder=folder, LTY=1)

abline(v=cutoff_rep, lty=1)

legend(x="topright", legend=c("true trees / all variants","Relate / typed variants", "GWAS", "ACAT-V", "local GRM", "local eGRM"), col=c("gray", "gray", org, "black", pin , blu), bty='n', pch=c(NA, NA, 15, 15, 15, 15), lty=c(2,1, NA, NA, NA, NA))

dev.off()

