library("xtable")

base_dir <- "/data/ARGWAS/power_sims/stdpopsim/"
setwd(base_dir)
run_acat <- TRUE
# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

plot_one <- function(dir, max_y, max_x, LTY, MAIN, tree_type){
  
  t <- matrix(nrow=4, ncol=7)
  colnames(t) <- c("tree_type", "causal_variants","method", "mean", "median", "q25", "q75")

  d_GWAS <- read.table(paste(dir,"/association_results_GWAS.txt", sep=''), header=TRUE)$distance_min_p_to_causal
  t[1,] <- c(tree_type, MAIN, "GWAS", round(mean(d_GWAS)), round(quantile(d_GWAS, probs=c(0.5, 0.25, 0.75))))
  
  d_acat <- read.table(paste(dir,"/association_results_acat.txt", sep=''), header=TRUE)$distance_min_p_to_causal
  t[2,] <- c(tree_type, MAIN, "ACAT-V", round(mean(d_acat)), round(quantile(d_acat, probs=c(0.5, 0.25, 0.75))))
  
  d_eGRM <- read.table(paste(dir,"/association_results_eGRM.txt", sep=''), header=TRUE)$distance_min_p_to_causal_REML
  t[3,] <- c(tree_type, MAIN, "local eGRM", round(mean(d_eGRM)), round(quantile(d_eGRM, probs=c(0.5, 0.25, 0.75))))
  
  d_GRM <- read.table(paste(dir,"/association_results_GRM.txt", sep=''), header=TRUE)$distance_min_p_to_causal_REML
  t[4,] <- c(tree_type, MAIN, "local GRM", round(mean(d_GRM)), round(quantile(d_GRM, probs=c(0.5, 0.25, 0.75))))
  
  
  #plot GWAS
  par(mar=c(5,7,3,3))

  plot(density(d_GWAS, na.rm = TRUE, bw="SJ"), col=org, 
       ylim=c(0,max_y), bty='n', las=2, xaxt='n', ylab='', xlim=c(0,max_x),
       xlab="Distance [kb]", main="", lty=LTY)
  axis(side=1, at=seq(0,max(max_x),by=100000), labels=seq(0,max(max_x),by=100000) / 1000)
  title(ylab = "Density", line = 5) 
  
  #add the rest relate
  lines(density(d_eGRM, na.rm = TRUE, bw="SJ"), col=blu, lty=LTY)
  lines(density(d_GRM, na.rm = TRUE, bw="SJ"), col=pin, lty=LTY)
  if(run_acat){
    lines(density(d_acat, na.rm = TRUE), col="black", lty=LTY)
  }
  return(t)
}

pdf("distance_causal_window.pdf", width=10, height=5)
layout(mat=matrix(c(1,2,3,4,5,6), ncol=3))


#AH

propCausal <- 0.1
relate_dir <- paste(base_dir, "relate_trees/oneRegion/eGRM_GRM/window_based/5k/tested5k/propCausal", propCausal, "/h0.02/", sep='')
true_dir <- paste(base_dir, "true_trees/oneRegion/eGRM_GRM/window_based/5k/tested5k/propCausal", propCausal, "/h0.02/", sep='')
relate_stats <- plot_one(dir=relate_dir, max_y=0.000005, max_x=600000, LTY=1, MAIN=paste("prop causal", propCausal), tree_type="relate_trees")
true_stats <- plot_one(dir=true_dir, max_y=0.00003, max_x=600000, LTY=2, MAIN=paste("prop causal", propCausal), tree_type="true_trees")
t <- rbind(relate_stats, true_stats)

propCausal <- 0.2
relate_dir <- paste(base_dir, "relate_trees/oneRegion/eGRM_GRM/window_based/5k/tested5k/propCausal", propCausal, "/h0.02/", sep='')
true_dir <- paste(base_dir, "true_trees/oneRegion/eGRM_GRM/window_based/5k/tested5k/propCausal", propCausal, "/h0.02/", sep='')
relate_stats <- plot_one(dir=relate_dir, max_y=0.000005, max_x=600000, LTY=1, MAIN=paste("prop causal", propCausal), tree_type="relate_trees")
true_stats <- plot_one(dir=true_dir, max_y=0.00003, max_x=600000, LTY=2, MAIN=paste("prop causal", propCausal), tree_type="true_trees")
t <- rbind(t, relate_stats, true_stats)

#rare Variant
relate_dir <- paste(base_dir, "/relate_trees/oneVariant/rareVariant/eGRM_GRM/window_based/10k/h0.02/", sep='')
true_dir <- paste(base_dir, "/true_trees/oneVariant/rareVariant/eGRM_GRM/window_based/10k/h0.02/", sep='')
relate_stats <- plot_one(dir=relate_dir, max_y=0.000006, max_x=600000, LTY=1, MAIN="rare variant", tree_type="relate_trees")
true_stats <- plot_one(dir=true_dir, max_y=0.00002, max_x=600000, LTY=2, MAIN="rare variant", tree_type="true_trees")
t <- rbind(t, relate_stats, true_stats)

#rename stuff  
t[which(t[,1] == "true_trees" & t[,3] != "local eGRM"),1] <- "all variants"
t[which(t[,1] == "true_trees" & t[,3] == "local eGRM"),1] <- "true trees"

t[which(t[,1] == "relate_trees" & t[,3] != "local eGRM"),1] <- "typed variants"
t[which(t[,1] == "relate_trees" & t[,3] == "local eGRM"),1] <- "Relate trees"

t[which(t[,2] == "rare variant"),2] <- "causal variant frequency 0.02"

t[which(t == "5k")] <- "5kb"
t[which(t == "10k")] <- "10kb"

legend(x="topright", legend=c("Relate / typed variants", "true trees / all variants", "local eGRM", "local GRM", "GWAS", "ACAT-V"), lty=c(1,2,NA, NA, NA, NA), pch=c(NA, NA, 15, 15, 15, 15), col=c("gray","gray",blu, pin, org, "black"), bty='n')
dev.off()

#write latex table
print(xtable(t, type = "latex"), file = paste(base_dir,"distance_causal_window.tex", sep=''))
