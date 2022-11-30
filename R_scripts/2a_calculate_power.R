setwd("/data/ARGWAS/power_sims/stdpopsim")
source("power_one_experiment_oneVariant.R")
library("pwr")
options(scipen = 100, digits = 4)
hs_all <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
region_type="window_based"

#--------------------------
#read low freq true trees
#--------------------------
power_results_trueTrees_lowFreq <- data.frame()
for(hs in hs_all){
  folder=paste("true_trees/oneVariant/rareVariant/eGRM_GRM/",region_type,"/", sep='')
  results_file <- paste(folder, "h", hs, "/power_results.txt", sep='')
  print(paste(hs, file.exists(results_file)))
  # if(file.exists(results_file) == FALSE){
    power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type="relate_trees", region_type=region_type)
    print("finished creating results file")
  # }
  t <- read.table(results_file, header=TRUE)
  power_results_trueTrees_lowFreq <- rbind(power_results_trueTrees_lowFreq, t)
}
rownames(power_results_trueTrees_lowFreq) = hs_all

#--------------------------
#read high freq true trees
#--------------------------
power_results_trueTrees_highFreq <- data.frame()
for(hs in hs_all){
  folder=paste("true_trees/oneVariant/commonVariant/eGRM_GRM/",region_type,"/", sep='')
  results_file <- paste(folder, "h", hs, "/power_results.txt", sep='')
  print(paste(hs, file.exists(results_file)))
  if(file.exists(results_file) == FALSE){
    power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type="true_trees", region_type=region_type)
    print("finished creating results file")
  }
  t <- read.table(results_file, header=TRUE)
  power_results_trueTrees_highFreq <- rbind(power_results_trueTrees_highFreq, t)
}
rownames(power_results_trueTrees_highFreq) = hs_all

#--------------------------
#read low freq relate trees
#--------------------------
power_results_relateTrees_lowFreq <- data.frame()
for(hs in hs_all){
  folder=paste("relate_trees/oneVariant/rareVariant/eGRM_GRM/",region_type,"/", sep='')
  results_file <- paste(folder, "h", hs, "/power_results.txt", sep='')
  print(paste(hs, file.exists(results_file)))
  # if(file.exists(results_file) == FALSE){
    power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type="relate_trees", region_type=region_type)
    print("finished creating results file")
  # }
  t <- read.table(results_file, header=TRUE)
  power_results_relateTrees_lowFreq <- rbind(power_results_relateTrees_lowFreq, t)
}
rownames(power_results_relateTrees_lowFreq) = hs_all

#--------------------------
#read high freq relate trees
#--------------------------
power_results_relateTrees_highFreq <- data.frame()
for(hs in hs_all){
  folder=paste("relate_trees/oneVariant/commonVariant/eGRM_GRM/",region_type,"/", sep='')
  results_file <- paste(folder, "h", hs, "/power_results.txt", sep='')
  print(paste(hs, file.exists(results_file)))
  # if(file.exists(results_file) == FALSE){
    power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type="relate_trees", region_type=region_type)
    print("finished creating results file")
  # }
  t <- read.table(results_file, header=TRUE)
  power_results_relateTrees_highFreq <- rbind(power_results_relateTrees_highFreq, t)
}
rownames(power_results_relateTrees_highFreq) = hs_all

#-----------------
#plot according to freq
#-----------------

plot_freq_error_bars <- function(power, color){
  points(y=power, x=x_pos-offset, pch=1, col=color)
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos-offset, y0=power - std, y1=power + std, col="orange2")
}

plot_freq <- function(freq, trueTreeResult, relateTreeResults){
  
  trueTreeResult <- trueTreeResult[rownames(trueTreeResult) %in% hs_all,]
  relateTreeResults <- relateTreeResults[rownames(relateTreeResults) %in% hs_all,]
  
  pdf(paste("power_results_all_trees_af", freq,"_", region_type, ".pdf", sep=''), height=10, width=10)
  par(mfrow=c(2,2))
  
  #all tree types
  plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
  axis(side=1, at=1:length(hs_all), labels=hs_all)
  
  lines(trueTreeResult$power_GWAS, col="orange2", lty=1, type = "b", pch=20)
  lines(trueTreeResult$power_REML_eGRM, col="dodgerblue", lty=1, type = "b", pch=20)
  lines(trueTreeResult$power_REML_GRM, col="maroon2", lty=1, type = "b", pch=20)
  # lines(power_results_lowFreq$power_HE_SD, col="orange2", lty=2)
  # lines(power_results_lowFreq$power_HE_CP, col="red", lty=2)
  
  lines(relateTreeResults$power_GWAS, col="orange2", lty=2, type = "b", pch=20)
  lines(relateTreeResults$power_REML_eGRM, col="dodgerblue", lty=2, type = "b", pch=20)
  lines(relateTreeResults$power_REML_GRM, col="maroon2", lty=2, type = "b", pch=20)
  # lines(power_results_highFreq$power_HE_SD, col="orange2")
  # lines(power_results_highFreq$power_HE_CP, col="red")
  
  legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))
  
  # #all tree types within region
  # plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="within 1000bp from causal pos significant")
  # axis(side=1, at=1:length(hs_all), labels=hs_all)
  # 
  # lines(trueTreeResult$power_GWAS_region, col="orange2", lty=1, type = "b", pch=20)
  # lines(trueTreeResult$power_REML_region, col="dodgerblue", lty=1, type = "b", pch=20)
  # 
  # lines(relateTreeResults$power_GWAS_region, col="orange2", lty=2, type = "b", pch=20)
  # lines(relateTreeResults$power_REML_region, col="dodgerblue", lty=2, type = "b", pch=20)
  
  # legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))

  #only relate trees
  plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
  axis(side=1, at=1:length(hs_all), labels=hs_all)
  
  lines(relateTreeResults$power_GWAS, col="orange2", lty=2, type = "b", pch=20)
  lines(relateTreeResults$power_REML_eGRM, col="dodgerblue", lty=2, type = "b", pch=20)
  lines(relateTreeResults$power_REML_GRM, col="maroon2", lty=2, type = "b", pch=20)
  # lines(power_results_highFreq$power_HE_SD, col="orange2")
  # lines(power_results_highFreq$power_HE_CP, col="red")
  
  legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))
  
  #only true trees
  plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
  axis(side=1, at=1:length(hs_all), labels=hs_all)
  
  lines(trueTreeResult$power_GWAS, col="orange2", lty=1, type = "b", pch=20)
  lines(trueTreeResult$power_REML_eGRM, col="dodgerblue", lty=1, type = "b", pch=20)
  lines(trueTreeResult$power_REML_GRM, col="maroon2", lty=1, type = "b", pch=20)
  
  dev.off()
}

plot_freq(freq=0.02, trueTreeResult=power_results_trueTrees_lowFreq, relateTreeResult=power_results_relateTrees_lowFreq)
plot_freq(freq=0.2, trueTreeResult=power_results_trueTrees_highFreq, relateTreeResult=power_results_relateTrees_highFreq)

#--------------------------
# plot power error bars
#--------------------------
setwd("/data/ARGWAS/power_sims/stdpopsim/")
nreps=200

hs <- c(0.02)
# hs <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
offset <- 0

pdf("power_oneVariant_erorBars.pdf", height=5, width=5)
plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, 3+offset), ylab="power", xlab='', xaxt='n')
for(h in hs){

  #GWAS relate
  power_results <- read.table(paste("/data/ARGWAS/power_sims/stdpopsim/relate_trees/oneVariant/commonVariant/eGRM_GRM/h", h, "/power_results.txt", sep=''), header=TRUE)
  points(y=power_results$power_GWAS, x=1-offset, pch=1, col="orange2")
  power <- power_results$power_GWAS
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=1-offset, y0=power - std, y1=power + std, col="orange2")	
  
  #REML eGRM relate
  power_results <- read.table(paste("/data/ARGWAS/power_sims/stdpopsim/relate_trees/oneVariant/commonVariant/eGRM_GRM/h", h, "/power_results.txt", sep=''), header=TRUE)
  points(y=power_results$power_REML_eGRM, x=2+offset, pch=1, col="dodgerblue")
  power <- power_results$power_REML_eGRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=2+offset, y0=power - std, y1=power + std, col="dodgerblue")
  
  #REML true tree
  power_results <- read.table(paste("/data/ARGWAS/power_sims/stdpopsim/true_trees/oneVariant/commonVariant/eGRM_GRM/h", h, "/power_results.txt", sep=''), header=TRUE)
  points(y=power_results$power_REML_GRM, x=3+offset, pch=19, col="dodgerblue")
  power <- power_results$power_REML_GRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=3+offset, y0=power - std, y1=power + std, col="dodgerblue")
  
}

# legend(x="bottomright", legend=c("local REML eGRM", "local REML GRM", "GWAS"), pch=c(1, 1), col=c("dodgerblue","black", "orange2"), bty='n')
dev.off()

#--------------------------
# power files allelic heterogeneity
#--------------------------

setwd("/data/ARGWAS/power_sims/stdpopsim")

source("/home/linkv/git/argwas/R_scripts/2b_calculate_power_one_experiment.R")
library("pwr")
options(scipen = 100, digits = 4)
hs_all <- c(0.02, 0.04,  0.06,  0.08, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
# hs_all <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 


#read low freq true trees
for(tree_type in c("relate_trees", "true_trees")){ #, 
  for(region_type in c("window_based")){
    for(ws in c("10k", "5k")){ #, "20k", "50k"
      power_results_aH <- data.frame()
      for(hs in hs_all){
        # folder=paste("/data/ARGWAS/power_sims/stdpopsim/high_mut_trees/oneTree/eGRM_and_GRM/", tree_type, "/", region_type, "/", ws, "/",sep="")
        folder=paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws, "/",sep="")
    
        results_file <- paste(folder, "h", hs, "/power_results.txt", sep='')
        print(paste(hs, file.exists(results_file)))
        if(file.exists(results_file) == FALSE){
          power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type=tree_type, region_type=region_type, window_size=ws)
          # power_one_experiment(hsquared = hs, REPS = 200, folder=folder, tree_type="high_mut_trees")
          # print("finished creating results file")
        }
        t <- read.table(results_file, header=TRUE)
        power_results_aH <- rbind(power_results_aH, t)
      }
    }
  }
}

