# setwd("/data/ARGWAS/power_sims/stdpopsim")
# source("~/git/argwas/R_scripts/2b_calculate_power_one_experiment.R")
# library("pwr")
# options(scipen = 100, digits = 4)
# hs_all <- c(0.02) #0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1       , 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
# region_type <- "window_based"
# window_size <- "10k"
# 
# #--------------------------
# #read low freq true trees
# #--------------------------
# power_results_trueTrees_lowFreq <- data.frame()
# for(hsquared in hs_all){
#   folder=paste("true_trees/oneVariant/rareVariant/eGRM_GRM/",region_type,"/",window_size, "/", sep='')
#   results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
#   print(paste(hsquared, file.exists(results_file)))
#   # if(file.exists(results_file) == FALSE){
#     power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type="true_trees", region_type=region_type, window_size = window_size)
#     print("finished creating results file")
#   # }
#   t <- read.table(results_file, header=TRUE)
#   power_results_trueTrees_lowFreq <- rbind(power_results_trueTrees_lowFreq, t)
# }
# rownames(power_results_trueTrees_lowFreq) = hs_all
# 
# #--------------------------
# #read high freq true trees
# #--------------------------
# power_results_trueTrees_highFreq <- data.frame()
# for(hsquared in hs_all){
#   folder=paste("true_trees/oneVariant/commonVariant/eGRM_GRM/",region_type,"/",window_size,"/", sep='')
#   results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
#   print(paste(hsquared, file.exists(results_file)))
#   # if(file.exists(results_file) == FALSE){
#     power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type="true_trees", region_type=region_type, window_size = window_size)
#     print("finished creating results file")
#   # }
#   t <- read.table(results_file, header=TRUE)
#   power_results_trueTrees_highFreq <- rbind(power_results_trueTrees_highFreq, t)
# }
# rownames(power_results_trueTrees_highFreq) = hs_all
# 
# #--------------------------
# #read low freq relate trees
# #--------------------------
# power_results_relateTrees_lowFreq <- data.frame()
# for(hsquared in hs_all){
#   folder=paste("relate_trees/oneVariant/rareVariant/eGRM_GRM/",region_type,"/", window_size, "/", sep='')
#   results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
#   print(paste(hsquared, file.exists(results_file)))
#   # if(file.exists(results_file) == FALSE){
#     power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type="relate_trees", region_type=region_type, window_size = window_size)
#     print("finished creating results file")
#   # }
#   t <- read.table(results_file, header=TRUE)
#   power_results_relateTrees_lowFreq <- rbind(power_results_relateTrees_lowFreq, t)
# }
# rownames(power_results_relateTrees_lowFreq) = hs_all
# 
# #--------------------------
# #read high freq relate trees
# #--------------------------
# power_results_relateTrees_highFreq <- data.frame()
# for(hsquared in hs_all){
#   folder=paste("relate_trees/oneVariant/commonVariant/eGRM_GRM/",region_type,"/", window_size, "/", sep='')
#   results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
#   print(paste(hsquared, file.exists(results_file)))
#   # if(file.exists(results_file) == FALSE){
#     power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type="relate_trees", region_type=region_type, window_size = window_size)
#     print("finished creating results file")
#   # }
#   t <- read.table(results_file, header=TRUE)
#   power_results_relateTrees_highFreq <- rbind(power_results_relateTrees_highFreq, t)
# }
# rownames(power_results_relateTrees_highFreq) = hs_all
# 
# #-----------------
# #plot according to freq
# #-----------------
# 
# plot_freq_error_bars <- function(power, color){
#   points(y=power, x=x_pos-offset, pch=1, col=color)
#   std <- sqrt(power*(1-power)/nreps)
#   segments(x0=x_pos-offset, y0=power - std, y1=power + std, col="orange2")
# }
# 
# plot_freq <- function(freq, trueTreeResult, relateTreeResults){
#   
#   trueTreeResult <- trueTreeResult[rownames(trueTreeResult) %in% hs_all,]
#   relateTreeResults <- relateTreeResults[rownames(relateTreeResults) %in% hs_all,]
#   
#   pdf(paste("power_results_all_trees_af", freq,"_", region_type, ".pdf", sep=''), height=10, width=10)
#   par(mfrow=c(2,2))
#   
#   #all tree types
#   plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
#   axis(side=1, at=1:length(hs_all), labels=hs_all)
#   
#   lines(trueTreeResult$power_GWAS, col="orange2", lty=1, type = "b", pch=20)
#   lines(trueTreeResult$power_REML_eGRM, col="dodgerblue", lty=1, type = "b", pch=20)
#   lines(trueTreeResult$power_REML_GRM, col="maroon2", lty=1, type = "b", pch=20)
#   # lines(power_results_lowFreq$power_HE_SD, col="orange2", lty=2)
#   # lines(power_results_lowFreq$power_HE_CP, col="red", lty=2)
#   
#   lines(relateTreeResults$power_GWAS, col="orange2", lty=2, type = "b", pch=20)
#   lines(relateTreeResults$power_REML_eGRM, col="dodgerblue", lty=2, type = "b", pch=20)
#   lines(relateTreeResults$power_REML_GRM, col="maroon2", lty=2, type = "b", pch=20)
#   # lines(power_results_highFreq$power_HE_SD, col="orange2")
#   # lines(power_results_highFreq$power_HE_CP, col="red")
#   
#   legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))
#   
#   # #all tree types within region
#   # plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="within 1000bp from causal pos significant")
#   # axis(side=1, at=1:length(hs_all), labels=hs_all)
#   # 
#   # lines(trueTreeResult$power_GWAS_region, col="orange2", lty=1, type = "b", pch=20)
#   # lines(trueTreeResult$power_REML_region, col="dodgerblue", lty=1, type = "b", pch=20)
#   # 
#   # lines(relateTreeResults$power_GWAS_region, col="orange2", lty=2, type = "b", pch=20)
#   # lines(relateTreeResults$power_REML_region, col="dodgerblue", lty=2, type = "b", pch=20)
#   
#   # legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))
# 
#   #only relate trees
#   plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
#   axis(side=1, at=1:length(hs_all), labels=hs_all)
#   
#   lines(relateTreeResults$power_GWAS, col="orange2", lty=2, type = "b", pch=20)
#   lines(relateTreeResults$power_REML_eGRM, col="dodgerblue", lty=2, type = "b", pch=20)
#   lines(relateTreeResults$power_REML_GRM, col="maroon2", lty=2, type = "b", pch=20)
#   # lines(power_results_highFreq$power_HE_SD, col="orange2")
#   # lines(power_results_highFreq$power_HE_CP, col="red")
#   
#   legend("bottomright", legend=c("GWAS", "REML eGRM", "REML GRM", "true trees", "relate trees"), bty='n', lty=c(1,1,1,1,2), col=c("orange2", "dodgerblue", "maroon2", "black", "black"))
#   
#   #only true trees
#   plot(0, type='n', ylim=c(0,1), xlim=c(1,length(hs_all)), xaxt='n', xlab="heritability", ylab="power", main="any distance from causal pos significant")
#   axis(side=1, at=1:length(hs_all), labels=hs_all)
#   
#   lines(trueTreeResult$power_GWAS, col="orange2", lty=1, type = "b", pch=20)
#   lines(trueTreeResult$power_REML_eGRM, col="dodgerblue", lty=1, type = "b", pch=20)
#   lines(trueTreeResult$power_REML_GRM, col="maroon2", lty=1, type = "b", pch=20)
#   
#   dev.off()
# }
# 
# plot_freq(freq=0.02, trueTreeResult=power_results_trueTrees_lowFreq, relateTreeResult=power_results_relateTrees_lowFreq)
# plot_freq(freq=0.2, trueTreeResult=power_results_trueTrees_highFreq, relateTreeResult=power_results_relateTrees_highFreq)

#--------------------------
# power files single variant
#--------------------------

setwd("/data/ARGWAS/power_sims/stdpopsim")

source("/home/linkv/git/argwas/R_scripts/2b_calculate_power_one_experiment.R")
library("pwr")
options(scipen = 100, digits = 4)
hs_all <- c(0.02) #, 0.04,  0.06,  0.08, 0.1     , 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
# hs_all <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
run_acat <- TRUE

for(tree_type in c("true_trees")){ #, "relate_trees", true_trees
  for(variant_freq in c("rareVariant", "commonVariant")){
    for(region_type in c("window_based")){
      for(ws_testing in c("10k")){ #, "20k", "50k" , "10k", ,"10k"
        for(ws_causal in  c("10k")){
          power_results_aH <- data.frame()
          for(hsquared in hs_all){
            folder=paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneVariant/",variant_freq, "/eGRM_GRM/", region_type, "/", ws_causal, "/" ,sep="")
            print(paste("analyzing folder", folder))
            if(run_acat){
              results_file <- paste(folder, "h", hsquared, "/power_results_acat.txt", sep='')
            } else {
              results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
            }
            print(paste(hsquared, file.exists(results_file)))
            # if(file.exists(results_file) == FALSE){
              if(ws_testing == ws_causal){
                pheno_file_dir <- folder
              } else {
                pheno_file_dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/" ,sep="")
              }
              power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type=tree_type, region_type=region_type, window_size_testing=ws_testing, window_size_causal=ws_causal, pheno_file_dir=pheno_file_dir, run_acat=run_acat)
            # }
            t <- read.table(results_file, header=TRUE)
            power_results_aH <- rbind(power_results_aH, t)
          }
        }
      }
    }
  }
}

#--------------------------
# power files allelic heterogeneity
#--------------------------

setwd("/data/ARGWAS/power_sims/stdpopsim")

source("/home/linkv/git/argwas/R_scripts/2b_calculate_power_one_experiment.R")
library("pwr")
options(scipen = 100, digits = 4)
hs_all <- c(0.02) #, 0.04,  0.06,  0.08, 0.1     , 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
# hs_all <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
run_acat <- TRUE

#read low freq true trees
for(propCausal in c(0.1, 0.2, 0.5, 0.8)){ #0.1, 
  for(tree_type in c("relate_trees")){ #, "relate_trees", true_trees
    for(region_type in c("window_based")){
      for(ws_testing in c("5k")){ #, "20k", "50k" , "10k", ,"10k"
        for(ws_causal in  c("5k")){
          power_results_aH <- data.frame()
          for(hsquared in hs_all){
            folder=paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/tested", ws_testing, "/propCausal", propCausal, "/" ,sep="")
            print(paste("analyzing folder", folder))
            if(run_acat){
              results_file <- paste(folder, "h", hsquared, "/power_results_acat.txt", sep='')
            } else {
              results_file <- paste(folder, "h", hsquared, "/power_results.txt", sep='')
            }
            print(paste(hsquared, file.exists(results_file)))
            # if(file.exists(results_file) == FALSE){
              if(ws_testing == ws_causal){
                pheno_file_dir <- folder
              } else {
                pheno_file_dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/tested", ws_causal, "/propCausal", propCausal, "/" ,sep="")
              }
              power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type=tree_type, region_type=region_type, window_size_testing=ws_testing, window_size_causal=ws_causal, pheno_file_dir=pheno_file_dir, run_acat=run_acat)
              # power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type="high_mut_trees")
              # print("finished creating results file")
            # }
            t <- read.table(results_file, header=TRUE)
            power_results_aH <- rbind(power_results_aH, t)
          }
        }
      }
    }
  }
}


# # with ACAT results
# source("/home/linkv/git/argwas/R_scripts/2b_calculate_power_one_experiment.R")
# library(ACAT)
# library("pwr")
# 
# tree_type <- "true_trees"
# region_type <- "window_based"
# ws_testing <- "10k"
# ws_causal <- "5k"
# hs_all <- c(0.02, 0.04,  0.06,  0.08, 0.1)
# propCausal <- 0.2
# 
# folder=paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/tested", ws_testing, "/propCausal", propCausal, "/" ,sep="")
# if(ws_testing == ws_causal){
#   pheno_file_dir <- folder
# } else {
#   pheno_file_dir <- paste("/data/ARGWAS/power_sims/stdpopsim/", tree_type, "/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/tested", ws_causal, "/propCausal", propCausal, "/" ,sep="")
# }
# for(hsquared in hs_all){
#   results_file <- paste(folder, "h", hsquared, "/power_results_acat.txt", sep='')
#   print(paste(hsquared, file.exists(results_file)))
#   # if(file.exists(results_file) == FALSE){
#     power_one_experiment(hsquared = hsquared, REPS = 200, folder=folder, tree_type=tree_type, region_type=region_type, window_size_testing=ws_testing, window_size_causal=ws_causal, pheno_file_dir=pheno_file_dir, run_acat=TRUE)
#   # }
# }