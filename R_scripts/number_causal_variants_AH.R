setwd("/data/ARGWAS/power_sims/stdpopsim")

hs_all <- c(0.02, 0.04,  0.06,  0.08, 0.1) #,      , 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005,
pc_all <- c(0.1, 0.2, 0.5, 0.8)
# hs_all <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005,
run_acat <- TRUE

#read low freq true trees
num_causal_table <- matrix(ncol=6, nrow=length(pc_all))
colnames(num_causal_table) <- c("Min", "Q25", "Median", "Mean", "Q75", "Max")
rownames(num_causal_table) <- paste("proportion causal",pc_all)

for(propCausal in pc_all){ #0.1, 0.2, 0.5, 0.8
  for(tree_type in c("relate_trees")){ #, "relate_trees", true_trees
    for(region_type in c("window_based")){
      for(ws_testing in c("5k")){ #, "20k", "50k" , "10k", ,"10k"
        for(ws_causal in  c("5k")){
          power_results_aH <- data.frame()
          for(hsquared in 0.02){ #all hsquared have same phenotypes
            num_causal <- numeric(length=200)
            for(rep in 1:200){
              pheno_file_dir <- paste("/data/ARGWAS/power_sims/stdpopsim/relate_trees/oneRegion/eGRM_GRM/", region_type, "/", ws_causal, "/tested", ws_causal, "/onlyUntyped/propCausal", propCausal, "/" ,sep="")
              df_pheno <- read.csv(paste(pheno_file_dir,"/h", hsquared, "/rep", rep,"/power_sims_", rep, "_pheno_causal_vars.csv", sep=''))
              # print(sum(df_pheno$causal == TRUE))
              num_causal[rep] <- sum(df_pheno$causal == TRUE)
            }
            print(paste("propCausal", propCausal))
            print(summary(num_causal))
            print(which(propCausal == pc_all))
            num_causal_table[which(propCausal == pc_all)] <- summary(num_causal)
          }
        }
      }
    }
  }
}

print(xtable(t, type = "latex"), file = paste("num_causal_variants_AH_causalWindow5kb.tex", sep=''), floating = FALSE)
