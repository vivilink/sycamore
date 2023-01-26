# sudo apt-get install build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev
# sudo apt install libfontconfig1-dev
# sudo apt install libtiff-dev
install.packages("systemfonts")
install.packages("textshaping")
install.packages("devtools")

library(devtools)
library(ACAT)
devtools::install_github("yaowuliu/ACAT")
setwd("/data/ARGWAS/power_sims/stdpopsim/true_trees/oneRegion/eGRM_GRM/window_based/5k/tested5k/propCausal0.2")

hs_all <- c(0.1)

for(hs in hs_all){
  t_variants <- read.table(paste("h", hs, "/rep1/power_sims_1_GWAS_variants_results.csv", sep=''), header=TRUE, sep=',')
  t_windows <- read.table(paste("h", hs, "/rep1/power_sims_1_eGRM_trees_REML_results.csv", sep=''), header=TRUE, sep=',')[which(!is.na(t_windows$V_e)),]
  
  for(r in 1:nrow(t_windows)){
    start_w <- t_windows$start[r]
    end_w <- t_windows$end[r]
    ps <- t_variants$p_value[t_variants$start >= start_w & t_variants$start < end_w]
    print(ACAT(ps))
  }
}

     