source("functions.R")
library("MASS")
setwd("~/ARGWAS/hawaiian/stratification_control")
args <- commandArgs(trailingOnly=TRUE)
rep <- as.numeric(args[1])
set.seed(rep)

individuals <- read.table("/home1/linkv/ARGWAS/hawaiian/relate.sample", header=TRUE, sep=' ')[-1,1]
n_inds <- length(individuals)


locoGRM_obj <- ReadGRMBin("/home1/linkv/ARGWAS/hawaiian/global_grm/global")
locoGRM <- locoGRM_obj$total_grm


v_locoGRM  <- mvrnorm(n = 1, mu = rep(0, n_inds), Sigma = locoGRM, tol = 1e-6, empirical = TRUE, EISPACK = FALSE)

v_identity <- mvrnorm(n = 1, mu = rep(0, n_inds), Sigma = diag(n_inds), tol = 1e-6, empirical = TRUE, EISPACK = FALSE)


phenotypes <- v_locoGRM + v_identity

df <- cbind(0, individuals, phenotypes)
write.table(df, sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE, file=paste("phenotypes_", rep, ".phen", sep=''))



