#! /usr/bin/Rscript
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.0/")

library("plinkFile")
args = commandArgs(trailingOnly=TRUE)

grm <- read.table(paste(args[1], "_GRM_covariance.txt", sep=''), header=F)
saveGRM(pfx=paste(args[1], "_GRM_covariance", sep=''), grm, vcm = NULL, fid = ".")

# m <-matrix(ncol=500, nrow=500)
# for(v in 1:length(grm)){
#   m[,v] <- grm[[v]]
# }


# 
# 
# m <- as.matrix(grm)
# m[10,100]
# m[100,10]
# 
# isSymmetric(m)
# 
# 
# pheno <- read.table()
