library("plinkFile")
args = commandArgs(trailingOnly=TRUE)

grm <- read.table(paste("/data/ARGWAS/argwas/GRM_covariance_", args[1], ".txt", sep=''), header=F)
saveGRM(pfx=paste("GRM_covariance_", args[1], sep=''), grm, vcm = NULL, fid = ".")

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