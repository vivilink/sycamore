library("plinkFile")
grm <- read.table("/data/ARGWAS/argwas/GRM_covariance.txt", header=F)
saveGRM(pfx="GRM_covariance", grm, vcm = NULL, fid = ".")

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