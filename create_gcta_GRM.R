library("plinkFile")
grm <- read.table("/data/ARGWAS/argwas/GRM_covariance.txt", header=F)
saveGRM(pfx="GRM_covariance", grm, vcm = NULL, fid = ".")

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