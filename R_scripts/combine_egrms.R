setwd("~/ARGWAS/hawaiian")
library("plinkFile")

ReadGRMBin=function(prefix, AllN=F, size=4){
  
  #part of script that is from gcta website
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }  else {N=readBin(NFile, n=1, what=numeric(0), size=size)}
  i=sapply(1:n, sum_i)
  
  #written by me, putting parts together
  diag=grm[i]
  off=grm[-i]
  m <- matrix(nrow=n, ncol=n)
  m[upper.tri(m, diag=FALSE)] <- off
  diag(m) <- diag
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  total_grm <- m
  
  return(list(diag=grm[i], off=grm[-i], id=id, N=N, total_grm=total_grm))

  close(BinFileName)
  close(NFileName)
  close(IDFileName)
}

files <- read.table("list_trees_global_egrm.txt")

# read first egrm to get dimensions
prefix <- strsplit(files$V1[1], split=".trees")[[1]]
prefix <- strsplit(prefix, split="/")[[1]][7]
grm_obj <- ReadGRMBin(prefix)


# prepare empty containers (need to name rows with individual IDs or else they are inferred wrongly)
m <- matrix(0, nrow(grm_obj$total_grm), ncol(grm_obj$total_grm))
total_N <- 0
sample_names <- read.table("relate.sample", header=T)
sample_names <- sample_names[-1,]
rownames(m) <- sample_names$ID_1
colnames(m) <- sample_names$ID_1

# add each partial egrm
for(f in files$V1){
	prefix <- strsplit(f, split=".trees")[[1]]
	prefix <- strsplit(prefix, split="/")[[1]][7]
	print(paste("adding grm",prefix))
	grm_obj <- ReadGRMBin(prefix)
	m <- m + grm_obj$total_grm * grm_obj$N
	total_N <- total_N + grm_obj$N
}

m <- m / total_N

print(paste("total_N", total_N))

#write global egrm to file
saveGRM(pfx="global", grm=m, vcm = total_N, fid = ".")





# #write global egrm to file
# BinFileName <- file("global.grm.bin", "wb")
# NFileName <- file("global.grm.N.bin", "wb")
# IDFileName <- file("global.grm.id", "w")
# 
# writeBin(as.vector(m[lower.tri(m, diag = TRUE)]), BinFileName)
# writeBin(total_N, NFileName)
# write(unlist(grm_obj$id), IDFileName)
# close(BinFileName)
# close(NFileName)
# close(IDFileName)

# grm_obj_read <- ReadGRMBin(prefix="global")
# print(sum(signif(grm_obj_read$total_grm,5) != signif(m,5))) 
# index <- which(signif(grm_obj_read$total_grm,4) != signif(m,4)) 
# 
# grm_obj_read$total_grm[index[1:10]]
# m[index[1:10]]
# 
# grm_obj_read$diag[1:10]
# diag(m)[1:10]
