setwd("/data/ARGWAS/experiments_population_split/tests/msprime")
# setwd("/home/vivian/postdoc_USC/AIM/debugging/two_pops/heatmap/msprime")

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
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
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
}


plot_hm <- function(egrm, title){
  color = function(x)rev(heat.colors(x))
  vbreaks=seq(range(egrm$total_grm)[1],range(egrm$total_grm)[2], by=0.1)
  print(vbreaks)
  png(paste(title,"_heatmap.png", sep=''), 1800, 1800)
  print(range(egrm$total_grm))
  h=heatmap(egrm$total_grm, main=title) #, Rowv=NA, Colv=NA , col=color(length(vbreaks)-1), breaks=vbreaks
  dev.off()
  return(h)
}

# --------------------
# testing if GRM is being put together correctly
# --------------------

globalGRM <- ReadGRMBin("testing_GRM")
plot_hm(egrm=globalGRM, title="testing_GRM")

globaleGRM <- ReadGRMBin("testing_eGRM")
plot_hm(egrm=globaleGRM, title="testing_eGRM")

globaleGRM <- ReadGRMBin("testing")
plot_hm(egrm=globaleGRM, title="testing")

GRM_direct <- as.matrix(read.csv("testing_GRM_covariance_matrix.csv", sep=',', header=FALSE))
heatmap(as.matrix(GRM_direct))
sum(round(GRM_direct, 4) != round(globalGRM$whole_grm, 4))

eGRM_direct <- as.matrix(read.csv("testing_eGRM_covariance_matrix.csv", sep=',', header=FALSE))
heatmap(as.matrix(eGRM_direct))
sum(round(eGRM_direct, 4) != round(globaleGRM$whole_grm, 4))

# --------------------
# comparing actual matrices
#---------------------


globalGRM <- ReadGRMBin("two_pops_40")
h_global <- plot_hm(globalGRM, "global")
range(h$rowInd[1000:2000])

notworking <- ReadGRMBin("two_pops_20_phenoRep1_eGRMrep1_eGRM")
h_not_working <- plot_hm(notworking, "not_working")

working <- ReadGRMBin("two_pops_6_phenoRep1_eGRMrep1_eGRM")
h_working <- plot_hm(working, "working")





#splittime20k
notworking97 <- ReadGRMBin("splitTime20k/two_pops_97_phenoRep1_eGRMrep1_eGRM")
h_not_working97 <- plot_hm(notworking97, "splitTime20k/not_working_97")

notworking96 <- ReadGRMBin("splitTime20k/two_pops_96_phenoRep1_eGRMrep1_eGRM")
h_not_working96 <- plot_hm(notworking96, "splitTime20k/not_working_96")

notworking95 <- ReadGRMBin("splitTime20k/two_pops_95_phenoRep1_eGRMrep1_eGRM")
h_not_working95 <- plot_hm(notworking95, "splitTime20k/not_working_95")

working99 <- ReadGRMBin("splitTime20k/two_pops_99_phenoRep1_eGRMrep1_eGRM")
h_working_99 <- plot_hm(working99, "splitTime20k/working_99")

working9 <- ReadGRMBin("splitTime20k/two_pops_9_phenoRep1_eGRMrep1_eGRM")
h_working_9 <- plot_hm(working9, "splitTime20k/working_9")

working84 <- ReadGRMBin("splitTime20k/two_pops_84_phenoRep1_eGRMrep1_eGRM")
h_working_84 <- plot_hm(working84, "splitTime20k/working_84")

working94 <- ReadGRMBin("splitTime20k/two_pops_94_phenoRep1_eGRMrep1_eGRM")
h_working_94 <- plot_hm(working94, "splitTime20k/working_94")

globalGRM <- ReadGRMBin("splitTime20k/two_pops_99")
h_global <- plot_hm(globalGRM, "splitTime20k/global_two_pops_99")
range(h_global$rowInd[1000:2000])

globalGRM <- ReadGRMBin("splitTime20k/two_pops_97")
h_global <- plot_hm(globalGRM, "splitTime20k/global_two_pops_97")
range(h_global$rowInd[1000:2000])

globalGRM <- ReadGRMBin("splitTime20k/two_pops_1")
h_global <- plot_hm(globalGRM, "splitTime20k/global_two_pops_1")
range(h_global$rowInd[1000:2000])
#------------
# hamming distance
#------------
# install.packages("cultevo")
# library("cultevo")

haps <- read.table("/data/ARGWAS/experiments_population_split/tests/msprime/testing_haplotypes.txt", sep=' ', header=FALSE)
m_haps <- as.matrix((haps))
# m_hamm <- hammingdists(m_haps)
m_hamm <- dist(m_haps, method="manhattan")
dst <- data.matrix(m_hamm)
dim <- ncol(dst)
heatmap(dst)




