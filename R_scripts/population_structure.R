setwd("/data/ARGWAS/experiments_population_split/association_tests/eGRM_GRM/true_trees/window_based/5k/")

plot_qq <- function(p_values, MAIN){
  unif <- runif(5000)
  qqplot(unif, p_values, xlim=c(0,1), ylim=c(0,1), main=MAIN, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', pch=20)
  axis(side=1, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1)) 
  axis(side=2, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1), las=2)
  title(ylab="p", line=2.8)
  title(xlab="Uniform(0,1)", line=2.2)
  abline(a=0, b=1)
}

pdf("p_value_qqplots.pdf", width=10, height=5)
par(mfrow=c(1,2))
n_reps <- 100
p_values <- c()

for(r in c(2:13,15:71,73:n_reps)){
  df <- read.csv(paste("no_strat_correction_highNe_splitTime10k/two_pops_",r,"_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}

plot_qq(p_values, MAIN="")

n_reps <- 100
p_values <- c()

for(r in 1:n_reps){
  df <- read.csv(paste("with_strat_correction_highNe_splitTime10k/two_pops_", r, "_phenoSeed_eGRMGlobalSeed_PCA_eGRM_trees_REML_results.csv", sep=''))
  p_values <- c(p_values, df$p_values)
}

plot_qq(p_values, MAIN="")

dev.off()


# 
# 
# pheno <- read.table("two_pops_98_eGRM_phenotypes.phen", header=FALSE)
# plot(density(pheno$V3[1:1000]))
# lines(density(pheno$V3[1001:2000]))
# mean(pheno$V3[1:1000])
# mean(pheno$V3[1001:2000])
# 
# 
# 
# pdf("mgrm_corrected.pdf", width=5, height=5)
# 
# n_reps <- 100
# p_values <- c()
# 
# for(r in c(1:20, 22:n_reps)){
#   df <- read.csv(paste("two_pops_", r, "_phenoSeed_eGRMGlobalSeed_mgrm_eGRM_trees_REML_results.csv", sep=''))
#   p_values <- c(p_values, df$p_values)
# }
# 
# qqplot(runif(1000), p_values, bty='n')
# abline(a=0, b=1)
# 
# dev.off()
# 

# 
# ReadGRMBin=function(prefix, AllN=F, size=4){
#   
#   #part of script that is from gcta website
#   sum_i=function(i){
#     return(sum(1:i))
#   }
#   BinFileName=paste(prefix,".grm.bin",sep="")
#   NFileName=paste(prefix,".grm.N.bin",sep="")
#   IDFileName=paste(prefix,".grm.id",sep="")
#   id = read.table(IDFileName)
#   n=dim(id)[1]
#   BinFile=file(BinFileName, "rb");
#   grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
#   NFile=file(NFileName, "rb");
#   if(AllN==T){
#     N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
#   }
#   else N=readBin(NFile, n=1, what=numeric(0), size=size)
#   i=sapply(1:n, sum_i)
#   
#   #written by me, putting parts together
#   diag=grm[i]
#   off=grm[-i]
#   m <- matrix(nrow=n, ncol=n)
#   m[upper.tri(m, diag=FALSE)] <- off
#   diag(m) <- diag
#   m[lower.tri(m)] <- t(m)[lower.tri(m)]
#   total_grm <- m
#   
#   return(list(diag=grm[i], off=grm[-i], id=id, N=N, total_grm=total_grm))
# }
# 
# 
# plot_hm <- function(egrm, title){
#   color = function(x)rev(heat.colors(x))
#   vbreaks=seq(-10,10, by=0.1)
#   # vbreaks=seq(range(egrm$total_grm)[1],range(egrm$total_grm)[2], by=0.1)
#   print(vbreaks)
#   png(paste(title,"_heatmap.png", sep=''), 1800, 1800)
#   print(range(egrm$total_grm))
#   h=heatmap(egrm$total_grm, main=title, col=color(length(vbreaks)-1), breaks=vbreaks) #, Rowv=NA, Colv=NA , 
#   dev.off()
#   return(h)
# }
# 
# setwd("/data/ARGWAS/experiments_population_split/association_tests/eGRM_GRM/true_trees/window_based/5k/with_strat_correction_highNe_splitTime1k")
# 
# globalGRM <- ReadGRMBin("two_pops_33")
# plot_hm(egrm=globalGRM, title="two_pops_33")
# 
# 




