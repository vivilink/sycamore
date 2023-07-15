remove_regions <- function(df_results, regions){
  df_results$start <- as.numeric(df_results$start)
  for(r in 1:nrow(regions)){
    region_start <- regions$start[r]
    region_end <- regions$end[r]
    #start in region
    if(length(df_results$start[which(df_results$start >= region_start & df_results$start <= region_end)]) > 0){
      df_results <- df_results[-which(df_results$start >= region_start & df_results$start <= region_end),]
    }
    #end in region
    if(length(df_results$end[which(df_results$end >= region_start & df_results$end <= region_end)]) > 0){
      df_results <- df_results[-which(df_results$end >= region_start & df_results$end <= region_end),]
    }
    
  }
  return(df_results)
}


plot_qq <- function(p_values, MAIN){
  unif <- runif(5000)
  qqplot(unif, p_values, xlim=c(0,1), ylim=c(0,1), main=MAIN, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', pch=20)
  axis(side=1, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1)) 
  axis(side=2, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1), las=2)
  title(ylab="p", line=2.8)
  title(xlab="Uniform(0,1)", line=2.2)
  abline(a=0, b=1)
}

plot_qq_REML <- function(p_values, MAIN){
  unif <- runif(5000)
  qqplot(unif, p_values, xlim=c(0,1), ylim=c(0,1), main=MAIN, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', pch=20)
  axis(side=1, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1)) 
  axis(side=2, at=seq(0,1,by=0.2), labels = c(0,seq(0.2,0.8,by=0.2),1), las=2)
  title(ylab="p", line=2.8)
  title(xlab="Uniform(0,1)", line=2.2)
  segments(x0=0, y0=0, x1=0.5, y1=0.5)
  segments(x0=0.5, y0=0.5, x1=1, y1=0.5)
}

plot_qq_log <- function(p_values, MAIN){
  unif <- runif(5000)
  qqplot(log(unif), p_values, xlim=c(-10,0), ylim=c(-10,0), main=MAIN, xlab="", ylab="", bty='n', xaxt='n', yaxt='n', pch=20)
  axis(side=1, at=seq(-10,0,by=1), labels = seq(-10,0,by=1)) 
  axis(side=2, at=seq(-10,0,by=1), labels = seq(-10,0,by=1), las=2)
  title(ylab="p", line=2.8)
  title(xlab="Uniform(0,1)", line=2.2)
  abline(a=0, b=1)
}

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
