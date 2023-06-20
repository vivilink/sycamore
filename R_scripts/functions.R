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