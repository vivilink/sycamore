setwd("/data/ARGWAS/hawaiians/CREBRF")

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

# region to plot
region_start <- 172000000
region_end <- 175000000



df_GWAS <- read.table("plink_GWAS/plink.assoc.linear", header=TRUE)
df_GWAS <- df_GWAS[df_GWAS$TEST == "ADD",]
df_GWAS <- df_GWAS[df_GWAS$BP >= region_start & df_GWAS$BP <= region_end,]


# GWAS hits

rs373863828_causal <- 173108771 #causal variant  
rs12513649_proxy <- 173044949 #proxy variant

other_GWAS <- read.table("/home/linkv/git/argwas/R_scripts/gwascat_hitsnearcrebrf.csv", sep=',', header=TRUE)
other_GWAS <- other_GWAS[-which(other_GWAS$gwasCatalog.name == "rs12513649"),]

do_annotation <- function(rs373863828_causal, rs12513649_proxy, other_GWAS){
  
  for(snp in other_GWAS$position_GR38){
    abline(v=snp, col="gray")
  }
  
  #CREBRF
  xmin <- 173056352
  xmax <- 173139284
  mygray <- col2rgb(pin)
  mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 50)
  polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)

  abline(v=rs373863828_causal,col=pin)
  abline(v=rs12513649_proxy,col="gray")
  

}

plot_association <- function(df, num_PCs){
  pdf(paste("CREBRF_PC", num_PCs, ".pdf", sep=''), width=8, height=4)
  
  par(mfrow=c(1,1))
  
  # # REML
  # plot(x=df$start, y=-log10(df$p_values), xaxt='n', las=2, xlab="window start [mb]", ylab=expression('log'[10]*'(p-value)'))
  # label_pos <- seq(172000000, 173500000, 100000)
  # axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  # do_annotation()
  # index_min_pvalue <- which(df$p_values == min(df$p_values))
  # print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828 - df$start[index_min_pvalue]) / 1000), "kb"))
  
  # both
  plot(x=df$start, y=-log10(df$p_values), xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), bty='n')
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [mb]", line=2.2)
  
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20)
  label_pos <- seq(region_start, region_end, 500000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  do_annotation(rs373863828_causal, rs12513649_proxy, other_GWAS)
  index_min_pvalue <- which(df$p_values == min(df$p_values))
  print(paste("min pvalue",min(df$p_values)))
  print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))
  abline(h=-log10(5*(10^-8)), col="gray", lty=2)
  legend(legend=c("local eGRM", "GWAS"), pch=20, col=c(blu, org), x="topleft", bty='n')
  
  dev.off()
}



# -----------------------
# PC10
# -----------------------

# read REML chunk 4
df <- read.table("chr5.part-04.CREBRF_pca10_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
df <- df[!is.na(df$p_values),]
df <- df[df$start >= region_start & df$start <= region_end,]

# read REML chunk 5
df2 <- read.table("chr5.part-05.CREBRF_pca10_limit_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE)
df2 <- df2[!is.na(df2$p_values),]
df2 <- df2[df2$start >= region_start & df2$start <= region_end,]

df3 <- read.table("chr5.part-05.CREBRF_pca10_limit3_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE)
df3 <- df3[!is.na(df3$p_values),]
df3 <- df3[df3$start >= region_start & df3$start <= region_end,]

# combine chunks
df <- rbind(df,df2, df3)

#plot
plot_association(df=df, num_PCs = 10)


# -----------------------
# PC12
# -----------------------

# read REML chunk 4
df <- read.table("chr5.part-04.CREBRF_pca12_eGRM_trees_REML_results.csv", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
df <- df[!is.na(df$p_values),]
df <- df[df$start >= region_start & df$start <= region_end,]

#plot
plot_association(df=df, num_PCs = 12)

# -----------------------
# PC20
# -----------------------

