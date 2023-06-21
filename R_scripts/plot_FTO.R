setwd("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr16")
source("/home1/linkv/ARGWAS/argwas/R_scripts/functions.R")

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

# region to plot
#53,703,963-54,114,467
region_start <- 52000000  #53000000
region_end <- 56000000



df_GWAS <- read.table("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr16/plink.assoc.linear_chr16", header=TRUE)
df_GWAS <- df_GWAS[df_GWAS$TEST == "ADD",]
df_GWAS <- df_GWAS[df_GWAS$BP >= region_start & df_GWAS$BP <= region_end,]


# GWAS hits

# rs373863828_causal <- 173108771 #causal variant  
# rs12513649_proxy <- 173044949 #proxy variant

other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/GWAScat_BMIweight_chr16.csv", sep=',', header=TRUE)
# other_GWAS <- other_GWAS[-which(other_GWAS$gwasCatalog.name == "rs12513649"),]

do_annotation <- function(other_GWAS){
  
  for(snp in other_GWAS$gwasCatalog.chromStart){
    abline(v=snp, col="snow2")
#	print(snp)
  }
  
  #FTO
  xmin <- 53703963
  xmax <- 54114467
  mygray <- col2rgb(pin)
  mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 50)
  polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)

  abline(h=-log10(5*(10^-8)), col=org, lty=2)
  abline(h=-log10(4.5*10^-7), col=blu, lty=2)

  # axes
  label_pos <- seq(region_start, region_end, 1000000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, num_PCs){
  pdf(paste("FTO_PC", num_PCs, ".pdf", sep=''), width=8, height=4)
  #png(paste("CREBRF_PC", num_PCs, "_GWAS.png", sep=''), width=8, height=4, units="in", res=1200)
  
  par(mfrow=c(1,1))
  
  # both
  plot(x=0, type='n', xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), ylim=c(0,8), bty='n') 
  do_annotation(other_GWAS)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20)
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20)
  points(x=df_GRM$start, y=-log10(df_GRM$p_values), col=pin, pch=20)

  index_min_pvalue <- which(df$p_values == min(df$p_values))
  # print(paste("min pvalue",min(df$p_values)))
  # print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))

  legend(legend=c("local eGRM","local GRM", "GWAS"), pch=20, col=c(blu,pin, org), x="topright", box.lwd = 0, box.col = "white", bg = "white")
  
  dev.off()
}

# -----------------------
# PC20
# -----------------------

df <- read.table("chr16_all_chunks_eGRM_pca20_results.csv", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
df <- df[!is.na(df$p_values),]
df <- df[df$start >= region_start & df$start <= region_end,]
df <- df[order(df$start, decreasing=FALSE),]

# remove encode regions
regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == "chr16",]
df <- remove_regions(df_results=df, regions=regions)

# read REML GRM 
df_GRM <- read.table("chr16_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
df_GRM <- df_GRM[!is.na(df$p_values),]
df_GRM <- df_GRM[df_GRM$start >= region_start & df_GRM$start <= region_end,]
df_GRM <- df_GRM[order(df_GRM$start, decreasing=FALSE),]

# remove encode regions
df_GRM <- remove_regions(df_results=df_GRM, regions=regions)


#plot
plot_association(df=df, num_PCs = 20)


# -----------------------
# PC12
# -----------------------

# read REML chunk 4
#df <- read.table("chr5.part-04.CREBRF_pca12_eGRM_trees_REML_results.csv", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
#df <- df[!is.na(df$p_values),]
#df <- df[df$start >= region_start & df$start <= region_end,]

# read REML chunk 5
#df2 <- read.table("chr5.part-05.CREBRF_pca12_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE) 
#df2 <- df2[!is.na(df2$p_values),]
#df2 <- df2[df2$start >= region_start & df2$start <= region_end,]

#df <- rbind(df, df2)

#plot
#plot_association(df=df, num_PCs = 12)

# -----------------------
# PC20
# -----------------------
# read REML chunk 4
#df <- read.table("chr5.part-04.CREBRF_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE) #cleaned just means the empty association tests (lines with ,,,,) are removed
#df <- df[!is.na(df$p_values),]
#df <- df[df$start >= region_start & df$start <= region_end,]

# read REML chunk 5
#df2 <- read.table("chr5.part-05.CREBRF_pca20_eGRM_trees_REML_results.csv", sep=',', header=TRUE) 
#df2 <- df2[!is.na(df2$p_values),]
#df2 <- df2[df2$start >= region_start & df2$start <= region_end,]

#df <- rbind(df, df2)

#plot
#plot_association(df=df, num_PCs = 20)


