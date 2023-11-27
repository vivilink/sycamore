# setwd("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr5")
# source("/home1/linkv/ARGWAS/argwas/R_scripts/functions.R")

setwd("/data/ARGWAS/hawaiians/association_all_chr/BLUP/local_eGRM")
source("~/git/sycamore/R_scripts/functions.R")
CLUSTER  <- FALSE

do_annotation <- function(other_GWAS){
  # 
  # for(snp in other_GWAS$gwasCatalog.chromStart){
  #   abline(v=snp, col="snow2")
  # }
  
  abline(h=-log10(5*(10^-8)), col=org, lty=2)
  abline(h=-log10(4.5*10^-7), col=blu, lty=2)
  
  # axes
  label_pos <- seq(region_start, region_end, 20000000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, df_GWAS, method="", CHROM){
  plot(type='n', x=100, xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), ylim=c(0,8.5), bty='n', main=CHROM) #ylim=c(0,max(-log10(df$p_values)))
  do_annotation(other_GWAS)
  points(x=df_GWAS$bp, y=-log10(df_GWAS$p), col=org, pch=20, lwd=0)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20, lwd=0)

  legend(legend=c("local eGRM", "GWAS"), pch=20, col=c(blu, org), x="topright", box.lwd=0, box.col = "white", bg = "white") #box.lwd=0, box.col = "white", bg = "white", , horiz = TRUE
}

# GWAS
df_GWAS <- read.table("MLMA_LOCO.loco.mlma", header=TRUE)

# centromere regions
if(CLUSTER==TRUE){
  regions_centro <- read.table("~/ARGWAS/sycamore/R_scripts/centromeres.bed", sep='\t', header=FALSE)
} else {
  regions_centro <- read.table("~/git/sycamore/R_scripts/centromeres.bed", sep='\t', header=FALSE)
}
colnames(regions_centro) <- c("chr", "start", "end", "type")

# encode regions
if(CLUSTER==TRUE){
  regions <- read.table("~/ARGWAS/sycamore/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
} else {
  regions <- read.table("~/git/sycamore/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
}
colnames(regions) <- c("chr", "start", "end", "type")

png(paste("hawaiians_BMI_all_chromosomes_residuals.png", sep=''), width=5, height=50, units="in", res=400)
par(mfrow=c(22,1))

for(CHROM in seq(1,22,1)){
  
  df_GWAS_chr <- df_GWAS[df_GWAS$Chr == CHROM,]

  if(CLUSTER==TRUE){
    df_BLUP_res <- read.table("", sep=',', header=TRUE) 
  } else{
    df_BLUP_res <- read.table(paste("chr", CHROM, "_association_results.csv", sep=''), sep=',', header=TRUE) 
  }
  df_BLUP_res <- df_BLUP_res[!is.na(df_BLUP_res$p_values),]
  df_BLUP_res <- df_BLUP_res[order(df_BLUP_res$start, decreasing=FALSE),]
  
  org <- "#E69F00"
  blu <- "#56B4E9"
  pin <- "#CC79A7"
  
  # region to plot
  region_start <- 0
  region_end <- 250000000
  
  # other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/GWAScat_BMIweight_chr5.csv", sep=',', header=TRUE)
  # other_GWAS <- read.table("~/git/argwas/R_scripts/GWAScat_BMIweight_chr5.csv", sep=',', header=TRUE)
  # other_GWAS <- other_GWAS[-which(other_GWAS$gwasCatalog.name == "rs12513649"),]
  
  # read local eGRM
  df_BLUP_res <- df_BLUP_res[df_BLUP_res$start >= region_start & df_BLUP_res$start <= region_end,]
 
  # remove encode regions
  regions_chr <- regions[regions$chr == paste("chr", CHROM, sep=''),]
  df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions_chr)

  # remove centromere regions
  regions_centro_chr <- regions_centro[regions_centro$chr == paste("chr", CHROM, sep=''),]
  df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions_centro_chr)

  # t <- table(as.factor(df$start))
  # print(paste("this position is present more than once", names(t)[as.integer(t) > 1]))
  
  #plot
  plot_association(df=df_BLUP_res, df_GWAS=df_GWAS_chr, CHROM=CHROM)
  
}
dev.off()

