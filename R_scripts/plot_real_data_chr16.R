CLUSTER <- FALSE

if(CLUSTER==TRUE){
  setwd("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr16")
  source("/home1/linkv/ARGWAS/argwas/R_scripts/functions.R")
  source("/home1/linkv/ARGWAS/argwas/R_scripts/read_data.R")
} else {
  setwd("~/postdoc_USC/AIM/correlations_p_values")
  source("~/git/argwas/R_scripts/functions.R")
  source("~/git/argwas/R_scripts/read_data.R")
}

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

# region to plot
region_start <- 0
region_end <- 96330374

if(CLUSTER){
  other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/GWAScat_BMIweight_chr16.csv", sep=',', header=TRUE)
} else {
  other_GWAS <- read.table("~/git/argwas/R_scripts/GWAScat_BMIweight_chr16.csv", sep=',', header=TRUE)
}

do_annotation <- function(other_GWAS){
  for(snp in other_GWAS$gwasCatalog.chromStart){
    abline(v=snp, col="snow2")
    #	print(snp)
  }
  
  #FTO
  xmin <- 53703963
  xmax <- 54114467
  mygray <- col2rgb(pin)
  mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 100)
  # polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)
  
  abline(h=-log10(5*(10^-8)), col=org, lty=2)
  abline(h=-log10(4.5*10^-7), col=blu, lty=2)
  abline(h=-log10(2.3*10^-7), col=pin, lty=2)
  
  # axes
  label_pos <- seq(region_start, region_end, 10000000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, df_GRM, df_GWAS){
  pdf(paste("hawaiians_BMI_chr16.pdf", sep=''), width=8, height=4)
  #png(paste("CREBRF_PC", num_PCs, "_GWAS.png", sep=''), width=8, height=4, units="in", res=1200)
  
  par(mfrow=c(1,1))
  
  plot(type='n', x=100, xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), ylim=c(0,8.5), bty='n') #ylim=c(0,max(-log10(df$p_values)))
  do_annotation(other_GWAS)
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20, lwd=0)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20, lwd=0)
  points(x=df_GRM$start, y=-log10(df_GRM$p_values), col=pin, pch=20)

  index_min_pvalue <- which(df$p_values == min(df$p_values))
  print(paste("min pvalue",min(df$p_values)))
  #print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))

  legend(legend=c("local eGRM", "GWAS", "local GRM"), pch=20, col=c(blu, org, pin), x="topleft", bty='n', horiz = TRUE)
# legend(legend=c("local eGRM", "GWAS"), pch=20, col=c(blu, org), x="topleft", bty='n', horiz=TRUE) #, box.lwd=0, box.col = "white", bg = "white"
  
  dev.off()
}
df_BLUP_res_16 <- df_BLUP_res_16[df_BLUP_res_16$start >= region_start & df_BLUP_res_16$start <= region_end,]
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$Chr==16,]
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$BP >= region_start & df_GWAS_GRM$BP <= region_end,]
df_GRM_PC100_16 <- df_GRM_PC100_16[df_GRM_PC100_16$start >= region_start & df_GRM_PC100_16$start <= region_end,]
df_GRM_residuals_16 <- df_GRM_residuals_16[df_GRM_residuals_16$start >= region_start & df_GRM_residuals_16$start <= region_end,]

# remove regions
regions <- regions[regions$chr == "chr16",]
regions_centro <- regions_centro[regions_centro$chr == "chr16",]
df_BLUP_res_16 <- remove_regions(df_results=df_BLUP_res_16, regions=regions)
df_BLUP_res_16 <- remove_regions(df_results=df_BLUP_res_16, regions=regions_centro)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions_centro)
df_GRM_PC100_16 <- remove_regions(df_results=df_GRM_PC100_16, regions=regions)
df_GRM_PC100_16 <- remove_regions(df_results=df_GRM_PC100_16, regions=regions_centro)
df_GRM_residuals_16 <- remove_regions(df_results=df_GRM_residuals_16, regions=regions)
df_GRM_residuals_16 <- remove_regions(df_results=df_GRM_residuals_16, regions=regions_centro)

#plot
pdf(paste("chr16.pdf", sep=''), width=8, height=4)
plot_association(df=df_BLUP_res_16, df_GRM=df_GRM_residuals_16, df_GWAS=df_GWAS_GRM)
dev.off()


