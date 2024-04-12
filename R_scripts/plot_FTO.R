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

library("grDevices")

# colors
#org <- rgb(red=230, green=159, blue=0, max = 255, alpha=50)
org <- rgb(red=230, green=159, blue=0, max=255)
blu <- "#56B4E9"
pin <- "#CC79A7"

# region to plot
#53,703,963-54,114,467
region_start <- 52000000  #53000000
region_end <- 56000000


# GWAS hits

# rs373863828_causal <- 173108771 #causal variant  
# rs12513649_proxy <- 173044949 #proxy variant
if(CLUSTER==TRUE){
  other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/GWAScat_BMIweight_chr16.csv", sep=',', header=TRUE)
} else {
  other_GWAS <- read.table("~/git/argwas/R_scripts/GWAScat_BMIweight_chr16.csv", sep=',', header=TRUE)
}
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
  abline(h=-log10(2.3*10^-7), col=pin, lty=2)
  
  # axes
  label_pos <- seq(region_start, region_end, 1000000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p value)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, df_GRM, df_GWAS){
  #png(paste("CREBRF_PC", num_PCs, "_GWAS.png", sep=''), width=8, height=4, units="in", res=1200)
  
  par(mfrow=c(1,1))
  
  # both
  plot(x=0, type='n', xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), ylim=c(0,8), bty='n') 
  do_annotation(other_GWAS)
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20, lwd=0)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20, lwd=0)
  points(x=df_GRM$start, y=-log10(df_GRM$p_values), col=pin, pch=20, lwd=0)

  index_min_pvalue <- which(df$p_values == min(df$p_values))
  # print(paste("min pvalue",min(df$p_values)))
  # print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))

  legend(legend=c("local eGRM","local GRM", "GWAS"), pch=20, col=c(blu,pin, org), x="topright", box.lwd = 0, box.col = "white", bg = "white")
  
}

df_BLUP_res_16 <- df_BLUP_res_16[df_BLUP_res_16$start >= region_start & df_BLUP_res_16$start <= region_end,]
df_GRM_PC100_16 <- df_GRM_PC100_16[df_GRM_PC100_16$start >= region_start & df_GRM_PC100_16$start <= region_end,]
df_GRM_GRM_16 <- df_GRM_PC100_16[df_GRM_PC100_16$start >= region_start & df_GRM_PC100_16$start <= region_end,]
df_GRM_residuals_16 <- df_GRM_residuals_16[df_GRM_residuals_16$start >= region_start & df_GRM_residuals_16$start <= region_end,]

df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$Chr==16,]
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$BP >= region_start & df_GWAS_GRM$BP <= region_end,]

# remove regions
regions <- regions[regions$chr == 16,]
regions_centro <- regions_centro[regions_centro$chr == 16,]
df_BLUP_res_16 <- remove_regions(df_results=df_BLUP_res_16, regions=regions)
df_BLUP_res_16 <- remove_regions(df_results=df_BLUP_res_16, regions=regions_centro)
df_GRM_PC100_16 <- remove_regions(df_results=df_GRM_PC100_16, regions=regions)
df_GRM_PC100_16 <- remove_regions(df_results=df_GRM_PC100_16, regions=regions_centro)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions_centro)
df_GRM_residuals_16 <- remove_regions(df_results=df_GRM_residuals_16, regions=regions)
df_GRM_residuals_16 <- remove_regions(df_results=df_GRM_residuals_16, regions=regions_centro)

#plot
pdf(paste("FTO.pdf", sep=''), width=8, height=4)
plot_association(df=df_BLUP_res_16, df_GRM=df_GRM_residuals_16, df_GWAS=df_GWAS_GRM)

dev.off()

print(paste("difference between cutoff and df_BLUP_res is", -log10(4.5*10^-7) - max(-log10(df_BLUP_res_16$p_values))))
print(paste("difference between cutoff and df_GWAS_GRM is", -log10(5*(10^-8)) - max(-log10(df_GWAS_GRM$P))))
print(paste("difference between cutoff and df_GRM_residuals_16 is", -log10(2.3*10^-7) - max(-log10(df_GRM_residuals_16$p_values))))

