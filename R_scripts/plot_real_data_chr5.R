# setwd("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr5")
# source("/home1/linkv/ARGWAS/argwas/R_scripts/functions.R")

setwd("~/postdoc_USC/AIM/correlations_p_values")
source("~/git/argwas/R_scripts/functions.R")
source("~/git/argwas/R_scripts/read_data.R")

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

# region to plot
region_start <- 0
region_end <- 182045439

# GWAS hits

rs373863828_causal <- 173108771 #causal variant  
rs12513649_proxy <- 173044949 #proxy variant

# other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/GWAScat_BMIweight_chr5.csv", sep=',', header=TRUE)
other_GWAS <- read.table("~/git/argwas/R_scripts/GWAScat_BMIweight_chr5.csv", sep=',', header=TRUE)
other_GWAS <- other_GWAS[-which(other_GWAS$gwasCatalog.name == "rs12513649"),]

do_annotation <- function(rs373863828_causal, rs12513649_proxy, other_GWAS){
  
  for(snp in other_GWAS$gwasCatalog.chromStart){
    abline(v=snp, col="snow2")
  }
  
  # CREBRF
  xmin <- 173056352
  xmax <- 173139284
  mygray <- col2rgb(pin)
  mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 50)
  polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)

  abline(v=rs373863828_causal,col=pin)
  abline(v=rs12513649_proxy,col="gray")
  

  abline(h=-log10(5*(10^-8)), col=org, lty=2)
  abline(h=-log10(4.5*10^-7), col=blu, lty=2)

  # axes
  label_pos <- seq(region_start, region_end, 20000000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, df_GWAS, df_GRM, method=""){
  plot(type='n', x=100, xaxt='n', las=2, xlab="", ylab="", col=blu, pch=20, xlim=c(region_start,region_end), ylim=c(0,max(-log10(df$p_values))), bty='n')
  do_annotation(rs373863828_causal, rs12513649_proxy, other_GWAS)
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20, bty='n')
#  points(x=df_GRM$start, y=-log10(df_GRM$p_values), col=pin, pch=20)

  index_min_pvalue <- which(df$p_values == min(df$p_values))
  print(paste("min pvalue",min(df$p_values)))
  print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))

#  legend(legend=c("local eGRM", "GWAS", "local GRM"), pch=20, col=c(blu, org, pin), x="topleft", bty='n')
  legend(legend=c("local eGRM", "GWAS"), pch=20, col=c(blu, org), x="topleft", bty='n')
  
}

# -----------------------
# PC20
# -----------------------
# read local eGRM
df_BLUP_res <- df_BLUP_res[df_BLUP_res$start >= region_start & df_BLUP_res$start <= region_end,]
df_PC100_egrm <- df_PC100_egrm[df_PC100_egrm$start >= region_start & df_PC100_egrm$start <= region_end,]

# read encode regions
if(CLUSTER==TRUE){
  regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
} else {
  regions <- read.table("~/git/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
}
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == "chr5",]
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions)

# read local GRM 
df_GRM_PC20 <- df_GRM_PC20[df_GRM_PC20$start >= region_start & df_GRM_PC20$start <= region_end,]

# read GWAS
df_GWAS_PC20 <- df_GWAS_PC20[df_GWAS_PC20$BP >= region_start & df_GWAS_PC20$BP <= region_end,]
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$BP >= region_start & df_GWAS_GRM$BP <= region_end,]

# remove encode regions
if(CLUSTER==TRUE){
  regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
} else {
  regions <- read.table("~/git/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
}
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == "chr5",]
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions)
df_GRM_PC20 <- remove_regions(df_results=df_GRM_PC20, regions=regions)

# t <- table(as.factor(df$start))
# print(paste("this position is present more than once", names(t)[as.integer(t) > 1]))

#plot
# pdf(paste("hawaiians_BMI_chr5.pdf", sep=''), width=8, height=4)
png(paste("hawaiians_BMI_chr5.png", sep=''), width=8, height=4, units="in", res=600)

par(mfrow=c(1,1))
plot_association(df=df_BLUP_res, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_PC20)
# plot_association(df=df_PC100_egrm, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_PC20)
dev.off()



