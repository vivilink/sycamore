CLUSTER <- FALSE

if(CLUSTER){
  setwd("~/ARGWAS/hawaiian/all_chr_for_review/chr5/")
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
region_start <- 172000000
region_end <- 175000000

# GWAS hits

rs373863828_causal <- 173108771 #causal variant  
rs12513649_proxy <- 173044949 #proxy variant

if(CLUSTER){
  other_GWAS <- read.table("/home1/linkv/ARGWAS/argwas/R_scripts/gwascat_hitsnearcrebrf.csv", sep=',', header=TRUE)
} else {
  other_GWAS <- read.table("~/git/argwas/R_scripts/gwascat_hitsnearcrebrf.csv", sep=',', header=TRUE)
}
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
  
  abline(h=-log10(5*(10^-8)), col=org, lty=2)
  abline(h=-log10(4.5*10^-7), col=blu, lty=2)
  abline(h=-log10(2.3*10^-7), col=pin, lty=2)

  # axes
  label_pos <- seq(region_start, region_end, 500000)
  axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
  title(ylab=expression("-log"[10]*"(p)"), line=2)
  title(xlab="genomic position [Mb]", line=2.2)
}

plot_association <- function(df, df_GWAS, df_GRM, method=""){
  plot(type='n', x=0, xaxt='n', las=2, xlab="", ylab="", ylim=c(0,8), xlim=c(region_start,region_end), bty='n')
  do_annotation(rs373863828_causal, rs12513649_proxy, other_GWAS)
  points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col=org, pch=20, lwd=0)
  points(x=df$start, y=-log10(df$p_values), col=blu, pch=20, lwd=0)
  points(x=df_GRM$start, y=-log10(df_GRM$p_values), col=pin, pch=20, lwd=0)

  index_min_pvalue <- which(df$p_values == min(df$p_values))
  print(paste("min pvalue",min(df$p_values)))
  print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828_causal - df$start[index_min_pvalue]) / 1000), "kb", "pvalue",-log10(df$p_values[index_min_pvalue])))

  # legend(legend=c("local eGRM", "GWAS"), pch=20, col=c(blu, org), horiz=FALSE, x="topright", box.lwd = 0, box.col = "white", bg = "white")
  legend(legend=c("local eGRM","local GRM", "GWAS"), horiz=FALSE, pch=20, col=c(blu,pin, org), x="topright", box.lwd = 0, box.col = "white", bg = "white")
#  legend(legend=c("GWAS"), pch=20, col=c(org), x="topleft", bty='n')
  
}

# read local eGRM
df_BLUP_res <- df_BLUP_res[df_BLUP_res$start >= region_start & df_BLUP_res$start <= region_end,]
df_PC100_egrm <- df_PC100_egrm[df_PC100_egrm$start >= region_start & df_PC100_egrm$start <= region_end,]

# read local GRM 
df_GRM_residuals <- df_GRM_residuals[df_GRM_residuals$start >= region_start & df_GRM_residuals$start <= region_end,]
df_GRM_PC100 <- df_GRM_PC100[df_GRM_PC100$start >= region_start & df_GRM_PC100$start <= region_end,]

# read GWAS
df_GWAS_PC20 <- df_GWAS_PC20[df_GWAS_PC20$BP >= region_start & df_GWAS_PC20$BP <= region_end,]
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$BP >= region_start & df_GWAS_GRM$BP <= region_end,]

# read encode regions
regions <- regions[regions$chr == "chr5",]
regions_centro <- regions_centro[regions_centro$chr == "chr5",]
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions)
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions_centro)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions_centro)
df_GWAS_GRM <- remove_regions_GWAS(df_result=df_GWAS_GRM, regions=regions)
df_GWAS_GRM <- remove_regions_GWAS(df_result=df_GWAS_GRM, regions=regions_centro)
df_GRM_residuals <- remove_regions(df_result=df_GRM_residuals, regions=regions)
df_GRM_residuals <- remove_regions(df_result=df_GRM_residuals, regions=regions_centro)
df_GRM_PC100 <- remove_regions(df_results=df_GRM_PC100, regions=regions)
df_GRM_PC100 <- remove_regions(df_results=df_GRM_PC100, regions=regions_centro)

#plot
pdf(paste("CREBRF.pdf", sep=''), width=8, height=4)
#png(paste("CREBRF_PC", num_PCs, "_GWAS.png", sep=''), width=8, height=4, units="in", res=1200)
  
par(mfrow=c(1,1))

plot_association(df=df_BLUP_res, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_residuals)
# plot_association(df=df_PC100_egrm, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_globalGRM)

dev.off()

print(paste("difference between cutoff and df_BLUP_res is", -log10(4.5*10^-7) - max(-log10(df_BLUP_res$p_values))))
print(paste("difference between cutoff and df_GWAS_GRM is", -log10(5*(10^-8)) - max(-log10(df_GWAS_GRM$P))))

