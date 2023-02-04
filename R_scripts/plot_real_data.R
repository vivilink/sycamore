setwd("/data/ARGWAS/hawaiians/CREBRF")
df <- read.table("chr5.part-04.CREBRF_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE)
df <- df[!is.na(df$p_values),]
df_GWAS <- read.table("plink_GWAS/plink.assoc.linear", header=TRUE)
df_GWAS <- df_GWAS[df_GWAS$TEST == "ADD",]

do_annotation <- function(){
  xmin <- 173056352
  xmax <- 173139284
  mygray <- col2rgb("gray")
  mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 80)
  polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)
  rs373863828_causal <- 173108771 #causal variant
  rs12513649_proxy <- 173044949 #proxy variant
  abline(v=rs373863828_causal,col="red")
  abline(v=rs12513649_proxy,col="black")
}



pdf("CREBRF.pdf", width=8, height=4)

par(mfrow=c(1,1))

# # REML
# plot(x=df$start, y=-log10(df$p_values), xaxt='n', las=2, xlab="window start [mb]", ylab=expression('log'[10]*'(p-value)'))
# label_pos <- seq(172000000, 173500000, 100000)
# axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
# do_annotation()
# index_min_pvalue <- which(df$p_values == min(df$p_values))
# print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828 - df$start[index_min_pvalue]) / 1000), "kb"))

# both
plot(x=df$start, y=-log10(df$p_values), xaxt='n', las=2, xlab="genomic position [mb]", ylab=expression('log'[10]*'(p-value)'), col="dodgerblue")
points(x=df_GWAS$BP, y=-log10(df_GWAS$P), col="orange2")
label_pos <- seq(172000000, 173500000, 100000)
axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
do_annotation()
index_min_pvalue <- which(df$p_values == min(df$p_values))
print(paste("distance rs373863828 and most significant REML hit:", round(abs(rs373863828 - df$start[index_min_pvalue]) / 1000), "kb"))
abline(h=-log10(5*(10^-8)), col="gray", lty=2)
legend(legend=c("local eGRM", "GWAS"), pch=1, col=c("dodgerblue", "orange2"), x="topleft", bty='n')

# # GWAS
# plot(x=df_GWAS$BP, y=-log10(df_GWAS$P), xaxt='n', las=2, xlab="window start [mb]", ylab=expression('log'[10]*'(p-value)'), col="orange2")
# label_pos <- seq(172000000, 173500000, 100000)
# axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
# do_annotation()
# index_min_pvalue <- which(df$p_values == min(df$p_values))
# print(paste("distance rs373863828 and most significant GWAS hit:", round(abs(rs373863828 - df_GWAS$BP[index_min_pvalue]) / 1000), "kb"))
# 





# add tree limits
# t_stats <- read.csv("chr5.part-04_trees_statistics.csv")
# segments(x0=df$start[index_min_pvalue], x1=rs373863828, y0=7, y1=7)
# for(e in t_stats$end){
#   abline(v=e, col="gray")
# }




dev.off()