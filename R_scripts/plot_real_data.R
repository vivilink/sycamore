setwd("/data/ARGWAS/hawaiians/CREBRF")
df <- read.table("chr5.part-04.CREBRF_eGRM_trees_REML_results.csv_cleaned", sep=',', header=TRUE)
df <- df[!is.na(df$p_values),]

t_stats <- read.csv("chr5.part-04_trees_statistics.csv")


pdf("CREBRF.pdf", width=8, height=4)
plot(x=df$start, y=-log10(df$p_values), xaxt='n', las=2, xlab="window start [mb]", ylab=expression('log'[10]*'(p-value)'))
label_pos <- seq(172000000, 173500000, 100000)
axis(1, at=label_pos, labels=label_pos / 1000000, las=1)
xmin <- 173056352
xmax <- 173139284
mygray <- col2rgb("gray")
mygray <- rgb(mygray[1], mygray[2], mygray[3], max = 255,  alpha = 80)
polygon(c(xmin,xmin, xmax, xmax), c(-100,100,100,-100), col = mygray, border=NA)

rs373863828 <- 173108771
abline(v=173108771,col="red")

index_min_pvalue <- which(df$p_values == min(df$p_values))
print(paste("distance rs373863828 and most significant hit:", round(abs(rs373863828 - df$start[index_min_pvalue]) / 1000), "kb"))
# segments(x0=df$start[index_min_pvalue], x1=rs373863828, y0=7, y1=7)

# for(e in t_stats$end){
#   abline(v=e, col="gray")
# }

dev.off()