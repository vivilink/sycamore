setwd("/data/ARGWAS/hawaiians/correlation_loadings_BMI")

# read eigenvec data
eigenvec <- read.table("chr5.part-05.CREBRF_pca50_eGRM.eigenvec", header=FALSE)[,-1]
colnames(eigenvec) <- c("ind", paste("PC", 1:50, sep=''))

# read phenotype data
pheno <- read.table("chr5.part-05.CREBRF_pca50_phenotypes.phen", header=FALSE)[,-1]
colnames(pheno) <- c("ind", "BMI")

# are they in same order?
eigenvec <- eigenvec[eigenvec$ind %in% pheno$ind,]
if(sum(eigenvec$ind != pheno$ind) > 0){
  stop("individuals don't match")
}

results <- matrix(nrow=50, ncol=3)
rownames(results) <- paste("PC", 1:50, sep='')
colnames(results) <- c("log10(p)", "R-squared", "slope")

for(i in 2:51){
  fit <- lm(pheno$BMI ~ eigenvec[,i])
  results[i-1, 1] <- log10(summary(fit)$coefficients[2,4])
  results[i-1, 2] <- (summary(fit)$r.squared)
  results[i-1, 3] <- summary(fit)$coefficients[2,1]
}

pdf("Supp_fig_corr_eigenvec_BMI.pdf", width=7, height=3)

par(mfrow=c(1,3))

plot(x=1:50, y=results[,1], main="", bty='n', xaxt='n', las=2, ylab="", xlab="", pch=20)
# points(x=12, y=results[12,1], col="red", pch=20)
axis(side=1, labels=TRUE)
title(ylab=expression("log"[10]*"(p)"), line=2.5)
title(xlab="PC", line=2.3)

plot(x=1:50, y=results[,2], main="", bty='n', xaxt='n', las=2, ylab="", xlab="", pch=20)
# points(x=12, y=results[12,2], col="red", pch=20)
axis(side=1, labels=TRUE)
title(ylab="R-squared", line=3.2)
title(xlab="PC", line=2.3)

plot(x=1:50, y=results[,3], main="", bty='n', xaxt='n', las=2, ylab="", xlab="", pch=20)
# points(x=12, y=results[12,3], col="red", pch=20)
axis(side=1, labels=TRUE)
title(ylab="slope", line=2.5)
title(xlab="PC", line=2.3)

dev.off()

print(xtable(results, type = "latex", digits=4), file = paste("CREBRF_correlation_PCA_loadings_BMI.tex", sep=''), floating = FALSE)

