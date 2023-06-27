setwd("~/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/association_with_GRM_correction")
source("~/ARGWAS/argwas/R_scripts/functions.R")

nreps <- 50

p_corrected <- c()
p_NOTcorrected <- c()

for(r in c(1:50)){
	result <- paste("CREBRF_test_region_withCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep='')
	if(file.exists(result)){
		print(r)
		t_corrected <- read.table(paste("CREBRF_test_region_withCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep=''), sep=',', header=TRUE)

		p_corrected <- c(p_corrected, t_corrected$p_values)
	#t_NOTcorrected <- read.table(paste("association_noCorrection/CREBRF_test_region_noCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep=''), sep=',', header=TRUE)
	#p_NOTcorrected <- c(p_NOTcorrected, t_NOTcorrected$p_values)
	}
}


pdf("stratification_hawaiians_p_values_GRM_correction.pdf", height=4, width=8)
par(mfrow=c(1,2))
plot_qq(p_corrected, MAIN="")
#plot_qq(p_NOTcorrected, MAIN="")
dev.off()

pdf("stratification_hawaiians_p_values_GRM_correction_log.pdf", height=4, width=8)
par(mfrow=c(1,2))
 unif <- runif(5000)
qqplot(log10(unif), log10(p_corrected), pch=20, ylim=c(-14,0), xlim=c(-14,0))
#qqplot(log10(unif), log10(p_NOTcorrected), main="", pch=20, ylim=c(-14,0), xlim=c(-14,0))
dev.off()



pdf("stratification_hawaiians_GRM_correction_individual_replicates.pdf", height=8*4, width=4*4)
par(mfrow=c(8,4))
for(r in c(1:50)){
#or(r in c(1, 3:6, 8:10, 12:13, 17:19, 21:22, 24:29, 35:43)){
	result <- paste("CREBRF_test_region_withCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep='')
	if(file.exists(result)){
		print(r)
		t_corrected <- read.table(paste("CREBRF_test_region_withCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep=''), sep=',', header=TRUE)

		p_corrected <- t_corrected$p_values

		#p_NOTcorrected <- read.table(paste("association_noCorrection/CREBRF_test_region_noCorrection_", r, "_eGRM_trees_GCTA_REML_results.csv", sep=''), sep=',', header=TRUE)$p_values
	
		plot_qq(p_corrected, MAIN=paste("rep", r))
#		#plot_qq(p_NOTcorrected, MAIN="")
	}
}
dev.off()



pdf("phenotype_distributions_replicates.pdf", height=8*4, width=4*4)
par(mfrow=c(8,4))
for(r in c(1:50)){
#or(r in c(1, 3:6, 8:10, 12:13, 17:19, 21:22, 24:29, 35:43)){

	phen <- read.table(paste("phenotypes_", r, ".phen", sep=''), sep=' ', header=FALSE)[,3]
	plot(density(phen), main=paste("rep",r, "mean", round(mean(phen), 2)))	

}
dev.off()
