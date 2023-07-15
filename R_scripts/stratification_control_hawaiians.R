setwd("/home1/linkv/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/")
source("~/ARGWAS/argwas/R_scripts/functions.R")

nreps <- 50
unif <- runif(5000)

plot_ps_eGRM <- function(directory, MAIN, LOG=FALSE){
	print(paste("collecting data for",MAIN))
	p_corrected <- c()
	n_reps_found <- 0
	for(r in c(1:50)){
		#print(r)
		result <- paste(directory,"", r, "_eGRM_trees_GCTA_REML_results.csv", sep='')
		if(file.exists(result)){
			n_reps_found <- n_reps_found + 1
			t_corrected <- read.table(result, sep=',', header=TRUE)
			p_corrected <- c(p_corrected, t_corrected$p_values)
			#print(p_corrected)
		} else {
			print(paste("rep", r, "files does not exist", result))
		}
	}
	plot_qq_REML(p_corrected, MAIN=paste(MAIN, n_reps_found))
	if(LOG==TRUE){
		plot_qq_log(log(p_corrected), MAIN=paste(MAIN, "log", n_reps_found))
	}

}

plot_ps_GWAS <- function(directory, MAIN, ending){
	p_corrected <- c()
	n_reps_found <- 0
	print(paste("collecting data for", MAIN))
	for(r in 1:50){
		#print(r)
		results <- paste(directory, "",r,ending, sep='')
		if(file.exists(results)){
			n_reps_found <- n_reps_found + 1
			t_corrected <- read.table(results, header=TRUE, comment.char="")
			p_corrected <- c(p_corrected, t_corrected$P)
		}else{
			print(paste("rep", r, "file does not exist", results))
		}
	}
	
	plot_qq(p_corrected, MAIN=paste(MAIN, n_reps_found))

}

plot_ps_GWAS_GRM <- function(directory, MAIN, ending){
	p_corrected <- c()
	n_reps_found <- 0
	print(paste("collecting data for", MAIN))
	for(r in 1:50){
		#print(r)
		results <- paste(directory, "",r,ending, sep='')
		if(file.exists(results)){
			n_reps_found <- n_reps_found + 1
			t_corrected <- read.table(results, header=TRUE, comment.char="")
			t_corrected <- t_corrected[t_corrected$Chr==5,]
			t_corrected <- t_corrected[which(t_corrected$bp >= 172750000 & t_corrected$bp <= 173000000),]
			p_corrected <- c(p_corrected, t_corrected$p)
			print(str(t_corrected))
		}else{
			print(paste("rep", r, "file does not exist", results))
		}
	}
	
	plot_qq(p_corrected, MAIN=paste(MAIN, n_reps_found))

}

pdf("stratification_experiments.pdf", width=20, height=4)
par(mfrow=c(2,10))
plot_ps_eGRM("association_with_100PC/CREBRF_test_region_withCorrection_", "100 PCs")
plot_ps_eGRM("association_with_eGRM_correction/CREBRF_test_region_withCorrection_", "global eGRM", LOG=TRUE)
plot_ps_eGRM("association_with_GRM_correction/CREBRF_test_region_withCorrection_", "global GRM")
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC20_","20 PC + global eGRM")
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC100_","100 PC + global eGRM", LOG=TRUE)
plot_ps_eGRM("BLUP_residuals/association_local_eGRM_residuals/CREBRF_test_region_withCorrection_", "BLUP residuals", LOG=TRUE)
plot_ps_eGRM("/home1/linkv/ARGWAS/hawaiian/stratification_control/multivariate_normal/association_with_grm_correction/CREBRF_test_region_withCorrection_", "MVN")
plot_ps_GWAS("GWAS_20PC/crebrf_region_", "GWAS 20 PCs", ending=".PHENO1.glm.linear")
plot_ps_GWAS("GWAS_100PC/crebrf_region_", "GWAS 100 PCs", ending=".PHENO1.glm.linear")
plot_ps_GWAS_GRM("GWAS_eGRM/crebrf_", "GWAS global eGRM", ending=".mlma")
plot_ps_GWAS_GRM("GWAS_GRM/crebrf_","GWAS global GRM", ending=".loco.mlma")
plot_ps_GWAS("GWAS_20PC/crebrf_region_", "GWAS 20 PCs", ending=".PHENO1.glm.linear")



dev.off()



