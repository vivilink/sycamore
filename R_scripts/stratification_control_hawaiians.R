setwd("/home1/linkv/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/")
source("~/ARGWAS/argwas/R_scripts/functions.R")

nreps <- 50
unif <- runif(5000)

plot_ps_eGRM <- function(directory, MAIN, LOG=FALSE, ending="_eGRM_trees_GCTA_REML_results.csv"){
	print(paste("collecting data for",MAIN))
	p_corrected <- c()
	n_reps_found <- 0
	for(r in c(1:50)){
		#print(r)
		result <- paste(directory,"", r, ending, sep='')
		if(file.exists(result)){
			n_reps_found <- n_reps_found + 1
			t_corrected <- read.table(result, sep=',', header=TRUE)
			p_corrected <- c(p_corrected, t_corrected$p_values)
			#print(p_corrected)
		} else {
			print(paste("rep", r, "files does not exist", result))
		}


		if(n_reps_found == 30){
			break
		}
	}
	if(LOG==FALSE){
		plot_qq(p_corrected, MAIN=paste(MAIN, "log"))
	} else {
		plot_qq_log(log(p_corrected), MAIN=paste(MAIN, "log"))
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

		if(n_reps_found == 30){
			break
		}
	}
	
	plot_qq(p_corrected, MAIN=paste(MAIN))

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

		if(n_reps_found == 30){
			break
		}
	}
	
	plot_qq(p_corrected, MAIN=paste(MAIN))

}


pdf("stratification_experiments.pdf", width=12, height=10)

m <- matrix(1:30, ncol=6, byrow = TRUE)
l <- layout(m)


plot_ps_eGRM("association_with_eGRM_correction/CREBRF_test_region_withCorrection_", "global eGRM", LOG=FALSE)
plot_ps_eGRM("association_with_GRM_correction/CREBRF_test_region_withCorrection_", "global GRM", LOG=FALSE)
plot_ps_eGRM("association_with_100PC/CREBRF_test_region_withCorrection_", "100 PCs", LOG=FALSE)
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC20_","20 PC + global eGRM", LOG=FALSE)
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC100_","100 PC + global eGRM", LOG=FALSE)
plot_ps_eGRM("BLUP_residuals/global_eGRM/association_local_eGRM_residuals/CREBRF_test_region_withCorrection_", "global eGRM residuals", LOG=FALSE)


plot_ps_eGRM("association_with_eGRM_correction/CREBRF_test_region_withCorrection_", "global eGRM", LOG=TRUE)
plot_ps_eGRM("association_with_GRM_correction/CREBRF_test_region_withCorrection_", "global GRM", LOG=TRUE)
plot_ps_eGRM("association_with_100PC/CREBRF_test_region_withCorrection_", "100 PCs", LOG=TRUE)
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC20_","20 PC + global eGRM", LOG=TRUE)
plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC100_","100 PC + global eGRM", LOG=TRUE)
plot_ps_eGRM("BLUP_residuals/global_eGRM/association_local_eGRM_residuals/CREBRF_test_region_withCorrection_", "global eGRM residuals", LOG=TRUE)

plot_ps_eGRM("/home1/linkv/ARGWAS/hawaiian/stratification_control/multivariate_normal/association_with_grm_correction/CREBRF_test_region_withCorrection_", "MVN")
plot.new()
plot.new()
plot.new()
plot.new()
plot.new()

plot_ps_GWAS_GRM("GWAS_eGRM/crebrf_", "GWAS global eGRM", ending=".mlma")
plot_ps_GWAS_GRM("GWAS_GRM/crebrf_","GWAS global GRM", ending=".loco.mlma")
plot_ps_GWAS("GWAS_100PC/crebrf_region_", "GWAS 100 PCs", ending=".PHENO1.glm.linear")
plot_ps_GWAS("GWAS_20PC/crebrf_region_", "GWAS 20 PCs", ending=".PHENO1.glm.linear")
plot.new()
plot.new()

plot_ps_eGRM("~/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/association_with_GRM_correction/CREBRF_test_region_withCorrection_", "GRM global GRM", ending="_GRM_trees_GCTA_REML_results.csv")
plot_ps_eGRM("/home1/linkv/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/BLUP_residuals/global_GRM/associations_residuals/CREBRF_test_region_withCorrection_", "global GRM residuals", ending="_GRM_trees_GCTA_REML_results.csv")
plot.new()
plot.new()
plot.new()
plot.new()

dev.off()

#plot_ps_eGRM("association_with_100PC/CREBRF_test_region_withCorrection_", "100 PCs")
#plot_ps_eGRM("association_with_eGRM_correction/CREBRF_test_region_withCorrection_", "global eGRM", LOG=TRUE)
#plot_ps_eGRM("association_with_GRM_correction/CREBRF_test_region_withCorrection_", "global GRM")
#plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC20_","20 PC + global eGRM")
#plot_ps_eGRM("association_with_PC_eGRM_correction/CREBRF_test_region_withCorrection_PC100_","100 PC + global eGRM", LOG=TRUE)
#plot_ps_eGRM("BLUP_residuals/global_eGRM/association_local_eGRM_residuals/CREBRF_test_region_withCorrection_", "global eGRM residuals", LOG=TRUE)
#plot_ps_eGRM("/home1/linkv/ARGWAS/hawaiian/stratification_control/multivariate_normal/association_with_grm_correction/CREBRF_test_region_withCorrection_", "MVN")
#plot_ps_GWAS("GWAS_20PC/crebrf_region_", "GWAS 20 PCs", ending=".PHENO1.glm.linear")
#plot_ps_GWAS("GWAS_100PC/crebrf_region_", "GWAS 100 PCs", ending=".PHENO1.glm.linear")
#plot_ps_GWAS_GRM("GWAS_eGRM/crebrf_", "GWAS global eGRM", ending=".mlma")
#plot_ps_GWAS_GRM("GWAS_GRM/crebrf_","GWAS global GRM", ending=".loco.mlma")
#plot_ps_GWAS("GWAS_20PC/crebrf_region_", "GWAS 20 PCs", ending=".PHENO1.glm.linear")
#plot_ps_eGRM("~/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/association_with_GRM_correction/CREBRF_test_region_withCorrection_", "GRM global GRM", ending="_GRM_trees_GCTA_REML_results.csv")
#plot_ps_eGRM("/home1/linkv/ARGWAS/hawaiian/stratification_control/variants1k_heritability0.45/BLUP_residuals/global_GRM/associations_residuals/CREBRF_test_region_withCorrection_", "global GRM residuals", ending="_GRM_trees_GCTA_REML_results.csv")


#dev.off()



