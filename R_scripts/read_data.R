CLUSTER=FALSE

# colors
org_full <- "#E69F00"
org <- rgb(230, 159, 0, maxColorValue=255, alpha=100)
blu_full <- "#56B4E9"
blu <- rgb(86, 180, 233, maxColorValue=255, alpha=100)
pin_full <- "#CC79A7"
pin <- rgb(204, 121, 167, maxColorValue=255, alpha=50)

# read GWAS
if(CLUSTER==TRUE){
  df_GWAS_PC20 <- read.table("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/plink.assoc.linear_chr5", header=TRUE)
} else{
  df_GWAS_PC20 <- data.frame()
}
df_GWAS_PC20 <- df_GWAS_PC20[df_GWAS_PC20$TEST == "ADD",]

if(CLUSTER==TRUE){
  df_GWAS_GRM <- read.table("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/mlma_loco/MLMA_LOCO.loco.mlma", header=TRUE)
} else {
  df_GWAS_GRM <- read.table("~/postdoc_USC/AIM/correlations_p_values/MLMA_LOCO.loco.mlma", header=TRUE)
}
names <- colnames(df_GWAS_GRM)
names[which(names=="bp")] <- "BP"
names[which(names=="p")] <- "P"
colnames(df_GWAS_GRM) <- names

# read local GRM 
if(CLUSTER==TRUE){
  df_GRM_PC20 <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr5/only_PC_correction/chr5_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) 
} else {
  df_GRM_PC20 <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) 
}
df_GRM_PC20 <- df_GRM_PC20[!is.na(df_GRM_PC20$p_values),]
df_GRM_PC20 <- df_GRM_PC20[order(df_GRM_PC20$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_GRM_PC100 <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr5/only_PC_correction/pca100/chr5_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) 
} else {
  df_GRM_PC100 <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_GRM_pca100_results.csv", sep=',', header=TRUE) 
}
df_GRM_PC100 <- df_GRM_PC100[!is.na(df_GRM_PC100$p_values),]
df_GRM_PC100 <- df_GRM_PC100[order(df_GRM_PC100$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_GRM_PC100_16 <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr16/only_PC_correction/pca100/chr5_all_chunks_GRM_pca100_results.csv", sep=',', header=TRUE) 
} else {
  df_GRM_PC100_16 <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr16_all_chunks_GRM_pca100_results.csv", sep=',', header=TRUE) 
}
df_GRM_PC100_16 <- df_GRM_PC100_16[!is.na(df_GRM_PC100_16$p_values),]
df_GRM_PC100_16 <- df_GRM_PC100_16[order(df_GRM_PC100_16$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_GRM_globalGRM <- read.table("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr5/only_grm_correction/GRM/chr5_all_chunks_GRM_globalGRM_results.csv", sep=',', header=TRUE) 
} else {
  df_GRM_globalGRM <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_GRM_globalGRM_results.csv", sep=',', header=TRUE) 
}
df_GRM_globalGRM <- df_GRM_globalGRM[!is.na(df_GRM_globalGRM$p_values),]
df_GRM_globalGRM <- df_GRM_globalGRM[order(df_GRM_globalGRM$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_GRM_residuals <- read.table("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr5/only_grm_correction/GRM/chr5_all_chunks_GRM_globalGRM_results.csv", sep=',', header=TRUE) 
} else {
  df_GRM_residuals <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_GRM_globalGRMResiduals_results.csv", sep=',', header=TRUE) 
}
df_GRM_residuals <- df_GRM_residuals[!is.na(df_GRM_residuals$p_values),]
df_GRM_residuals <- df_GRM_residuals[order(df_GRM_residuals$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_GRM_residuals_16 <- read.table("/home1/linkv/ARGWAS/hawaiian/all_chr_for_review/chr16/", sep=',', header=TRUE) 
} else {
  df_GRM_residuals_16 <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr16_all_chunks_GRM_BLUP_residuals_correction_results.csv", sep=',', header=TRUE) 
}
df_GRM_residuals_16 <- df_GRM_residuals_16[!is.na(df_GRM_residuals_16$p_values),]
df_GRM_residuals_16 <- df_GRM_residuals_16[order(df_GRM_residuals_16$start, decreasing=FALSE),]


# read local eGRM
if(CLUSTER==TRUE){
  df_BLUP_res <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr5/BLUP/association_BLUP_residuals/chr5_all_chunks_eGRM_eGRMBLUP_residuals_results.csv", sep=',', header=TRUE) 
} else{
  df_BLUP_res <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_eGRM_eGRMBLUP_residuals_results.csv", sep=',', header=TRUE) 
}
df_BLUP_res <- df_BLUP_res[!is.na(df_BLUP_res$p_values),]
df_BLUP_res <- df_BLUP_res[order(df_BLUP_res$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_BLUP_res_16 <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr16/BLUP/association_BLUP_residuals/chr16_all_chunks_eGRM_BLUP_residuals_correction_results.csv", sep=',', header=TRUE) 
} else{
  df_BLUP_res_16 <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr16_all_chunks_eGRM_BLUP_residuals_correction_results.csv", sep=',', header=TRUE) 
}
df_BLUP_res_16 <- df_BLUP_res_16[!is.na(df_BLUP_res_16$p_values),]
df_BLUP_res_16 <- df_BLUP_res_16[order(df_BLUP_res_16$start, decreasing=FALSE),]

if(CLUSTER==TRUE){
  df_PC100_egrm <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr5/egrm_and_pca_correction/chr5_all_chunks_eGRM_pca100_egrm_and_pca_correction_results.csv", sep=',', header=TRUE)
} else{
  df_PC100_egrm <- read.table("~/postdoc_USC/AIM/correlations_p_values/chr5_all_chunks_eGRM_pca100_egrm_and_pca_correction_results.csv", sep=',', header=TRUE)
}
df_PC100_egrm <- df_PC100_egrm[!is.na(df_PC100_egrm$p_values),]
df_PC100_egrm <- df_PC100_egrm[order(df_PC100_egrm$start, decreasing=FALSE),]

# centromere regions
if(CLUSTER==TRUE){
  regions_centro <- read.table("~/ARGWAS/argwas/R_scripts/centromeres.bed", sep='\t', header=FALSE)
} else {
  regions_centro <- read.table("~/git/argwas/R_scripts/centromeres.bed", sep='\t', header=FALSE)
}
colnames(regions_centro) <- c("chr", "start", "end", "type")

# encode regions
if(CLUSTER==TRUE){
  regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
} else {
  regions <- read.table("~/git/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
}
colnames(regions) <- c("chr", "start", "end", "type")
