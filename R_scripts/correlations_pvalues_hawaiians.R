library("psych")

setwd("~/ARGWAS/hawaiian/all_chr_for_review/chr5")
chrom <- 5
source("~/ARGWAS/argwas/R_scripts/functions.R")

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"


# region to plot
region_start <- 0
region_end <- 182045439

df_GWAS_PC20 <- read.table(paste("~/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr",chrom,"/plink.assoc.linear_chr", chrom, sep=''), header=TRUE)
df_GWAS_PC20 <- df_GWAS_PC20[df_GWAS_PC20$TEST == "ADD",]
df_GWAS_PC20 <- df_GWAS_PC20[df_GWAS_PC20$BP >= region_start & df_GWAS_PC20$BP <= region_end,]

df_GWAS_GRM <- read.table("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/mlma_loco/MLMA_LOCO.loco.mlma", header=TRUE
)
df_GWAS_GRM <- df_GWAS_GRM[df_GWAS_GRM$bp >= region_start & df_GWAS_GRM$bp <= region_end,]
names <- colnames(df_GWAS_GRM)
names[which(names=="bp")] <- "BP"
names[which(names=="p")] <- "P"
colnames(df_GWAS_GRM) <- names


# read REML 
df_PC20 <- read.table(paste("~/ARGWAS/hawaiian/all_chr_for_review/chr", chrom, "/chr", chrom, "_all_chunks_eGRM_pca20_results.csv", sep=''), sep=',', header=TRUE) 
df_PC20 <- df_PC20[!is.na(df_PC20$p_values),]
df_PC20 <- df_PC20[df_PC20$start >= region_start & df_PC20$start <= region_end,]
df_PC20 <- df_PC20[order(df_PC20$start, decreasing=FALSE),]

df_egrm_pc100 <- read.table(paste("~/ARGWAS/hawaiian/all_chr_for_review/chr5/egrm_and_pca_correction/pc100/chr5_all_chunks_eGRM_pca100_egrm_and_pca_correction_results.csv", sep=''), header=TRUE)
df_egrm_pc100 <- df_egrm_pc100[!is.na(df_egrm_pc100$p_values),]
df_egrm_pc100 <- df_egrm_pc100[df_egrm_pc100$start >= region_start & df_egrm_pc100$start <= region_end,]
df_egrm_pc100 <- df_egrm_pc100[order(df_egrm_pc100$start, decreasing=FALSE),]

# remove encode regions
regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == paste("chr", chrom, paste=''),]
df_PC20 <- remove_regions(df_results=df_PC20, regions=regions)
df_egrm_pc100 <- remove_regions(df_results=df_egrm_pc100, regions=regions)

# remove centromere
regions <- read.table("~/ARGWAS/argwas/R_scripts/centromeres.bed", sep='\t', header=FALSE)
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == paste("chr", chrom, paste=''),]
df_PC20 <- remove_regions(df_results=df_PC20, regions=regions)
df_egrm_pc100 <- remove_regions(df_results=df_egrm_pc100, regions=regions)

# read REML GRM 
df_GRM <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr",chrom,"/chr",chrom,"_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) 
df_GRM <- df_GRM[!is.na(df$p_values),]
df_GRM <- df_GRM[df_GRM$start >= region_start & df_GRM$start <= region_end,]
df_GRM <- df_GRM[order(df_GRM$start, decreasing=FALSE),]

eGRM <- log10(df$p_values)

make_plots <- function(df_eGRM, df_GWAS, MAIN){
  gwas_mean_arithmetic <- numeric(length=nrow(df))
  gwas_mean_geometric <- numeric(length=nrow(df))
  gwas_mean_harmonic <- numeric(length=nrow(df))
  gwas_min <- numeric(length=nrow(df))

  for(w in 1:nrow(df)){
    window_start <- df$start[w]
    window_end <- df$end[w]
  
    p_GWAS_in_window <- (df_GWAS$P[which(df_GWAS$BP >= window_start & df_GWAS$BP < window_end)])
    gwas_mean_arithmetic[w] <- mean(p_GWAS_in_window, na.rm=TRUE)
    gwas_mean_geometric[w] <- geometric.mean(p_GWAS_in_window, na.rm=TRUE)
    gwas_mean_harmonic[w] <- harmonic.mean(p_GWAS_in_window,na.rm=TRUE)
    gwas_min[w] <- min(p_GWAS_in_window, na.rm=TRUE)
  }

  pdf(paste("p_value_correlations", MAIN, ".pdf", sep=''), width=8, height=8)
  par(mfrow=c(2,2))
  plot_correlation(eGRM, gwas=gwas_mean_arithmetic, XLAB="arithmetic mean")
  plot_correlation(eGRM, gwas=gwas_mean_geometric, XLAB="geometric mean")
  plot_correlation(eGRM, gwas=gwas_mean_harmonic, XLAB="harmonic mean")
  plot_correlation(eGRM, gwas=gwas_min, XLAB="min")
  dev.off()
}

plot_correlation <- function(eGRM, gwas, XLAB){
  c <- cor(eGRM[!is.nan(gwas)], gwas[!is.nan(gwas)])
  plot(eGRM[!is.na(gwas)], gwas[!is.na(gwas)], main=c)
}


make_plots(df_eGRM=df_egrm_PC100, df_GWAS=df_GWAS_GRM, MAIN="egrmPC100_GRM")






