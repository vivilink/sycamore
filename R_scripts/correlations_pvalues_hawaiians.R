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

df_GWAS <- read.table(paste("~/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr",chrom,"/plink.assoc.linear_chr", chrom, sep=''), header=TRUE)
df_GWAS <- df_GWAS[df_GWAS$TEST == "ADD",]
df_GWAS <- df_GWAS[df_GWAS$BP >= region_start & df_GWAS$BP <= region_end,]


# read REML 
df <- read.table(paste("~/ARGWAS/hawaiian/all_chr_for_review/chr", chrom, "/chr", chrom, "_all_chunks_eGRM_pca20_results.csv", sep=''), sep=',', header=TRUE) 
df <- df[!is.na(df$p_values),]
df <- df[df$start >= region_start & df$start <= region_end,]
df <- df[order(df$start, decreasing=FALSE),]

# remove encode regions
regions <- read.table("~/ARGWAS/argwas/R_scripts/encode_blacklist.bed", sep='\t', header=FALSE)
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == paste("chr", chrom, paste=''),]
df <- remove_regions(df_results=df, regions=regions)

# remove centromere
regions <- read.table("~/ARGWAS/argwas/R_scripts/centromeres.bed", sep='\t', header=FALSE)
colnames(regions) <- c("chr", "start", "end", "type")
regions <- regions[regions$chr == paste("chr", chrom, paste=''),]
df <- remove_regions(df_results=df, regions=regions)

# read REML GRM 
df_GRM <- read.table("~/ARGWAS/hawaiian/all_chr_for_review/chr",chrom,"/chr",chrom,"_all_chunks_GRM_pca20_results.csv", sep=',', header=TRUE) 
df_GRM <- df_GRM[!is.na(df$p_values),]
df_GRM <- df_GRM[df_GRM$start >= region_start & df_GRM$start <= region_end,]
df_GRM <- df_GRM[order(df_GRM$start, decreasing=FALSE),]

eGRM <- log10(df$p_values)
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

plot_correlation <- function(eGRM, gwas, XLAB){
  c <- cor(eGRM[!is.nan(gwas)], gwas[!is.nan(gwas)])
  plot(eGRM[!is.na(gwas)], gwas[!is.na(gwas)], main=c)
}

pdf("p_value_correlations.pdf", width=8, height=8)
par(mfrow=c(2,2))
plot_correlation(eGRM, gwas=gwas_mean_arithmetic, XLAB="arithmetic mean")
plot_correlation(eGRM, gwas=gwas_mean_geometric, XLAB="geometric mean")
plot_correlation(eGRM, gwas=gwas_mean_harmonic, XLAB="harmonic mean")
plot_correlation(eGRM, gwas=gwas_min, XLAB="min")
dev.off()











