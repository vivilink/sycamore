library("psych")

# setwd("~/ARGWAS/hawaiian/all_chr_for_review/chr5")
# source("~/ARGWAS/argwas/R_scripts/functions.R")

setwd("~/postdoc_USC/AIM/correlations_p_values")
source("~/git/argwas/R_scripts/functions.R")
source("~/git/argwas/R_scripts/read_data.R")
num_windows = 15266
# region to plot
region_start <- 0
region_end <- 182045439
# region_start <- 172000000
# region_end <- 175000000
chrom <- 5

# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

plot_correlation <- function(eGRM, gwas, eGRM_name=expression("-log"[10]*"(local eGRM)"), GWAS_name){
  eGRM <- log10(eGRM)
  gwas[gwas==0] <- .Machine$double.xmin
  gwas <- log10(gwas)
  # eGRM <- eGRM[is.finite(gwas)]
  # gwas <- gwas[is.finite(gwas)]
  c <- cor(eGRM, gwas, method="spearman", use="complete.obs")
  print(paste("correlation for",GWAS_name,c))
  plot(x=eGRM, y=gwas , xlab=eGRM_name, ylab=GWAS_name, yaxt='n', bty='n', main=paste("Spearman correlation", round(c, 3))) #
  axis(side=2, las=2)
}


# remove encode regions
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions)
df_GRM_globalGRM <- remove_regions(df_results=df_GRM_globalGRM, regions=regions)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions)
df_GWAS_PC20 <- remove_regions_GWAS(df_results=df_GWAS_PC20, regions=regions)

# remove centromere regions
df_BLUP_res <- remove_regions(df_results=df_BLUP_res, regions=regions_centro)
df_PC100_egrm <- remove_regions(df_result=df_PC100_egrm, regions=regions_centro)
df_GRM_globalGRM <- remove_regions(df_results=df_GRM_globalGRM, regions=regions_centro)
df_GWAS_GRM <- remove_regions_GWAS(df_results=df_GWAS_GRM, regions=regions_centro)
df_GWAS_PC20 <- remove_regions_GWAS(df_results=df_GWAS_PC20, regions=regions_centro)


make_plots <- function(df_eGRM, df_GWAS, df_GRM, MAIN){
  df_eGRM <- df_eGRM[df_eGRM$start >= region_start & df_eGRM$start <= region_end,]
  df_eGRM <- remove_regions(df_results=df_eGRM, regions=regions_centro)
  df_GRM <- df_GRM[df_GRM$start >= region_start & df_GRM$start <= region_end,]
  df_GRM <- remove_regions(df_results=df_GRM, regions=regions_centro)
  df_GWAS <- df_GWAS[df_GWAS$BP >= region_start & df_GWAS$BP <= region_end,]
  df_GWAS <- df_GWAS[df_GWAS$BP >= region_start & df_GWAS$BP <= region_end,]
  
  df_eGRM <- remove_regions(df_result=df_eGRM, regions=regions)
  df_eGRM <- remove_regions(df_result=df_eGRM, regions=regions_centro)
  
  # df_GRM <- remove_regions(df_result=df_GRM, regions=regions)
  # df_GRM <- remove_regions(df_result=df_GRM, regions=regions_centro)
  
  df_GWAS <- remove_regions_GWAS(df_result=df_GWAS, regions=regions)
  df_GWAS <- remove_regions_GWAS(df_result=df_GWAS, regions=regions_centro)
  
  num_windows <- min(num_windows, length(df_eGRM$p_values))
  
  gwas_mean_arithmetic <- rep(NA, length=num_windows)
  gwas_mean_geometric <- rep(NA, length=num_windows)
  gwas_mean_harmonic <- rep(NA, length=num_windows)
  gwas_min <- rep(NA, length=num_windows)

  eGRM <- (df_eGRM$p_values)[1:num_windows]
  
  # for(w in 1:nrow(df_eGRM)){
  for(w in 1:num_windows){
      
    if(w %% 1000 == 0){
      print(paste("at window", w, "of", nrow(df_eGRM)))
    }
    window_start <- df_eGRM$start[w]
    window_end <- df_eGRM$end[w]
  
    p_GWAS_in_window <- (df_GWAS$P[which(df_GWAS$BP >= window_start & df_GWAS$BP < window_end)])
    if(length(p_GWAS_in_window) > 0){
      gwas_mean_arithmetic[w] <- mean(p_GWAS_in_window, na.rm=TRUE)
      gwas_mean_geometric[w] <- geometric.mean(p_GWAS_in_window, na.rm=TRUE)
      gwas_mean_harmonic[w] <- harmonic.mean(p_GWAS_in_window,na.rm=TRUE)
      gwas_min[w] <- min(p_GWAS_in_window, na.rm=TRUE)
      print(paste(w, min(p_GWAS_in_window, na.rm=TRUE)))
    }
  }

  pdf(paste("p_value_correlations", MAIN, ".pdf", sep=''), width=8, height=8)
  par(mfrow=c(2,2))
  plot_correlation(eGRM, gwas=gwas_mean_arithmetic, GWAS_name=expression("-log"[10]*"(GWAS arithmetic mean)"))
  plot_correlation(eGRM, gwas=gwas_mean_geometric, GWAS_name=expression("-log"[10]*"(GWAS geometric mean)"))
  plot_correlation(eGRM, gwas=gwas_mean_harmonic, GWAS_name=expression("-log"[10]*"(GWAS harmonic mean)"))
  plot_correlation(eGRM, gwas=gwas_min, GWAS_name=expression("-log"[10]*"(GWAS min)"))
  dev.off()
}



make_plots(df_eGRM=df_PC100_egrm, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_PC20, MAIN="egrmPC100_GRM")
make_plots(df_eGRM=df_BLUP_res, df_GWAS=df_GWAS_GRM, df_GRM=df_GRM_PC20, MAIN="residuals")

df_BLUP_res <- remove_regions(df_result=df_BLUP_res, regions=regions)
df_BLUP_res <- remove_regions(df_result=df_BLUP_res, regions=regions_centro)
df_GRM_globalGRM <- remove_regions(df_result=df_GRM_globalGRM, regions=regions)
df_GRM_globalGRM <- remove_regions(df_result=df_GRM_globalGRM, regions=regions_centro)
df_GRM_residuals <- remove_regions(df_result=df_GRM_residuals, regions=regions)
df_GRM_residuals <- remove_regions(df_result=df_GRM_residuals, regions=regions_centro)
eGRM <- log10(df_BLUP_res$p_values[1:num_windows])
GRM <- log10(df_GRM_residuals$p_values[1:num_windows])
c <- cor(eGRM, GRM, method="spearman", use="complete.obs")

pdf(paste("p_value_correlations_eGRM_GRMResiduals.pdf", sep=''), width=4, height=4)
plot(eGRM, GRM, ylab=expression("-log"[10]*"(local GRM)"), xlab=expression("-log"[10]*"(local eGRM)"),yaxt='n', bty='n', main=paste("Spearman correlation", round(c, 3))) #
axis(side=2, las=2)
dev.off()



