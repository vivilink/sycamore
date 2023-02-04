#--------------------------
# plot power error bars
#--------------------------
setwd("/data/ARGWAS/power_sims/stdpopsim/")
nreps=200

# hs <- c(0.02)
hs <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005,
offset <- 0.1
variant <- "commonVariant"
ws <- "10k"

pdf(paste("power_", variant, "_errorBars.pdf", sep=''), height=5, width=5)
# tiff(paste("power_", variant, "_errorBars_onlyTyped.tiff", sep=''), units="in", width=5, height=5, res=300)

plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="power", xlab='local heritability', xaxt='n', main=ws, bty='n', las=2)
axis(side=1, at=1:length(hs), labels=hs)

for(h in hs){
  x_pos <- which(hs == h)
  
  power_results <- read.table(paste("relate_trees/oneVariant/",variant,"/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
  
  #GWAS
  points(y=power_results$power_GWAS, x=x_pos-offset, pch=1, col="orange2")
  power <- power_results$power_GWAS
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos-offset, y0=power - std, y1=power + std, col="orange2")	
  
  #REML eGRM
  points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=1, col="dodgerblue")
  power <- power_results$power_REML_eGRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col="dodgerblue")

  #REML GRM
  points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=1, col="maroon2")
  power <- power_results$power_REML_GRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos+offset, y0=power - std, y1=power + std, col="maroon2")
  
  
  #true
  power_results <- read.table(paste("true_trees/oneVariant/", variant, "/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)

  # REML GWAS
  points(y=power_results$power_GWAS, x=x_pos-offset, pch=19, col="orange2")
  power <- power_results$power_GWAS
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos-offset, y0=power - std, y1=power + std, col="orange2")

  # REML eGRM
  points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=19, col="dodgerblue")
  power <- power_results$power_REML_eGRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col="dodgerblue")

  # REML GRM
  points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=19, col="maroon2")
  power <- power_results$power_REML_GRM
  std <- sqrt(power*(1-power)/nreps)
  segments(x0=x_pos+offset, y0=power - std, y1=power + std, col="maroon2")
  
}

# legend(x="topleft", legend=c("local eGRM", "local GRM", "GWAS"), pch=c(15, 15, 15), col=c("dodgerblue", "maroon2", "orange2"), bty='n')
legend(x="topleft", legend=c("true trees / all variants","Relate / typed variants", "local eGRM", "local GRM", "GWAS"), pch=c(19, 1, 15, 15, 15), col=c("black","black","dodgerblue", "maroon2", "orange2"), bty='n')

dev.off()
