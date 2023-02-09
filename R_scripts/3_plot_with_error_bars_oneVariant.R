#--------------------------
# plot power error bars
#--------------------------
setwd("/data/ARGWAS/power_sims/stdpopsim/")
nreps=200
# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

# hs <- c(0.02)
hs <- c(0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005,
offset <- 0.1
variant <- 
ws <- "10k"

pdf(paste("power_", variant, "_errorBars.pdf", sep=''), height=5, width=10)
par(mfrow=c(1,2))
# tiff(paste("power_", variant, "_errorBars_onlyTyped.tiff", sep=''), units="in", width=5, height=5, res=300)


for(variant in c("rareVariant", "commonVariant")){
  plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="", xlab='', xaxt='n', yaxt='n', main=variant, bty='n', las=2)
  axis(side=1, at=1:length(hs), labels=hs)
  axis(side=2, at=seq(0,1,by=0.2), labels=c(0,seq(0.2,0.8,by=0.2), 1), las=2)
  title(ylab=expression("power"), line=2.5)
  title(xlab="local heritability", line=2.2)
  
  for(h in hs){
    x_pos <- which(hs == h)
    
    power_results <- read.table(paste("relate_trees/oneVariant/",variant,"/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    
    #GWAS
    points(y=power_results$power_GWAS, x=x_pos-offset, pch=19, col=org)
    power <- power_results$power_GWAS
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=org)	
    
    #REML eGRM
    points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=19, col=blu)
    power <- power_results$power_REML_eGRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col=blu)
    
    #REML GRM
    points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=19, col=pin)
    power <- power_results$power_REML_GRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=pin)
    
    
    #true
    power_results <- read.table(paste("true_trees/oneVariant/", variant, "/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    
    # REML GWAS
    points(y=power_results$power_GWAS, x=x_pos-offset, pch=1, col=org)
    power <- power_results$power_GWAS
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=org)
    
    # REML eGRM
    points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=1, col=blu)
    power <- power_results$power_REML_eGRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col=blu)
    
    # REML GRM
    points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=1, col=pin)
    power <- power_results$power_REML_GRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=pin)
    
  }
}


# legend(x="bottomright", legend=c("local eGRM", "local GRM", "GWAS"), pch=c(19,19,19), col=c(blu, pin, org), bty='n')
legend(x="bottomright", legend=c("true trees / all variants","Relate / typed variants", "local eGRM", "local GRM", "GWAS"), pch=c(1, 19, 15, 15, 15), col=c("gray","gray",blu, pin, org), bty='n')

dev.off()
