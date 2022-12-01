#scp linkv@RRI403-PC04.dts.usc.edu:/data/ARGWAS/power_sims/stdpopsim/high_mut_trees/oneTree/eGRM/h0.1/power_*txt .

# setwd("/data/ARGWAS/power_sims/stdpopsim/high_mut_trees/oneTree/eGRM_and_GRM/")
setwd("/data/ARGWAS/power_sims/stdpopsim/")
nreps=200

hs <- c(0.02,0.04,0.06,0.08, 0.1)
# hs <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
offset <- 0.08

pdf("power_aH_erorBars_window_based.pdf", height=5, width=10)
par(mfrow=c(1,2))
for(ws in c("5k", "10k")){  #, "20k", "50k"
  plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="power", xlab='local heritability', xaxt='n', main=ws, bty='n', las=2)
  axis(side=1, at=1:length(hs), labels=hs)
  for(h in hs){
    print(paste("h",h))
    x_pos <- which(hs == h)
    power_results <- read.table(paste("relate_trees/oneRegion/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    # power_results <- read.table(paste(ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    
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
    power_results <- read.table(paste("true_trees/oneRegion/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    
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
  legend(x="bottomright", legend=c("true trees / all variants","Relate / typed variants", "local eGRM", "local GRM", "GWAS"), pch=c(19, 1, 15, 15, 15), col=c("black","black","dodgerblue", "maroon2", "orange2"), bty='n')
}

dev.off()


colors <- c("dodgerblue", "maroon2") #, "orange2"
window_sizes <- c("5k","10k") #, "20k", "50k"

pdf("power_aH_erorBars_window_based_methods.pdf", height=5, width=20)
par(mfrow=c(1,4))
plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="power", xlab='local heritability', xaxt='n', main="eGRM")
for(ws in window_sizes){
  axis(side=1, at=1:length(hs), labels=hs)
  for(h in hs){
    x_pos <- which(hs == h)
    power_results <- read.table(paste("relate_trees/oneRegion/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)

    #REML eGRM
    points(y=power_results$power_REML_eGRM, x=x_pos+offset, pch=1, col=colors[which(window_sizes == ws)])
    power <- power_results$power_REML_eGRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=colors[which(window_sizes == ws)])
    print(power_results$power_REML_eGRM)
    
  }
}
legend(x="topleft", legend=window_sizes, pch=c(1, 1, 1), col=colors, bty='n')


plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="power", xlab='local heritability', xaxt='n', main="GWAS")
for(ws in window_sizes){
  axis(side=1, at=1:length(hs), labels=hs)
  print(paste("ws", ws))
  for(h in hs){
    x_pos <- which(hs == h)
    power_results <- read.table(paste("relate_trees/oneRegion/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    print(paste("h", h, power_results$power_GWAS))
    
    #GWAS
    points(y=power_results$power_GWAS, x=x_pos-offset, pch=1, col=colors[which(window_sizes == ws)])
    power <- power_results$power_GWAS
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=colors[which(window_sizes == ws)])	
    
  }
}
legend(x="topleft", legend=window_sizes, pch=c(1, 1, 1), col=colors, bty='n')


plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+offset), ylab="power", xlab='local heritability', xaxt='n', main="GRM")
for(ws in window_sizes){
  axis(side=1, at=1:length(hs), labels=hs)
  for(h in hs){
    x_pos <- which(hs == h)
    power_results <- read.table(paste("relate_trees/oneRegion/eGRM_GRM/window_based/", ws, "/h", h, "/power_results.txt", sep=''), header=TRUE)
    
    #REML GRM
    points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=1, col=colors[which(window_sizes == ws)])
    power <- power_results$power_REML_GRM
    std <- sqrt(power*(1-power)/nreps)
    segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=colors[which(window_sizes == ws)])
    
  }
}
legend(x="topleft", legend=window_sizes, pch=c(1, 1, 1), col=colors, bty='n')


dev.off()
