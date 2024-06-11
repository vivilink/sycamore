setwd("/data/ARGWAS/power_sims/stdpopsim/")
nreps=100
run_acat <- NA
run_arg_needle <- NA
allowTyped <- "allowTyped"   #allowTyped onlyUntyped
# colors
org <- "#E69F00"
blu <- "#56B4E9"
pin <- "#CC79A7"

hs <- c(0.02,0.04,0.06,0.08, 0.1) #
# hs <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1) #, 0.07, 0.04, 0.0025, , 0.2 0.001, 0.0001, 0.0002, 0.0005, 
offset <- 0.08

for(propCausal in c(0.5)){ #0.1,0.2,0.5,
pdf(paste("power_binary_errorBars_window_based_propCausal", propCausal, ".pdf", sep=''), height=5, width=5)

par(mfrow=c(1,1))
  for(ws_testing in c("5k")){  # ,"10k"
    for(ws_causal in c("5k")){
      plot(0, type='n', ylim=c(0, 1), xlim=c(min(1)-offset, length(hs)+2*offset), ylab="", xlab='', yaxt='n', xaxt='n', main=ws_testing, bty='n', las=2)
      axis(side=1, at=1:length(hs), labels=hs)
      axis(side=2, at=seq(0,1,by=0.2), labels=c(0,seq(0.2,0.8,by=0.2), 1), las=2)
      
      title(ylab=expression("power"), line=2.5)
      title(xlab="local heritability", line=2.2)
      
      for(h in hs){
        print(paste("h",h))
        x_pos <- which(hs == h)
        
      #  # relate trees
      #   power_results <- read.table(paste("relate_trees/oneRegion/eGRM_GRM/window_based/", ws_causal,"/tested",ws_testing, "/", allowTyped, "/propCausal", propCausal, "/h", h, "/power_results.txt", sep=''), header=TRUE)
      # 
      #   # #GWAS
      #   # points(y=power_results$power_GWAS, x=x_pos-offset, pch=19, col=org)
      #   # power <- power_results$power_GWAS
      #   # std <- sqrt(power*(1-power)/nreps)
      #   # segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=org)
      #   # 
      #   #REML eGRM
      #   points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=19, col=blu)
      #   power <- power_results$power_REML_eGRM
      #   std <- sqrt(power*(1-power)/nreps)
      #   segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col=blu)
      # 
      #   # #REML GRM
      #   # points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=19, col=pin)
      #   # power <- power_results$power_REML_GRM
      #   # std <- sqrt(power*(1-power)/nreps)
      #   # segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=pin)
      #   # 
      #   # if(run_acat){
      #   #   points(y=power_results$power_acat, x=x_pos+2*offset, pch=19, col="black")
      #   #   power <- power_results$power_acat
      #   #   std <- sqrt(power*(1-power)/nreps)
      #   #   segments(x0=x_pos+2*offset, y0=power - std, y1=power + std, col="black")
      #   # }
      #   # 
      #   # if(run_arg_needle){
      #   #   power_arg_needle <- read.table(paste("arg_needle_trees/oneRegion/eGRM_GRM/window_based/", ws_causal,"/tested",ws_testing, "/", allowTyped, "/propCausal", propCausal,"/h", h, "/power_results.txt", sep=''), header=TRUE)$power
      #   #   points(y=power_arg_needle, x=x_pos+3*offset, pch=17, col=blu)
      #   #   power <- power_arg_needle
      #   #   std <- sqrt(power*(1-power)/nreps)
      #   #   segments(x0=x_pos+3*offset, y0=power - std, y1=power + std, col=blu)
      #   # }
      #   # 
      #   # 
  
        #true trees
        power_results <- read.table(paste("binary/true_trees/oneRegion/eGRM_GRM/window_based/", ws_causal,"/tested",ws_testing, "/", allowTyped, "/propCausal", propCausal, "/h", h, "/power_results.txt", sep=''), header=TRUE)
      
        
        # GWAS
        points(y=power_results$power_GWAS, x=x_pos-offset, pch=1, col=org)
        power <- power_results$power_GWAS
        std <- sqrt(power*(1-power)/nreps)
        segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=org)
        # 
        # REML eGRM
        points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=1, col=blu)
        power <- power_results$power_REML_eGRM
        std <- sqrt(power*(1-power)/nreps)
        segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col=blu)
        
        # # REML GRM
        # points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=1, col=pin)
        # power <- power_results$power_REML_GRM
        # std <- sqrt(power*(1-power)/nreps)
        # segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=pin)
  
        # #relate trees all variants
        # allowTyped <- "onlyUntyped"   #allowTyped
        # 
        # if(run_acat){
        #   power_results <- read.table(paste("relate_trees_allVariants/oneRegion/eGRM_GRM/window_based/", ws_causal,"/tested",ws_testing, "/", allowTyped, "/propCausal", propCausal,  "/h", h, "/power_results_acat.txt", sep=''), header=TRUE)
        # } else {
        #   power_results <- read.table(paste("relate_trees_allVariants/oneRegion/eGRM_GRM/window_based/", ws_causal,"/tested",ws_testing, "/", allowTyped, "/propCausal", propCausal, "/h", h, "/power_results.txt", sep=''), header=TRUE)
        # }
        # 
        # # REML GWAS
        # points(y=power_results$power_GWAS, x=x_pos-offset, pch=2, col=org)
        # power <- power_results$power_GWAS
        # std <- sqrt(power*(1-power)/nreps)
        # segments(x0=x_pos-offset, y0=power - std, y1=power + std, col=org)
        # 
        # # REML eGRM
        # points(y=power_results$power_REML_eGRM, x=x_pos+0*offset, pch=2, col=blu)
        # power <- power_results$power_REML_eGRM
        # std <- sqrt(power*(1-power)/nreps)
        # segments(x0=x_pos+0*offset, y0=power - std, y1=power + std, col=blu)
        # 
        # # REML GRM
        # points(y=power_results$power_REML_GRM, x=x_pos+offset, pch=2, col=pin)
        # power <- power_results$power_REML_GRM
        # std <- sqrt(power*(1-power)/nreps)
        # segments(x0=x_pos+offset, y0=power - std, y1=power + std, col=pin)
        # 
        # # acat
        # if(run_acat){
        #   points(y=power_results$power_acat, x=x_pos+2*offset, pch=2, col="black")
        #   power <- power_results$power_acat
        #   std <- sqrt(power*(1-power)/nreps)
        #   segments(x0=x_pos+2*offset, y0=power - std, y1=power + std, col="black")
        # }
      }
  
      legend(x="bottomright", legend=c("true trees / all variants","Relate / typed variants", "local eGRM", "local GRM", "GWAS"), pch=c(1, 19, 15, 15, 15), col=c("gray","gray", blu, pin, org), bty='n')
      
    }
  }
}
  dev.off()



