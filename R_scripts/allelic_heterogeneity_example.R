library("amap")
library("ade4")
library("vegan")
library("RAINBOWR")

# data
pheno <- c(0,1,1,2,2,0,0)
h1 <- c(0,1,1,0,0,0,0)
h2 <- c(0,1,1,1,1,1,1)
h3 <- c(0,0,0,0,0,1,1)

# pairwise phenotype distance matrix
pheno_dist <- Dist(pheno, method="manhattan")

# pairwise genotype allele sharing
h1_as <- h1 %*% t(h1)
h2_as <- h2 %*% t(h2)
h3_as <- h3 %*% t(h3)

# GRM

# mantel for individual allele sharing matrices

mantel(pheno_dist, h1_as)
mantel(pheno_dist, h2_as)
mantel(pheno_dist, h3_as)




# correlations

cor(pheno, h1)
cor(pheno, h2)
cor(pheno, h3)