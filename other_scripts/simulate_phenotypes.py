import numpy as np
from pandas_plink import read_plink1_bin
import random

# read genotypes
G = read_plink1_bin("~/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.bed", "~/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.bim", "~/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.fam", verbose=False)

num_variants = len(G.variant)
print("there are",num_variants, "variants")
num_samples = len(G.sample)
print("there are", num_samples, "samples")
heritability = 1

# choose causal variants
num_causal_variants = 100
#variants_causal = random.choices(G.variant, k=num_causal_variants) #variant objects
#variants_causal = random.choices(G.variant.values, k=num_causal_variants) #only variant names
i_variants_causal = random.choices(list(range(num_variants)), k=num_causal_variants)
print("the variants with the following indeces are causal", i_variants_causal)

effect_sizes = np.zeros(num_variants)

# simulate effect sizes per variant
print("simulating effect sizes for causal variants")
for v in i_variants_causal:
    
    variant_name = G.variant.values[v]

    #allele frequency
    geno_var = G.sel(variant=variant_name).values
    num_inds_with_data = np.count_nonzero(~np.isnan(geno_var))
    af = np.nansum(geno_var) / (2 * num_inds_with_data) 

    #effect size
    sd = (heritability / num_causal_variants) / np.sqrt(2 * af * (1 - af))
    beta = np.random.normal(0, sd, 1)[0]
    effect_sizes[v] = beta

# get phenotype for each individual and write to file
print("writing phenotypes to file for each individual")
pheno_file = open("phenotypes.phen", 'w')
for sample_name in G.sample.values:
    genotypes = G.sel(sample=sample_name).values
    phenotype = np.nansum(np.multiply(effect_sizes, genotypes))
    pheno_file.write("0 " + sample_name + " " + str(phenotype) + "\n")

pheno_file.close()



