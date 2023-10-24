import numpy as np
from pandas_plink import read_plink1_bin
from numpy.random import RandomState
import sys
import time

replicate = int(sys.argv[1])
random = RandomState(seed=int(replicate))


# read genotypes
G = read_plink1_bin("/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.bed", "/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.bim", "/home1/linkv/ARGWAS/hawaiian/plink_files_analysis_chromosomes/chr5/All_MEGA_GWAS.and.CIDR_GDA_GWAS.nh.intersect1M.drop5k.chr1-22.asVCF.rmvoutoforder.geno05.hwe1e-6.ld50-5-.8.prune.in.maf01.fam", verbose=False)

# remove variants from chromosome 5
G = G.where((G.chrom == '1') | (G.chrom == '2') |  (G.chrom == '3'), drop=True)

num_variants = len(G.variant)
print("there are",num_variants, "variants")
num_samples = len(G.sample)
print("there are", num_samples, "samples")
heritability = 0.45

# choose causal variants
num_causal_variants = int(sys.argv[2])
print("simulating", num_causal_variants, "variants to be causal. This corresponds to", num_causal_variants / num_variants, "of all variants")
#variants_causal = random.choices(G.variant, k=num_causal_variants) #variant objects
#variants_causal = random.choices(G.variant.values, k=num_causal_variants) #only variant names
i_variants_causal = random.choice(a=list(range(num_variants)), size=num_causal_variants)
print("the variants with the following indeces are causal", i_variants_causal)

effect_sizes = np.zeros(num_variants)

# simulate effect sizes per variant
print("simulating effect sizes for causal variants")
start = time.time()
for s in range(len(i_variants_causal)):

    v = i_variants_causal[s]
    
    variant_name = G.variant.values[v]

    #allele frequency
    geno_var = G.sel(variant=variant_name).values
    num_inds_with_data = np.count_nonzero(~np.isnan(geno_var))
    af = np.nansum(geno_var) / (2 * num_inds_with_data) 

    #effect size
    sd = (heritability / num_causal_variants) / np.sqrt(2 * af * (1 - af))
    beta = random.normal(loc=0, scale=sd, size=1)[0]
    effect_sizes[v] = beta

    if s % 1000 == 0:
        end = time.time()
        print("simulated effect sizes for", s, "variants in", end-start, "s")

# get phenotype for each individual and write to file
print("writing phenotypes to file for each individual")
pheno_file = open("phenotypes_" + str(replicate) + ".phen", 'w')
start = time.time()
#for i in range(10):
phenotypes = np.zeros(num_samples)
for i in range(len(G.sample.values)):
    sample_name = G.sample.values[i]
    genotypes = G.sel(sample=sample_name).values
    phenotypes[i] += np.nansum(np.multiply(effect_sizes[i_variants_causal], genotypes[i_variants_causal]))  

    if i % 100 == 0:
        end = time.time()
        print("simulated phenotye for", i, "individuals in", end - start, "s")

var_genetic = np.var(phenotypes)
print("var_genetic", var_genetic)
var_env = var_genetic * (1 - heritability) / heritability
print("var_env", var_env)
for i in range(len(phenotypes)):
    p += random.normal(loc=0, scale=np.sqrt(var_env), size=1)[0]
    pheno_file.write("0 " + G.sample.values[i] + " " + str(phenotypes[i]) + "\n")


pheno_file.close()



