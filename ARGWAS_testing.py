#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""


import tskit
import matplotlib.pyplot as plt
import TPhenotypes as pt
import TGWAS as gwas
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import datetime
import argparse
import logging
import os
import sys
import numpy as np

os.chdir(os.path.dirname(sys.argv[0]))


from numpy.random import RandomState
class randomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)
        

tim = datetime.datetime.now()
seed = tim.hour*10000+tim.minute*100+tim.second
r = randomGenerator(seed)
r.random.get_state()[1][0]

logger = logging.getLogger()


#if you check end of last tree you still get 249,250,621, the length of chr1 in hg19....
# trees.aslist()[trees.num_trees - 1].interval.right
#A site defines a particular location along the genome in which we are interested in observing the allelic state. So I guess that means sites are only defined where there are mutations
# print(trees.tables.sites) #the position of the last site is very close to the end of the interval we kept above

#-----------------------
# Simulate
#-----------------------

simulator = tsim.TSimulatorStdPopsim()
trees = simulator.run_simulation("default", r)

#-----------------------
# Read trees
#-----------------------

directory = "/data/ARGWAS/experiments_N500/"

trees = tskit.load(directory + "test_2.trees")

samp_ids = trees.samples()
N = len(samp_ids)

trees_relate = tskit.load(directory + "relate_test_2_propTyped1/test_2_propTyped1.trees")
variants = tvar.TVariantsFiltered(trees, samp_ids, 0, 1, prop_typed_variants = 1, pos_int = True, random = r, logfile = logger,
                                  filtered_variants_file = directory + "/relate_test_2_propTyped1/test_2_propTyped1_filtered_sample_variants.csv")
print(variants.number_typed)
len(list(trees_relate.variants(samples=samp_ids)))

trees_relate = tskit.load(directory + "relate_test_2_propTyped0.01/test_2_propTyped0.01.trees")
variants = tvar.TVariantsFiltered(trees, samp_ids, 0, 1, prop_typed_variants = 1, pos_int = True, random = r, logfile = logger,
                                  filtered_variants_file = directory + "relate_test_2_propTyped0.01/test_2_propTyped0.01_filtered_sample_variants.csv")
print(variants.number_typed)
len(list(trees_relate.variants(samples=samp_ids)))


trees_relate = tskit.load(directory + "relate_test_2_propTyped0.05/test_2_propTyped0.05.trees")
variants = tvar.TVariantsFiltered(trees, samp_ids, 0, 1, prop_typed_variants = 1, pos_int = True, random = r, logfile = logger,
                                  filtered_variants_file = directory + "relate_test_2_propTyped0.05/test_2_propTyped0.05_filtered_sample_variants.csv")
print(variants.number_typed)
len(list(trees_relate.variants(samples=samp_ids)))

#-----------------------
# simulate pheno based on all variants of a tree
#-----------------------
causal_tree = trees.at(49500000)
len(list(causal_tree.mutations()))
for m in causal_tree.mutations():
    print(m)

for m in causal_tree.sites():
    print(m)

#-----------------------
# create diploids and variants
#-----------------------

inds = tind.Individuals(2, N)
variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, prop_typed_variants = 1, pos_int = True, random = r, logfile = logger)
variants.number_typed

#see variants from tree directly
list(trees.variants(samples=samp_ids))[0]
len(list(trees.variants(samples=samp_ids)))

variants_relate = tvar.TVariantsFiltered(trees_relate, samp_ids, 0.01, 1, prop_typed_variants = 1, pos_int = True, random = r, logfile = logger)
len(list(trees_relate.variants(samples=samp_ids)))
variants_relate.number_typed
#-----------------------
# GRM
#------------------------
M_sum = np.zeros(shape=(10000, 10000))  

for v in range(100):
    # gt = np.array([1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0,1,1,0,1,1,1,1,0,0,0])
    gt = np.random.binomial(1, 0.1, size=10000)
    af = np.sum(gt) / len(gt)   
    first = np.array([gt - af]).T
    second = np.array([gt - af])
    M = np.dot(first, second)
    M = M / (af * (1 - af))
    len(gt)
    M_sum += (M / (100))
np.trace(M_sum)


#--------------------------
# create phenotypes
#-----------------------

# phenotypes with genetic influence
sd_environmental_noise = 1
prop_causal_mutations = 0.002 #this is only for variants found in sampled haplotypes
sd_beta_causal_mutations = 1
pheno_unif = pt.Phenotypes("uniform distr. of causal SNPs", variants, N)
pheno_unif.simulateEnvNoise(sd_environmental_noise, r)
pheno_unif.simulateUniform(variants, prop_causal_mutations=prop_causal_mutations, sd_beta_causal_mutations=sd_beta_causal_mutations, random=r)

# phenotypes with genetic influence, no noise
pheno_unif_noNoise = pt.Phenotypes("uniform distr., no noise", variants, N)
pheno_unif_noNoise.simulateFixed(pheno_unif.causal_variants, pheno_unif.causal_variant_indeces, pheno_unif.causal_betas)
pheno_unif_noNoise.findCausalTrees(trees)

# random phenotypes
sd_environmental_noise = 1
pheno_random = pt.Phenotypes("random", variants, N)
pheno_random.simulateEnvNoise(sd_environmental_noise, r)

# fixed causal variant
sd_environmental_noise = 0
pheno_fixed = pt.Phenotypes("fixed beta -0.67, no noise", variants, N)
pheno_fixed.simulateFixed([pheno_unif.causal_variants[2]], pheno_unif.causal_variant_indeces[2], [-0.67])

# fixed causal variant with high allele freq
sd_environmental_noise = 0
index = np.where(variants.allele_frequencies > 0.4)[0][int(np.floor(len(np.where(variants.allele_frequencies > 0.4)[0]) / 2))]
pheno_fixed_hp = pt.Phenotypes("fixed high freq beta -0.67, no noise", variants, N)
pheno_fixed_hp.simulateFixed([variants.variants[index]], index, [-0.67])

# fixed causal variant with high allele freq with noise
sd_environmental_noise = 1
pheno_fixed_hp_wn = pt.Phenotypes("fixed high freq beta -0.67, with noise", variants, N)
pheno_fixed_hp_wn.simulateEnvNoise(sd_environmental_noise, r)
pheno_fixed_hp_wn.simulateFixed([variants.variants[index]], index, [-0.67])





#-----------------------
# run association tests and plot
#-----------------------

pGWAS_unif = gwas.TpGWAS(phenotypes=pheno_unif)
pGWAS_unif.OLS(variants)
pGWAS_unif_noNoise = gwas.TpGWAS(phenotypes=pheno_unif_noNoise)
pGWAS_unif_noNoise.OLS(variants)

pGWAS_random = gwas.TpGWAS(phenotypes=pheno_random)
pGWAS_random.OLS(variants)

pGWAS_fixed = gwas.TpGWAS(phenotypes=pheno_fixed)
pGWAS_fixed.OLS(variants)

pGWAS_fixed_hp = gwas.TpGWAS(phenotypes=pheno_fixed_hp)
pGWAS_fixed_hp.OLS(variants)
pGWAS_fixed_hp_wn = gwas.TpGWAS(phenotypes=pheno_fixed_hp_wn)
pGWAS_fixed_hp_wn.OLS(variants)


fig, ax = plt.subplots(7,figsize=(30,30))
pGWAS_unif.manhattan_plot(variants.positions, ax[0])
# ax[0].axhline(y=30, color="black", lw=0.5)

pGWAS_unif_noNoise.manhattan_plot(variants.positions, ax[1])

pGWAS_random.manhattan_plot(variants.positions, ax[2])

pGWAS_fixed.manhattan_plot(variants.positions, ax[3])

pGWAS_fixed_hp.manhattan_plot(variants.positions, ax[4])

pGWAS_fixed_hp_wn.manhattan_plot(variants.positions, ax[5])

pGWAS_fixed_hp.manhattan_plot_subset(variants.positions, ax[6], pheno_fixed_hp.causal_variant_indeces-200, pheno_fixed_hp.causal_variant_indeces+200, size=1.5)

fig.tight_layout()
fig.set_size_inches(30, 30)
fig.show()
fig.savefig('sims/sims_13_africans.png', bbox_inches='tight')# 

#-----------------------
# Mantel
#-----------------------
pheno_fixed_hp.findCausalTrees(trees)
tGWAS_fixed_hp = gwas.TtGWAS(trees, pheno_fixed_hp)
tGWAS_fixed_hp.runMantel(trees, pheno_fixed, N)

fig, ax = plt.subplots(2,figsize=(30,30))
pGWAS_fixed_hp.manhattan_plot(variants.positions, ax[0])
tGWAS_fixed_hp.manhattan_plot(range(trees.num_trees), ax[1])

fig.tight_layout()
fig.set_size_inches(30, 30)
fig.show()
fig.savefig('sims/sims_13_randomSeq.png', bbox_inches='tight')# 

#-----------------------
# GCTA HE
#-----------------------

pheno_fixed_hp.findCausalTrees(trees)
tGWAS_fixed_hp = gwas.TtGWAS(trees, pheno_fixed_hp)
tGWAS_fixed_hp.runCGTA_HE(trees, N)

fig, ax = plt.subplots(5,figsize=(30,30))
pGWAS_fixed_hp.manhattan_plot(variants.positions, ax[0])
tGWAS_fixed_hp.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp.p_values_HECP_Jackknife, ax[1], title_supplement = "HECP_Jackknife")
tGWAS_fixed_hp.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp.p_values_HECP_OLS, ax[2], title_supplement = "HECP_OLS")
tGWAS_fixed_hp.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp.p_values_HESD_Jackknife, ax[3], title_supplement = "HESD_Jackknife")
tGWAS_fixed_hp.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp.p_values_HESD_OLS, ax[4], title_supplement = "HESD_OLS")

fig.tight_layout()
fig.set_size_inches(30, 30)
fig.show()
fig.savefig('sims/sims_16_HE.png', bbox_inches='tight')# 

# start = time.time()
# tGWAS_fixed_hp_test = gwas.TtGWAS(trees, pheno_fixed_hp)
# tree = trees.at_index(1141)
# tGWAS_fixed_hp_test.runCGTA_HE_one_tree(trees.at_index(1141), N, start)

# tGWAS_fixed_hp_test.runCGTA_HE(trees, N)

# fig, ax = plt.subplots(5,figsize=(30,30))
# pGWAS_fixed_hp.manhattan_plot(variants.positions, ax[0])
# tGWAS_fixed_hp_test.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_test.p_values_HECP_Jackknife, ax[1], title_supplement = "HECP_Jackknife")
# tGWAS_fixed_hp_test.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_test.p_values_HECP_OLS, ax[2], title_supplement = "HECP_OLS")
# tGWAS_fixed_hp_test.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_test.p_values_HESD_Jackknife, ax[3], title_supplement = "HESD_Jackknife")
# tGWAS_fixed_hp_test.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_test.p_values_HESD_OLS, ax[4], title_supplement = "HESD_OLS")

# fig.tight_layout()
# fig.set_size_inches(30, 30)
# fig.show()
# fig.savefig('sims/sims_15_HE.png', bbox_inches='tight')# 

#-----------------------
# GCTA REML
#-----------------------

pheno_fixed_hp.findCausalTrees(trees)
tGWAS_fixed_hp_reml = gwas.TtGWAS(trees, pheno_fixed_hp)
tGWAS_fixed_hp_reml.runGCTA_REML(trees, N)

fig, ax = plt.subplots(5,figsize=(30,30))
tGWAS_fixed_hp_reml.manhattan_plot(variants.positions, ax[0])
tGWAS_fixed_hp_reml.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_reml.p_values_HECP_Jackknife, ax[1], title_supplement = "HECP_Jackknife")
tGWAS_fixed_hp_reml.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_reml.p_values_HECP_OLS, ax[2], title_supplement = "HECP_OLS")
tGWAS_fixed_hp_reml.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_reml.p_values_HESD_Jackknife, ax[3], title_supplement = "HESD_Jackknife")
tGWAS_fixed_hp_reml.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS_fixed_hp_reml.p_values_HESD_OLS, ax[4], title_supplement = "HESD_OLS")

fig.tight_layout()
fig.set_size_inches(30, 30)
fig.show()
fig.savefig('sims/sims_15_HE.png', bbox_inches='tight')# 




# #-----------------------
# # limix
# #-----------------------

# trees_obj = tt.TTrees(trees)

# def solving_function(array):   
#     covariance = trees_obj.TMRCA(trees_obj.trees[0], N)
#     inv = np.linalg.inv(covariance)
#     tmp = np.dot(inv, array)
#     # print("shape of my dot product",np.shape(tmp))
#     return(tmp)

# y = pheno_random.y
# k = 1
# m = 2
# E = random.randn(N,k)
# N = 500
# # F = sp.concatenate([sp.ones((N,1)), random.randn(N,1)], 1)
# F = sp.zeros(N)


# lmm = LMMCore(y.reshape(N,1), F.reshape(N,1), solving_function)

# S = trees_obj.trees[0].interval.right - trees_obj.trees[0].interval.left #this produces a float, but we want number of sites. there is a sites iterator, but dont know how to use it to find number of sites. why are genomic positions floats?
# G = 0.*(random.rand(N,1)<0.2)
# Inter = random.randn(N, m)

# lmm.process(G, Inter) #this needs step param to produce only one p-value per tree. for this i need the number of sites per tree, or just use 1?
# pv = lmm.getPv()
# beta = lmm.getBetaSNP()
# beta_ste = lmm.getBetaSNPste()
# lrt = lmm.getLRT() #likelihood ratio

# tGWAS_fixed = gwas.TtGWAS(trees, pGWAS_fixed)
# F = np.zeros(N)
# tGWAS_fixed.runLimix(trees, N, pheno_fixed.y.reshape(N,1), F.reshape(N,1), random)


#-----------------------
# plot p-values
#-----------------------

# num_bins = 20
# fig, ax = plt.subplots(4,figsize=(15,15))
# # fig, ax = plt.subplots(1,figsize=(15,15))

# pGWAS_unif.p_value_dist(ax[0], num_bins)

# pGWAS_unif_noNoise.p_value_dist(ax[1], num_bins)

# pGWAS_random.p_value_dist(ax[2], num_bins)
# pGWAS_random.chiSquared(num_bins)

# pGWAS_fixed.p_value_dist(ax[3], num_bins)

# fig.tight_layout()
# fig.set_size_inches(10, 20)
# fig.savefig('sims_pvalues_random_pt_2.png', bbox_inches='tight')# 


#-----------------------
# plot p-values for random
#-----------------------
# num_bins = 20
# fig, ax = plt.subplots(5,figsize=(15,15))

# for i in range(5):
#     sd_environmental_noise = 1
#     pheno_random = pt.Phenotypes("random",trees)
#     pheno_random.simulateEnvNoise(sd_environmental_noise)
#     pGWAS_random = gwas.TpGWAS(ts_object=trees, phenotypes=pheno_random)
#     pGWAS_random.OLS()
    
#     x= pGWAS_random.chiSquared(num_bins)
#     pGWAS_random.p_value_dist(ax[i], num_bins)
#     ax[i].set(title=x)

    
#     fig.tight_layout()
#     fig.set_size_inches(10, 20)
# fig.savefig('sims_pvalues_random_pt.png', bbox_inches='tight')#


# #cumulative p-value distribution
# c_steps = np.arange(0, 1, 0.005)
# c_probs =np.empty(len(c_steps))
# fig, ax = plt.subplots(1,figsize=(15,15))
# for s,c in enumerate(c_steps):
#     c_probs[s] = len(pGWAS_random.p_values[pGWAS_random.p_values <c]) / len(pGWAS_random.p_values)
#     # print(str(c) + ": " + str(len(pGWAS_random.p_values[pGWAS_random.p_values <c]) / len(pGWAS_random.p_values)))
# ax.scatter(c_steps, c_probs)
# ax.set(xlabel='cutoff', ylabel='P(p-value < cutoff)', title="cumulative dist. p-values")
# fig.savefig('sims_pvalues_random_cumulative.png', bbox_inches='tight')#


# fig, ax = plt.subplots(10,figsize=(15,30))
# for i,index in enumerate(range(20000,21000, 100)):
#     pGWAS_random.manhattan_plot_subset(variants.positions, ax[i], index, index+100)
#     ax[i].axhline(y=8, color="black", lw=0.5)
# fig.savefig('manhattan_zooms_random.png', bbox_inches='tight')#


#-----------------------
# test if more causal SNPs bring down p-values
#-----------------------

# # allele_freq = [0.005,0.01, 0.02, 0.05, 0.1, 0.2,0.5]
# num_freqs = 5
# fig, ax = plt.subplots(num_freqs,figsize=(15,15))

# offset = np.floor(len(variants.allele_frequencies)/num_freqs)
# ordered_allele_freq = np.argsort(variants.allele_frequencies)
# for i in range(num_freqs):
#     index = int(offset + i*offset)
#     variant = list(trees.variants(samples=samp_ids))[ordered_allele_freq[index]]
#     freq = variants.allele_frequencies[ordered_allele_freq[index]]
#     pheno = pt.Phenotypes("allele freq=" + str(freq) + ", beta=0.79, with noise" , trees)
#     pheno.simulateEnvNoise(sd_environmental_noise=1)
#     pheno.simulateFixed([variant], [0.79])

#     pGWAS = gwas.TpGWAS(ts_object=trees, phenotypes=pheno)
#     pGWAS.OLS()
    
#     pGWAS.manhattan_plot(variants.variant_positions, ax[i])


# fig.tight_layout()
# fig.set_size_inches(10, 20)
# fig.show()
# fig.savefig('sims_alleleFreq_africans_withNoise.png', bbox_inches='tight')# 
















# # 
# # 
# # y = np.random.normal(scale=environmental_noise_sd, size=N)

# # #add phenotypic effect to mutations that are uniformly distributed
# # betas = [0] * num_variants
# # causal_positions = []
# # variant_positions = []
# # for v, var in enumerate(trees.variants(samples=samp_ids)):  
    
# #     ## Debugging:
# #     ##------------
# #     # print(var.site.mutations[0].node) # the site ID goes from 0 to len(list(trees.variants(samples=samp_ids)))-1
# #     # break

# #     # for v in trees.variants(samples=samp_ids):
# #     #     print(v)   
    
# #     # len(list(trees.variants(samples=samp_ids)))
    
# #     variant_positions.append(var.site.position)

# #     #causal mutation
# #     r = random.uniform(0,1)
# #     if(r < causal_mutations_prop):
        
# #         #define beta
# #         beta = random.normalvariate(0, causal_mutation_beta_sd)
# #         betas[v] = beta
        
# #         #simulate phenotype
# #         y[var.genotypes == 1] += beta
        
# #         #save causal position
# #         causal_positions.append(var.site.position)
        

        
        
# #-------------------
# # associating
# #-------------------
        
# #test for associations
# # diffs = ut.diff(y, N)
# # p_variants = []
# # for v in trees.variants(samples=samp_ids):
# #   p_variants.append(sm.OLS(v.genotypes, y).fit().pvalues[0])

# p_trees = []
# for tree in trees.trees():
#     if tree.total_branch_length == 0: 
#         continue
#     tmrca = ut.TMRCA(tree, N)
#     p_trees.append(ut.mantel(ut.make(tmrca), diffs))
    

# # #--------------------
# # # manhattan plot
# # #--------------------
# # q_variants = -np.log10(p_variants)
# # plt.scatter(variant_positions, q_variants, s=0.5)
# # for pos in causal_positions:
# #     plt.axvline(x=pos, color="red", lw=0.5)
# # plt.show()
    
    
# #---------------------
# # plot ROC curve
# #---------------------
# """
# https://docs.eyesopen.com/toolkits/cookbook/python/plotting/roc.html
# """

# # if(len(p_variants) != len(causal_positions) != len(p_trees) != len(causal_tree_indeces)):
# #     raise ValueError("lengths of p-values and causal positions/trees are not identical. Cannot create ROC curve.")

# TPR = []
# FPR = []
# cutoffs = np.linspace(start=1, stop=0, num=1000)
# for c in cutoffs:
    
#     #variants
#     true_positives = np.where((np.array(betas) != 0) & (np.array(p_variants) < c))[0] #predicted to be positive and the actual value is also positive
#     TPR.append(float(len(true_positives)) / float(len(causal_positions)))
    
#     false_negatives = np.where((np.array(betas) != 0) & (np.array(p_variants) >= c))[0] # predicted to be negative but the actual value is positive
#     FPR.append(float(len(false_negatives)) / float(num_variants - len(causal_positions)))
    

# plt.plot(FPR,TPR)
# plt.ylabel("FPR")
# plt.xlabel("TPR")
# plt.show()

# df = pd.DataFrame(columns = ["FPR", "TPR"])
# df["FPR"] = FPR
# df["TPR"] = TPR

# df.to_csv("ROC.csv",sep=",", header=True, index=False)
# """
# all inds have phenotye = 0 + noise

# - for each mutation:
#     - choose a beta from N(0,sigma) with certain prob given by prop * 1/ num mutations
#     - check which individuals have the mutation
#     - add the beta to all inds with the mutation
#     - save mutation and the beta associated
#     - keep track of which tree the mutation is located in?
    
# - for each SNP:
#     - perform OLS
    
# - for each tree:
#     - perform Mantel
    
# - plot ROC curve?:
#     - need true positives and false negatives etc.: need to keep track of which trees contain causal mutation 

# """

# N = trees.num_samples
# M = trees.num_mutations
# y = np.random.normal(scale=1, size=N)

