#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
import numpy as np
import stdpopsim
# import statsmodels.api as sm
# import pickle
# import tqdm
# import tskit
import matplotlib.pyplot as plt
# from itertools import takewhile
# import random
import TPhenotypes as pt
import TGWAS as gwas
import TVariants as tvar
import TTree as tt
import scipy as sp
from limix_lmm.lmm_core import LMMCore


# import msprime

from numpy.random import RandomState
random = RandomState(1)

# simulate 500 haplotypes of chr1 of individuals from Europe, keep only 5 Mb. The species, contig, model, samples, engine and trees are all objects of stdpopsim
# where do we specify ploidy?
species = stdpopsim.get_species("HomSap")
contig = species.get_contig("chr1")
model = species.get_demographic_model("OutOfAfrica_3G09")
samples = model.get_samples(500, 0, 0) # Returns a list of msprime.Sample objects, with the number of samples from each population determined by the positional arguments.
engine = stdpopsim.get_engine("msprime") #returns an engine with a "simulate" method
trees_full = engine.simulate(model, contig, samples) #this runs "msprime.sim_ancestry", default ploidy = 2. Extra arguments passed to simulate are passed to msprime.sim_ancestry

trees = trees_full.keep_intervals([[0,10e6]], simplify=True) 
#if you check end of last tree you still get 249,250,621, the length of chr1 in hg19....
# trees.aslist()[trees.num_trees - 1].interval.right
#A site defines a particular location along the genome in which we are interested in observing the allelic state. So I guess that means sites are only defined where there are mutations
# print(trees.tables.sites) #the position of the last site is very close to the end of the interval we kept above

samp_ids = trees.samples()
N = len(samp_ids)
variants = tvar.TVariantsSamples(trees, samp_ids, 0.01, 1)



# trees_msprime = msprime.sim_ancestry(samples=250, ploidy=2, sequence_length=100000, recombination_rate=0.5)
# trees_msprime_mut = msprime.sim_mutations(trees_msprime, rate=0.01)
# samp_ids = trees_msprime_mut.samples()
# N = len(samp_ids)
# variants = tvar.TVariantsSamples(trees_msprime_mut, samp_ids, 0.01, 1)

# trees = trees_msprime_mut



    
#-----------------------
# create phenotypes
#-----------------------

# phenotypes with genetic influence
sd_environmental_noise = 1
prop_causal_mutations = 0.00002 #this is only for variants found in sampled haplotypes
sd_beta_causal_mutations = 1
pheno_unif = pt.Phenotypes("uniform distr. of causal SNPs", variants, N)
pheno_unif.simulateEnvNoise(sd_environmental_noise)
pheno_unif.simulateUniform(variants, prop_causal_mutations=prop_causal_mutations, sd_beta_causal_mutations=sd_beta_causal_mutations)
# pheno_unif.findCausalTrees(trees)

# phenotypes with genetic influence, no noise
# sd_environmental_noise = 0.0
# prop_causal_mutations = 0.001 #this is only for variants found in sampled haplotypes
# sd_beta_causal_mutations = 1
pheno_unif_noNoise = pt.Phenotypes("uniform distr., no noise", variants, N)
pheno_unif_noNoise.simulateFixed(pheno_unif.causal_variants, pheno_unif.causal_variant_indeces, pheno_unif.causal_betas)
pheno_unif_noNoise.findCausalTrees(trees)

# random phenotypes
sd_environmental_noise = 1
pheno_random = pt.Phenotypes("random", variants, N)
pheno_random.simulateEnvNoise(sd_environmental_noise)

# fixed causal variant
sd_environmental_noise = 0
# index = np.where(allele_frequencies > 0.4)[0][1000]
pheno_fixed = pt.Phenotypes("fixed beta -0.67, no noise", variants, N)
# pheno_fixed.simulateFixed([list(trees.variants(samples=samp_ids))[index]], [0.01])
pheno_fixed.simulateFixed([pheno_unif.causal_variants[0]], pheno_unif.causal_variant_indeces[0], [-0.67])

# fixed causal variant with high allele freq
sd_environmental_noise = 0
index = np.where(variants.allele_frequencies > 0.4)[0][int(np.floor(len(np.where(variants.allele_frequencies > 0.4)[0]) / 2))]
pheno_fixed_hp = pt.Phenotypes("fixed high freq beta -0.67, no noise", variants, N)
pheno_fixed_hp.simulateFixed([variants.variants[index]], index, [-0.67])

# fixed causal variant with high allele freq with noise
sd_environmental_noise = 1
pheno_fixed_hp_wn = pt.Phenotypes("fixed high freq beta -0.67, with noise", variants, N)
pheno_fixed_hp_wn.simulateEnvNoise(sd_environmental_noise)
pheno_fixed_hp_wn.simulateFixed([variants.variants[index]], index, [-0.67])


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






# fig, ax = plt.subplots(2,figsize=(30,30))
# pGWAS_random.manhattan_plot(variants.positions, ax[0])
# tGWAS_random.manhattan_plot(range(trees.num_trees), ax[1])

# fig.tight_layout()
# fig.set_size_inches(30, 30)
# fig.show()
# fig.savefig('sims/sims_13_randomSeq.png', bbox_inches='tight')# 

#-----------------------
# run association tests and plot
#-----------------------
F = np.zeros(N)

pGWAS_unif = gwas.TpGWAS(phenotypes=pheno_unif)
pGWAS_unif.OLS(variants)
pGWAS_unif_noNoise = gwas.TpGWAS(phenotypes=pheno_unif_noNoise)
pGWAS_unif_noNoise.OLS(variants)

pGWAS_random = gwas.TpGWAS(phenotypes=pheno_random)
pGWAS_random.OLS(variants)
tGWAS_random = gwas.TtGWAS(trees, pheno_random)
tGWAS_random.runLimix(trees, N, pheno_random.y.reshape(N,1), F.reshape(N,1), random)

pGWAS_fixed = gwas.TpGWAS(phenotypes=pheno_fixed)
pGWAS_fixed.OLS(variants)
tGWAS_fixed = gwas.TtGWAS(trees, pGWAS_fixed)
tGWAS_fixed.runLimix(trees, N, pheno_fixed.y.reshape(N,1), F.reshape(N,1), random)

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
fig.savefig('sims/sims_12_africans.png', bbox_inches='tight')# 



#-----------------------
# plot p-values
#-----------------------

num_bins = 20
fig, ax = plt.subplots(4,figsize=(15,15))
# fig, ax = plt.subplots(1,figsize=(15,15))

pGWAS_unif.p_value_dist(ax[0], num_bins)

pGWAS_unif_noNoise.p_value_dist(ax[1], num_bins)

pGWAS_random.p_value_dist(ax[2], num_bins)
pGWAS_random.chiSquared(num_bins)

pGWAS_fixed.p_value_dist(ax[3], num_bins)

fig.tight_layout()
fig.set_size_inches(10, 20)
fig.savefig('sims_pvalues_random_pt_2.png', bbox_inches='tight')# 


#-----------------------
# plot p-values for random
#-----------------------
num_bins = 20
fig, ax = plt.subplots(5,figsize=(15,15))

for i in range(5):
    sd_environmental_noise = 1
    pheno_random = pt.Phenotypes("random",trees)
    pheno_random.simulateEnvNoise(sd_environmental_noise)
    pGWAS_random = gwas.TpGWAS(ts_object=trees, phenotypes=pheno_random)
    pGWAS_random.OLS()
    
    x= pGWAS_random.chiSquared(num_bins)
    pGWAS_random.p_value_dist(ax[i], num_bins)
    ax[i].set(title=x)

    
    fig.tight_layout()
    fig.set_size_inches(10, 20)
fig.savefig('sims_pvalues_random_pt.png', bbox_inches='tight')#


#cumulative p-value distribution
c_steps = np.arange(0, 1, 0.005)
c_probs =np.empty(len(c_steps))
fig, ax = plt.subplots(1,figsize=(15,15))
for s,c in enumerate(c_steps):
    c_probs[s] = len(pGWAS_random.p_values[pGWAS_random.p_values <c]) / len(pGWAS_random.p_values)
    # print(str(c) + ": " + str(len(pGWAS_random.p_values[pGWAS_random.p_values <c]) / len(pGWAS_random.p_values)))
ax.scatter(c_steps, c_probs)
ax.set(xlabel='cutoff', ylabel='P(p-value < cutoff)', title="cumulative dist. p-values")
fig.savefig('sims_pvalues_random_cumulative.png', bbox_inches='tight')#


fig, ax = plt.subplots(10,figsize=(15,30))
for i,index in enumerate(range(20000,21000, 100)):
    pGWAS_random.manhattan_plot_subset(variants.positions, ax[i], index, index+100)
    ax[i].axhline(y=8, color="black", lw=0.5)
fig.savefig('manhattan_zooms_random.png', bbox_inches='tight')#


#-----------------------
# test if more causal SNPs bring down p-values
#-----------------------

# allele_freq = [0.005,0.01, 0.02, 0.05, 0.1, 0.2,0.5]
num_freqs = 5
fig, ax = plt.subplots(num_freqs,figsize=(15,15))

offset = np.floor(len(variants.allele_frequencies)/num_freqs)
ordered_allele_freq = np.argsort(variants.allele_frequencies)
for i in range(num_freqs):
    index = int(offset + i*offset)
    variant = list(trees.variants(samples=samp_ids))[ordered_allele_freq[index]]
    freq = variants.allele_frequencies[ordered_allele_freq[index]]
    pheno = pt.Phenotypes("allele freq=" + str(freq) + ", beta=0.79, with noise" , trees)
    pheno.simulateEnvNoise(sd_environmental_noise=1)
    pheno.simulateFixed([variant], [0.79])

    pGWAS = gwas.TpGWAS(ts_object=trees, phenotypes=pheno)
    pGWAS.OLS()
    
    pGWAS.manhattan_plot(variants.variant_positions, ax[i])


fig.tight_layout()
fig.set_size_inches(10, 20)
fig.show()
fig.savefig('sims_alleleFreq_africans_withNoise.png', bbox_inches='tight')# 
















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
        
# # #find causal tree
# # causal_tree_indeces = []
# # p = 0 
# # t = 0
# # tree = trees.first()
# # while p < len(causal_positions):    
# #     ## Debugging:
# #     ##------------
# #     # print("p: " + str(p))
# #     # print("tree index " + str(t))
# #     # print("causal_positions[p] + " + str(causal_positions[p]))
# #     # print("tree.interval.left " + str(tree.interval.left))
# #     # print("tree.interval.right " + str(tree.interval.right)) 
# #     # print("trees.at(var.site.position).get_index() " + str(trees.at(var.site.position).get_index()))
        
# #     if tree.interval.left <= causal_positions[p] <= tree.interval.right:        
# #         #save causal tree
# #         causal_tree_indeces.append(tree.get_index())
# #         p += 1
        
# #     elif causal_positions[p] < tree.interval.left:
# #         p += 1        
    
# #     elif causal_positions[p] > tree.interval.right:
# #         tree.next()
# #         t += 1
        
        
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

