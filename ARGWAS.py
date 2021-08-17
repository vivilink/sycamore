#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
import numpy as np
import stdpopsim
import utils as ut
import statsmodels.api as sm
import pickle
import tqdm
import tskit
import matplotlib.pyplot as plt
from itertools import takewhile
import pandas as pd
import random
import TPhenotypes as pt







# simulate 500 haplotypes of chr1 of individuals from Europe, keep only 5 Mb. The species, contig, model, samples, engine and trees are all objects of stdpopsim
# where do we specify ploidy?
species = stdpopsim.get_species("HomSap")
contig = species.get_contig("chr1")
model = species.get_demographic_model("OutOfAfrica_3G09")
samples = model.get_samples(0, 500, 0) # Returns a list of msprime.Sample objects, with the number of samples from each population determined by the positional arguments.
engine = stdpopsim.get_engine("msprime") #returns an engine with a "simulate" method
trees = engine.simulate(model, contig, samples) #this runs "msprime.sim_ancestry", default ploidy = 2. Extra arguments passed to simulate are passed to msprime.sim_ancestry

trees = trees.keep_intervals([[0,5e6]], simplify=True)

#-----------------------
# phenotype parameters
#-----------------------

sd_environmental_noise = 0.0
prop_causal_mutations = 0.001 #this is only for variants found in sampled haplotypes
sd_beta_causal_mutations = 1

#simulate phenotypes for all inds or haplotypes? trees.num_individuals = 0, but num_samples = 500
N = trees.num_samples
samp_ids = trees.samples()
num_variants = len(list(trees.variants(samples=samp_ids)))
pheno = pt.Phenotypes(N=N, num_variants = num_variants)
pheno.simulateUniform(trees.variants(samples=samp_ids), prop_causal_mutations, sd_environmental_noise, sd_beta_causal_mutations)

# 
# 
# 
# y = np.random.normal(scale=environmental_noise_sd, size=N)

# #add phenotypic effect to mutations that are uniformly distributed
# betas = [0] * num_variants
# causal_positions = []
# variant_positions = []
# for v, var in enumerate(trees.variants(samples=samp_ids)):  
    
#     ## Debugging:
#     ##------------
#     # print(var.site.mutations[0].node) # the site ID goes from 0 to len(list(trees.variants(samples=samp_ids)))-1
#     # break

#     # for v in trees.variants(samples=samp_ids):
#     #     print(v)   
    
#     # len(list(trees.variants(samples=samp_ids)))
    
#     variant_positions.append(var.site.position)

#     #causal mutation
#     r = random.uniform(0,1)
#     if(r < causal_mutations_prop):
        
#         #define beta
#         beta = random.normalvariate(0, causal_mutation_beta_sd)
#         betas[v] = beta
        
#         #simulate phenotype
#         y[var.genotypes == 1] += beta
        
#         #save causal position
#         causal_positions.append(var.site.position)
        
#find causal tree
causal_tree_indeces = []
p = 0 
t = 0
tree = trees.first()
while p < len(causal_positions):    
    ## Debugging:
    ##------------
    # print("p: " + str(p))
    # print("tree index " + str(t))
    # print("causal_positions[p] + " + str(causal_positions[p]))
    # print("tree.interval.left " + str(tree.interval.left))
    # print("tree.interval.right " + str(tree.interval.right)) 
    # print("trees.at(var.site.position).get_index() " + str(trees.at(var.site.position).get_index()))
        
    if tree.interval.left <= causal_positions[p] <= tree.interval.right:        
        #save causal tree
        causal_tree_indeces.append(tree.get_index())
        p += 1
        
    elif causal_positions[p] < tree.interval.left:
        p += 1        
    
    elif causal_positions[p] > tree.interval.right:
        tree.next()
        t += 1
        
        
#-------------------
# associating
#-------------------
        
#test for associations
diffs = ut.diff(y, N)
p_variants = []
for v in trees.variants(samples=samp_ids):
  p_variants.append(sm.OLS(v.genotypes, y).fit().pvalues[0])

p_trees = []
for tree in trees.trees():
    if tree.total_branch_length == 0: 
        continue
    tmrca = ut.TMRCA(tree, N)
    p_trees.append(ut.mantel(ut.make(tmrca), diffs))
    

#--------------------
# manhattan plot
#--------------------
q_variants = -np.log10(p_variants)
plt.scatter(variant_positions, q_variants, s=0.5)
for pos in causal_positions:
    plt.axvline(x=pos, color="red", lw=0.5)
plt.show()
    
    
#---------------------
# plot ROC curve
#---------------------
"""
https://docs.eyesopen.com/toolkits/cookbook/python/plotting/roc.html
"""

# if(len(p_variants) != len(causal_positions) != len(p_trees) != len(causal_tree_indeces)):
#     raise ValueError("lengths of p-values and causal positions/trees are not identical. Cannot create ROC curve.")

TPR = []
FPR = []
cutoffs = np.linspace(start=1, stop=0, num=1000)
for c in cutoffs:
    
    #variants
    true_positives = np.where((np.array(betas) != 0) & (np.array(p_variants) < c))[0] #predicted to be positive and the actual value is also positive
    TPR.append(float(len(true_positives)) / float(len(causal_positions)))
    
    false_negatives = np.where((np.array(betas) != 0) & (np.array(p_variants) >= c))[0] # predicted to be negative but the actual value is positive
    FPR.append(float(len(false_negatives)) / float(num_variants - len(causal_positions)))
    

plt.plot(FPR,TPR)
plt.ylabel("FPR")
plt.xlabel("TPR")
plt.show()

df = pd.DataFrame(columns = ["FPR", "TPR"])
df["FPR"] = FPR
df["TPR"] = TPR

df.to_csv("ROC.csv",sep=",", header=True, index=False)
"""
all inds have phenotye = 0 + noise

- for each mutation:
    - choose a beta from N(0,sigma) with certain prob given by prop * 1/ num mutations
    - check which individuals have the mutation
    - add the beta to all inds with the mutation
    - save mutation and the beta associated
    - keep track of which tree the mutation is located in?
    
- for each SNP:
    - perform OLS
    
- for each tree:
    - perform Mantel
    
- plot ROC curve?:
    - need true positives and false negatives etc.: need to keep track of which trees contain causal mutation 

"""

N = trees.num_samples
M = trees.num_mutations
y = np.random.normal(scale=1, size=N)

