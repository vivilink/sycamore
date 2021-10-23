#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
# import numpy as np
import tskit
import matplotlib.pyplot as plt
import random
import TPhenotypes as pt
import TGWAS as gwas
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import datetime
import argparse
# import statsmodels.api as sm
# import pickle
# import tqdm
# import TTree as tt
# import scipy as sp
# from limix_lmm.lmm_core import LMMCore
# import time
# import sys


#-----------------------------
# initialize arguments
#-----------------------------

parser = argparse.ArgumentParser(description='Running GWAS on variants and trees.')
parser.add_argument('--task', dest="task", required=True,
                    help = 'The task to be executed (simulate or associate)')
parser.add_argument('--out', dest="out", required=True,
                    help = 'Prefix of all output files')
parser.add_argument('--seed', dest="seed", default = datetime.datetime.now().hour*10000+datetime.datetime.now().minute*100+datetime.datetime.now().second, 
                    help='Set seed of random generator. Default is time stamp.')
parser.add_argument('--tree_file', dest = "tree_file", 
                    help = "File of simulated trees in tskit format")
parser.add_argument('--name', dest = "name", 
                    help = "Name of phenotype and GWAS object, will be used for headers in plots")


#simulating trees
parser.add_argument('--sim.tree_simulator', dest = "sim.tree_simulator")


#simulating phenotypes
pty = parser.add_argument_group('phenotypes')
pty.add_argument('--pt.sd_envNoise', dest = "pty.sd_envNoise", default = 0, 
                    help = "Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
pty.add_argument('--pt.sim_method', dest = "pty.sim_method", choices=['uniform', 'fixed'],
                    help = "Phenotype simulations method")
pty.add_argument('--pt.prop_causal_mutations', dest = "pty.prop_causal_mutations", default = 0, 
                    help = "Proportion of causal mutations to simulate at uniformly distributed positions if pt.sim_method is set to 'uniform'. If set to 0, there will be no causal mutations simulated randomly")
pty.add_argument('--pt.sd_beta_causal_mutations', dest = "pty.sd_beta_causal_mutations", 
                    help = "Std. dev. for betas of causal mutations if pt.sim_method is set to 'uniform'.")
pty.add_argument('--pt.fixed_betas', dest = "pty.fixed_betas", 
                    help = "Fixed betas of causal mutations if pt.sim_method is set to 'fixed_variants'.")
pty.add_argument('--pt.fixed_variant_indeces', dest = "pty.fixed_variant_indeces", 
                    help = "Indeces of variants that should be simulated as causal if pt.sim_method is set to 'fixed'.")


#run associations
assoc = parser.add_argument_group('associations')
assoc.add_argument('--ass.method', dest = "--ass.method", choices = ["GWAS", "ARGWAS", "both"], 
                   help = "Either run only GWAS, ARGWAS or both")

args = parser.parse_args()

#-----------------------------
# initialize random generator
#-----------------------------
    
from numpy.random import RandomState
class randomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)

r = randomGenerator(args.seed)
print("randomGenerator seed set to:",r.random.get_state()[1][0])


#-----------------------
# Simulate
#-----------------------

if args.task == "simulate":
    if args.sim.tree_simulator == "stdPopsim":
        simulator = tsim.TSimulatorStdPopsim()
        trees = simulator.run_simulation(args.out, r)
        samp_ids = trees.samples()
        N = len(samp_ids)
        inds = tind.Individuals(2, N)
        variants = tvar.TVariantsSamples(trees, samp_ids, min_allele_freq = 0.01, max_allele_freq = 1)
        variants.writeAlleleFreq(args.out)
        
    else:
        raise ValueError("use of any simulator besides stdPopSim not tested")


#-----------------------
# Read simulation
#-----------------------

if args.task == "associate":
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()
    N = len(samp_ids)

    #--------------------------------
    # create diploids and variants
    #--------------------------------
    
    inds = tind.Individuals(2, N)
    variants = tvar.TVariantsSamples(trees, samp_ids, 0.01, 1)
    # variants.fill_diploidGenotypes(samp_ids)

  
    #--------------------------------
    # create phenotypes
    #--------------------------------
    
    pheno = pty.Phenotypes(args.name, variants, N)
    pheno.simulateEnvNoise(args.pt.sd_envNoise, r)
    if args.pty.sim_method == 'uniform':
        pheno.simulateUniform(variants, prop_causal_mutations=args.pty.prop_causal_mutations, sd_beta_causal_mutations=args.pty.sd_beta_causal_mutations, random=r)
    elif args.pty.sim_method == 'fixed':
        pheno.simulateFixed(variants, args.fixed_variant_indeces, args.pty.fixed_betas)

    #--------------------------------
    # run association tests and plot
    #--------------------------------
    
    if args.ass.method == "GWAS" or args.ass.method == "both":
        pGWAS = gwas.TpGWAS(phenotypes=pheno)
        pGWAS.OLS(variants)
        pGWAS.writeToFile(trees, args.out)
        
        fig, ax = plt.subplots(1,figsize=(10,10))
        pGWAS.manhattan_plot(variants.positions, ax[0])
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(args.out + '_OLS_GWAS.png', bbox_inches='tight')# 

    if args.ass.method == "ARGWAS" or args.ass.method == "both":
        tGWAS = gwas.TtGWAS(trees, pheno)
        tGWAS.runCGTA_HE(trees, N)
        tGWAS.writeToFile(trees, args.out)

        fig, ax = plt.subplots(5,figsize=(30,30))
        tGWAS.manhattan_plot(variants.positions, ax[0])
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_Jackknife, ax[1], title_supplement = "HECP_Jackknife")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_OLS, ax[2], title_supplement = "HECP_OLS")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_Jackknife, ax[3], title_supplement = "HESD_Jackknife")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_OLS, ax[4], title_supplement = "HESD_OLS")
        
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.show()
        fig.savefig('sims/sims_16_HE.png', bbox_inches='tight')# 
 

