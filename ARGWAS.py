#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
# import numpy as np
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
# import statsmodels.api as sm
# import pickle
# import tqdm
# import TTree as tt
# import scipy as sp
# from limix_lmm.lmm_core import LMMCore
# import time
# import sys

os.chdir(os.path.dirname(sys.argv[0]))

#-----------------------------
# initialize arguments
#-----------------------------

parser = argparse.ArgumentParser(description='Running GWAS on variants and trees.')
parser.add_argument('--task', dest="task", required=True,
                    help = 'The task to be executed (simulate or associate)')
parser.add_argument('--out', dest="out", required=True,
                    help = 'Prefix of all output files')
parser.add_argument('--seed', dest="seed", type=int, default = datetime.datetime.now().hour*10000+datetime.datetime.now().minute*100+datetime.datetime.now().second, 
                    help='Set seed of random generator. Default is time stamp.')
# parser.add_argument('--verbose', dest="verbose", 
#                     help="Write output to screen")
parser.add_argument('--tree_file', dest = "tree_file", 
                    help = "File of simulated trees in tskit format")

#simulating trees
parser.add_argument('--sim_tree_simulator', dest = "sim_tree_simulator", default = "stdPopsim", choices=["stdPopsim"],
                    help = "Method used for simulating. stdPopsim is real human chromosome")
parser.add_argument('--pos_int', type=float, default = True,
                    help = "Should the positions of the variants be transformed into integers. Msprime simulates a continuous genome, so if pos_int is true, the simulated positions are rounded and if one position overlaps the previous, it is moved to the next position in the genome.")


#simulating phenotypes
pty = parser.add_argument_group('phenotypes')
# if args.prox and (args.lport is None or args.rport is None):
#     parser.error("--prox requires --lport and --rport.")
parser.add_argument('--name', dest = "name", 
                    help = "Name of phenotype and GWAS object, will be used for headers in plots")
pty.add_argument('--pty_sd_envNoise', type=float, dest = "pty_sd_envNoise", default = 0, 
                    help = "Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
pty.add_argument('--pty_sim_method', dest = "pty_sim_method", choices=['uniform', 'fixed', 'singleTyped', 'singleUntyped'],
                    help = "Phenotype simulations method")
pty.add_argument('--pty_prop_causal_mutations', type=float, dest = "pty_prop_causal_mutations", default = 0, 
                    help = "Proportion of causal mutations to simulate at uniformly distributed positions if pt.sim_method is set to 'uniform'. If set to 0, there will be no causal mutations simulated randomly")
pty.add_argument('--pty_sd_beta_causal_mutations', type=float, dest = "pty_sd_beta_causal_mutations", 
                    help = "Std. dev. for betas of causal mutations if pty_sim_method is set to 'uniform'.")
parser.add_argument('--pty_fixed_betas', nargs='+', type=float, 
                    help = "Fixed betas of causal mutations if pt.sim_method is set to 'fixed_variants'.")
parser.add_argument('--pty_fixed_variant_indeces', nargs='+', type=int,
                    help = "Indeces of variants that should be simulated as causal if pty_sim_method is set to 'fixed'.")
parser.add_argument('--single_variant_af', type=float,
                    help = "Simulate a single, central causal variant that is typed. If there is no such variant in the given range, will search for one with an allele frequency that is close.")

#run associations
assoc = parser.add_argument_group('associations')
assoc.add_argument('--ass_method', dest = "ass_method", choices = ["GWAS", "ARGWAS", "both"], 
                   help = "Either run only GWAS, ARGWAS or both")


#limit data
parser.add_argument('--min_allele_freq', type=float, default=0.01, 
                    help = "Minimum frequency an allele needs to have to be typed")
parser.add_argument('--max_allele_freq', type=float, default=1, 
                    help = "Maximum frequency an allele needs to have to be typed")
parser.add_argument('--prop_typed_variants', type=float, default = 1, 
                    help = "Proportion of variants that are typed (out of the ones that pass the frequency filter).")
# parser.add_argument('--tree_interval', nargs='+',
#                     help = "Indeces of variants that should be simulated as causal if pty_sim_method is set to 'fixed'.")


args = parser.parse_args()


#-----------------------------
# initialize logfile
#-----------------------------

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logging.basicConfig(
    level = logging.INFO,
    format = '%(asctime)s %(levelname)s %(message)s',
    filename = args.out + ".log",
    filemode = 'w' #this makes the log file be not append I think
)

logger.info("writing output files with prefix '" + args.out + "'")

file_handler = logging.FileHandler(args.out + ".log")
logger.addHandler(logging.StreamHandler())
logger.addHandler(file_handler)

logger.info(str(args))


#-----------------------------
# initialize random generator
#-----------------------------
    
from numpy.random import RandomState
class randomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)

r = randomGenerator(args.seed)

logger.info("randomGenerator seed set to " + str(r.random.get_state()[1][0]))


#-----------------------
# Simulate
#-----------------------

if args.task == "simulate":
    if args.sim_tree_simulator == "stdPopsim":
        simulator = tsim.TSimulatorStdPopsim()
        trees = simulator.run_simulation(args.out, r, logger)
        samp_ids = trees.samples()
        N = len(samp_ids)
        inds = tind.Individuals(2, N)
        variants = tvar.TVariants(trees, samp_ids)
        variants.writeVariantInfo(trees, samp_ids, args.out)
        
    else:
        logger.error("use of any simulator besides stdPopSim not tested")
        raise ValueError("use of any simulator besides stdPopSim not tested")


#-----------------------
# Read simulation
#-----------------------

if args.task == "associate":
    
    logger.info("Reading simulations from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()
    N = len(samp_ids)

    #--------------------------------
    # create diploids and variants
    #--------------------------------
    
    inds = tind.Individuals(1, N)
    inds.writeShapeit2(args.out)
    variants = tvar.TVariantsFiltered(trees, samp_ids, args.min_allele_freq, args.max_allele_freq, args.prop_typed_variants, args.pos_int, r)
    # variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, 0.5, r)
    variants.writeVariantInfo(args.out)
    variants.writeShapeit2(args.out, N)
    # variants.fill_diploidGenotypes(samp_ids)
    

  
    #--------------------------------
    # create phenotypes
    #--------------------------------
    
    pheno = pt.Phenotypes(args.name, variants, N, logger)
    pheno.simulateEnvNoise(args.pty_sd_envNoise, r)
    logger.info("Simulating environmental noise with sd " + str(args.pty_sd_envNoise))
    if args.pty_sim_method == 'uniform':
        logger.info("Simulating phenotypes based on uniformly chosen variants with prop_causal_mutations: " + str(args.pty_prop_causal_mutations) + " and sd_beta_causal_mutations: " + str(args.pty_sd_beta_causal_mutations)) 
        pheno.simulateUniform(variants, prop_causal_mutations=args.pty_prop_causal_mutations, sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r)
 
    elif args.pty_sim_method == 'fixed':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'singleUntyped'")

        logger.info("- Simulating phenotypes based on the following indeces: " + str(args.pty_fixed_variant_indeces) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants, args.pty_fixed_variant_indeces, args.pty_fixed_betas, logger)
        
    elif args.pty_sim_method == 'singleTyped':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'singleUntyped'")
            
        var_index = variants.findVariant(typed=True, freq = args.single_variant_af, interval = [49461796, 49602827] , logfile=logger)
        logger.info("- Simulating a phenotypes based on the following typed variant index: " + str(var_index) + " at position " +  str(variants.info['position'][var_index]) + " with allele freq " + str(variants.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants, [var_index], args.pty_fixed_betas, logger)

        
    elif args.pty_sim_method == 'singleUntyped':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'singleUntyped'")
            
        var_index = variants.findVariant(typed=False, freq = args.single_variant_af, interval = [49461796, 49602827], logfile=logger)   
        logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants.info['position'][var_index]) + " with allele freq " + str(variants.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants, [var_index], args.pty_fixed_betas, logger)


    #--------------------------------
    # run association tests and plot
    #--------------------------------
    
    if args.ass_method == "both":
        logger.info("Running both GWAS and ARGWAS for associating")
    else:
        logger.info("Running " + args.ass_method + " for associating")
    
    if args.ass_method == "GWAS" or args.ass_method == "both":
        pGWAS = gwas.TpGWAS(phenotypes=pheno)
        pGWAS.OLS(variants, logger)
        pGWAS.writeToFile(variants, args.out, logger)
        
        fig, ax = plt.subplots(1,figsize=(10,10))
        pGWAS.manhattan_plot(variants.info['position'], ax, logger)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(args.out + '_OLS_GWAS.png', bbox_inches='tight')# 

    if args.ass_method == "ARGWAS" or args.ass_method == "both":
        
        pheno.findCausalTrees(trees)
        
        tGWAS = gwas.TtGWAS(trees, pheno)
        tGWAS.runCGTA_HE(trees, N, args.out, logger)
        tGWAS.writeToFile(trees, args.out, logger)

        fig, ax = plt.subplots(5,figsize=(30,30))
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_Jackknife, subplot=ax[1], logfile=logger, title_supplement = "HECP_Jackknife")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_OLS, subplot=ax[2],  logfile=logger, title_supplement = "HECP_OLS")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_Jackknife, subplot=ax[3], logfile=logger, title_supplement = "HESD_Jackknife")
        tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_OLS, ax[4], logfile=logger, title_supplement = "HESD_OLS")
        
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.show()
        fig.savefig(args.out + '_HE_ARGWAS.png', bbox_inches='tight')# 
 

