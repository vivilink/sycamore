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
import pandas as pd
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

parser = argparse.ArgumentParser(description='Running association tests on variants and trees.')

# general arguments
parser.add_argument('--task', required=True, choices=['simulate', 'associate', 'downsampleVariants'],
                    help = 'The task to be executed (simulate or associate)')
parser.add_argument('--out', required=True,
                    help = 'Prefix of all output files')
parser.add_argument('--seed', type=int, default = datetime.datetime.now().hour*10000+datetime.datetime.now().minute*100+datetime.datetime.now().second, 
                    help='Set seed of random generator. Default is time stamp.')
# parser.add_argument('--verbose', dest="verbose", 
#                     help="Write output to screen")
parser.add_argument('--tree_file',  
                    help = "File of trees to be used for association tests in tskit format")
parser.add_argument('--variants_file', dest = "variants_file",
                    help = "File of variant info. Useful if they have been downsampled and you want to use the same typed variants")
parser.add_argument('--tree_file_simulated',
                    help = "File of simulated trees to be used for phenotype simulation in tskit format")

#simulating trees
parser.add_argument('--sim_tree_simulator', dest = "sim_tree_simulator", default = "stdPopsim", choices=["stdPopsim"],
                    help = "Method used for simulating. stdPopsim is real human chromosome")
parser.add_argument('--pos_int', type=float, default = True,
                    help = "Should the positions of the variants be transformed into integers. Msprime simulates a continuous genome, so if pos_int is true, the simulated positions are rounded and if one position overlaps the previous, it is moved to the next position in the genome.")


#simulating phenotypes
pty = parser.add_argument_group('phenotypes')
parser.add_argument('--name', 
                    help = "Name of phenotype and GWAS object, will be used for headers in plots")
pty.add_argument('--pty_sd_envNoise', type=float, default = 0, 
                    help = "Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
pty.add_argument('--pty_sim_method', choices=['uniform', 'fixed', 'singleTyped', 'singleUntyped', "allelicHetero"],
                    help = "Phenotype simulations method")
pty.add_argument('--pty_prop_causal_mutations', type=float, default = 0, 
                    help = "Proportion of causal mutations to simulate at uniformly distributed positions if pt.sim_method is set to 'uniform'. If set to 0, there will be no causal mutations simulated randomly")
pty.add_argument('--pty_sd_beta_causal_mutations', type=float, 
                    help = "Std. dev. for betas of causal mutations if pty_sim_method is set to 'uniform'.")
pty.add_argument('--pty_fixed_betas', nargs='+', type=float, 
                    help = "Fixed betas of causal mutations if pt.sim_method is set to 'fixed_variants'.")
pty.add_argument('--pty_fixed_variant_indeces', nargs='+', type=int,
                    help = "Indeces of variants that should be simulated as causal if pty_sim_method is set to 'fixed'.")
pty.add_argument('--single_variant_af', type=float,
                    help = "Simulate a single, central causal variant that is typed. If there is no such variant in the given range, will search for one with an allele frequency that is close.")
pty.add_argument('--allelic_hetero_file',
                     help = "txt file with columns 'freq', 'typed', 'beta'")

#run associations
assoc = parser.add_argument_group('associations')
assoc.add_argument('--ass_method', choices = ["GWAS", "AIM", "both"], 
                   help = "Either run only GWAS, AIM or both")



#limit data
parser.add_argument('--min_allele_freq', type=float, default = 0, 
                    help = "Minimum frequency an allele needs to have to be typed")
parser.add_argument('--max_allele_freq', type=float, default = 1, 
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
# Downsample variants
#-----------------------
if args.task == "downsampleVariants":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariant'")
    logger.info("- TASK: Downsampling variants")
    logger.info("- Reading tree simulations from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()
    N = len(samp_ids)

    #--------------------------------
    # create diploids and variants
    #--------------------------------
    
    inds = tind.Individuals(1, N)
    inds.writeShapeit2(args.out, logger)
    variants = tvar.TVariantsFiltered(trees, samp_ids, args.min_allele_freq, args.max_allele_freq, args.prop_typed_variants, args.pos_int, r, logger)
    # variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, 0.5, r)
    variants.writeVariantInfo(args.out, logger)
    variants.writeShapeit2(args.out, N, logger)



#-----------------------
# Read simulation
#-----------------------

if args.task == "associate":
    
    logger.info("- TASK: Associate")
    
    logger.info("- Reading tree from " + args.tree_file)
    trees_orig = tskit.load(args.tree_file_simulated)
    logger.info("- Reading simulated tree from " + args.tree_file_simulated)
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()
    N = len(samp_ids)

    #--------------------------------
    # create diploids and variants
    #--------------------------------
    
    inds = tind.Individuals(1, N)
    # TODO: find way to save variants in their tskit format without needing to read the original tree. I only need original tree in association task for this. It would be nice if the only tree that needs to be read would be estimated tree
    variants = tvar.TVariantsFiltered(trees, samp_ids, args.min_allele_freq, args.max_allele_freq, args.prop_typed_variants, args.pos_int, r, logger)

    #--------------------------------
    # create phenotypes
    #--------------------------------
    variants_orig = tvar.TVariantsFiltered(trees_orig, samp_ids, 0, 1, 1, args.pos_int, r, logger, args.variants_file)
    print("len(variants_orig.variants) after reading variants_orig",len(variants_orig.variants))
    pheno = pt.Phenotypes(args.name, variants_orig, N, logger)
    pheno.simulateEnvNoise(args.pty_sd_envNoise, r)
    logger.info("- Simulating environmental noise with sd " + str(args.pty_sd_envNoise))
    if args.pty_sim_method == 'uniform':
        logger.info("- Simulating phenotypes based on uniformly chosen variants with prop_causal_mutations: " + str(args.pty_prop_causal_mutations) + " and sd_beta_causal_mutations: " + str(args.pty_sd_beta_causal_mutations)) 
        pheno.simulateUniform(variants_orig, prop_causal_mutations=args.pty_prop_causal_mutations, sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r)
        pheno.write_to_file(variants_orig, args.out, logger)
 
    elif args.pty_sim_method == 'fixed':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'fixed'")

        logger.info("- Simulating phenotypes based on the following indeces: " + str(args.pty_fixed_variant_indeces) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants_orig, args.pty_fixed_variant_indeces, args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)

    elif args.pty_sim_method == 'singleTyped':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'singleTyped'")
            
        var_index = variants_orig.findVariant(typed=True, freq = args.single_variant_af, interval = [49461796, 49602827] , logfile=logger)
        logger.info("- Simulating a phenotypes based on the following typed variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants_orig, [var_index], args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)

        
    elif args.pty_sim_method == 'singleUntyped':
        if args.pty_fixed_betas == None:
            raise ValueError("No beta values provided for phenotype 'singleUntyped'")
            
        var_index = variants_orig.findVariant(typed=False, freq = args.single_variant_af, interval = [49461796, 49602827], logfile=logger)   
        logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        #to know which variants are untyped you need variants from simulated tree, not estimated tree
        if args.variants_file is None:
            raise ValueError("Must provide file with untyped variants to simulate phenotype with 'singleUntyped' model")
        pheno.simulateFixed(variants_orig, [var_index], args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)
        
    elif args.pty_sim_method == 'allelicHetero':
        if args.allelic_hetero_file == None:
            raise ValueError("No instruction file provided for allelic heterogeneity simulation")
        if args.single_variant_af != None or args.pty_fixed_betas != None:
            raise ValueError("Provided allele frequency or beta value as well as instruction file for allelic heterogeneity. Can accept only one type of instructions.")
        
        ah_info = pd.read_csv(args.allelic_hetero_file, delimiter = "\t")
        variant_indeces = []
        fixed_betas = []
        for index, row in ah_info.iterrows():
            var_index = variants_orig.findVariant(typed=False, freq = row['freq'], interval = [49461796, 49602827], logfile=logger)   
            variant_indeces.append(var_index)
            fixed_betas.append(row['beta'])
        logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants_orig, variant_indeces, fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)
            
        # var_index = variants_orig.findVariant(typed=False, freq = args.single_variant_af, interval = [49461796, 49602827], logfile=logger)   
        # print("number of orig variatns", variants_orig.number)
        # logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        # #to know which variants are untyped you need variants from simulated tree, not estimated tree
        # if args.variants_file is None:
        #     raise ValueError("Must provide file with untyped variants to simulate phenotype with 'singleUntyped' model")
        # pheno.simulateFixed(variants_orig, [var_index], args.pty_fixed_betas, logger)
        # pheno.write_to_file(variants_orig, args.out, logger)

    #--------------------------------
    # run association tests and plot
    #--------------------------------
    
    if args.ass_method == "both":
        logger.info("- Running both GWAS and AIM for associating")
    else:
        logger.info("- Running " + args.ass_method + " for associating")
    
    if args.ass_method == "GWAS" or args.ass_method == "both":
        pGWAS = gwas.TpGWAS(phenotypes = pheno, num_typed_variants = variants.number_typed)
        print("variants.number_typed in ARGWAS 2", variants.number_typed)
        pGWAS.OLS(variants, logger)
        pGWAS.writeToFile(variants, args.out, logger)
        
        fig, ax = plt.subplots(1,figsize=(10,10))
        pGWAS.manhattan_plot(variants.info['position'], ax, logger)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(args.out + '_OLS_GWAS.png', bbox_inches='tight')# 

    if args.ass_method == "AIM" or args.ass_method == "both":
        
        logger.info("- Reading tree estimations for tree-based association from " + args.tree_file)
        
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
        fig.savefig(args.out + '_HE_AIM.png', bbox_inches='tight')# 
 

