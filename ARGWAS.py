#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
# import numpy as np
import tskit
import numpy as np
import matplotlib.pyplot as plt
import TPhenotypes as pt
import TGWAS as gwas
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import TTree as tt
import pandas as pd
import datetime
import argparse
from python_log_indenter import IndentedLoggerAdapter
import logging
import os
import sys
# import statsmodels.api as sm
# import pickle
# import tqdm
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
parser.add_argument('--task', required=True, choices=['simulate', 'associate', 'downsampleVariants', 'ARGStatistics', 'getTreeAtPosition'],
                    help = 'The task to be executed (simulate or associate)')
parser.add_argument('--out', required=True, type=str,
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
parser.add_argument('--N', type=int, 
                    help =  "Number of haploid individuals to simulate")
parser.add_argument('--pos_int', type=float, default = True,
                    help = "Should the positions of the variants be transformed into integers. Msprime simulates a continuous genome, so if pos_int is true, the simulated positions are rounded and if one position overlaps the previous, it is moved to the next position in the genome.")


#simulating phenotypes
pty = parser.add_argument_group('phenotypes')
parser.add_argument('--name', default = "default",
                    help = "Name of phenotype and GWAS object, will be used for headers in plots")
pty.add_argument('--pty_sd_envNoise', type=float, default = 0, 
                    help = "Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
pty.add_argument('--pty_sim_method', choices=['uniform', 'fixed', 'singleTyped', 'singleUntyped', "allelicHetero", "oneTree"],
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
pty.add_argument('--single_variant_interval', nargs='+', type=int, default=[49461796.0, 49602827.0],
                    help = "Simulate a single, causal variant that is typed within the here given range. If there is no such variant in the given range, will search for one with an allele frequency that is close.")
pty.add_argument('--allelic_hetero_file',
                     help = "txt file with columns 'freq', 'typed', 'beta'")
pty.add_argument('--causal_tree_pos', type=int,
                     help = "genomic position of causal tree")

#run associations
assoc = parser.add_argument_group('associations')
assoc.add_argument('--ass_method', choices = ["GWAS", "AIM", "both"], 
                   help = "Either run only GWAS, AIM or both")
assoc.add_argument('--AIM_method', nargs='+', #choices = ["HE", "REML"],
                   help = "Use either Haseman-Elston or REML to test trees for association")
assoc.add_argument('--test_only_tree_at', type=float, #choices = ["HE", "REML"],
                   help = "Only test tree that is overlapping the given position for association")
assoc.add_argument('--covariance_scaled', type=bool, default=True,
                   help = "scale the variance-covariance matrix by N/trace")

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
file_handler = logging.FileHandler(args.out + ".log")
logger.addHandler(logging.StreamHandler())
logger.addHandler(file_handler)

logging.basicConfig(
    format = '%(asctime)s %(levelname)s %(message)s',
    filename = args.out + ".log",
    filemode = 'w' #this makes the log file be not append I think
)
logger = IndentedLoggerAdapter(logger)
logger.setLevel(logging.INFO)



#-----------------------------
# Header
#-----------------------------

logger.info("---------------------")
logger.info("AIM 0.1")
logger.info("---------------------")


#print arguments to logfile
logger.info("- The following parameters were passed: " + str(args))
logger.info("- Writing output files with prefix '" + str(args.out) + "'")
logger.info("- Adding plots to the following directory '" + str(args.out) + "_plots'")

plots_dir = args.out + "_plots/"
if not os.path.exists(plots_dir):
    os.mkdir(plots_dir)
    

#-----------------------------    
# initialize random generator
#-----------------------------
    
from numpy.random import RandomState
class randomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)

r = randomGenerator(args.seed)

logger.info("- randomGenerator seed is set to " + str(r.random.get_state()[1][0]))


#-----------------------
# Simulate
#-----------------------

if args.task == "simulate":
    logger.info("- TASK: simulate")
    
    if args.N is None:
        raise ValueError("Must provide sample size with argument '--N'")

    if args.sim_tree_simulator == "stdPopsim":
        simulator = tsim.TSimulatorStdPopsim()
        trees = simulator.run_simulation(args.N, args.out, r, logger)
        samp_ids = trees.samples()
        N = len(samp_ids)
        if args.N != N:
            logger.warning("Number of samples in tree does not match number of samples in arguments")
        inds = tind.Individuals(2, N)
        variants = tvar.TVariants(trees, samp_ids)
        variants.writeVariantInfo(trees, samp_ids, args.out)
        
    else:
        logger.error("use of any simulator besides stdPopSim not tested")
        raise ValueError("use of any simulator besides stdPopSim not tested")
        
#-----------------------
# ARG statistics
#-----------------------   
if args.task == "ARGStatistics":

    logger.info("- TASK: ARGStatistics")    
    logger.info("- Reading tree from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    trees_class = tt.TTrees(trees)
    trees_class.writeStats(trees, args.out, logger)

#-----------------------
# Output single tree
#-----------------------

if args.task == "getTreeAtPosition":
    
    logger.info("- TASK: getTreeAtPosition")    
    logger.info("- Reading tree from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    trees_class = tt.TTrees(trees)
    trees_class.extract_single_tree(trees, args.out, logger, position = args.test_only_tree_at)  #49027865

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



#----------------------------------------------------------------
# Read simulation to simulate phenotypes and perform association
#----------------------------------------------------------------

if args.task == "associate":
    
    if args.tree_file_simulated == None or args.tree_file == None:
        raise ValueError("Both the simulated and estimated trees need to be provided with 'tree_file_simulated' and 'tree_file'.")
    
    logger.info("- TASK: Associate")
    logger.info("- Reading simulated tree used for simulating phenotypes from " + args.tree_file_simulated)
    trees_orig = tskit.load(args.tree_file_simulated)
    logger.info("- Reading tree used for Aim's association testing, and for defining variants to be tested by GWAS, from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()
    N = len(samp_ids)

    #--------------------------------
    # create diploids and variants
    #--------------------------------
    
    inds = tind.Individuals(1, N)
 
    # TODO: find way to save variants in their tskit format without needing to read the original tree. I only need original tree in association task for this. It would be nice if the only tree that needs to be read would be estimated tree
    # do not provide variant file here but have it estimated from tree, otherwise variants and tree won't match (tree only contains typed variants). The variant file is only useful for simulating phenotypes to be able to keep track of untyped variants
    variants = tvar.TVariantsFiltered(trees, samp_ids, args.min_allele_freq, args.max_allele_freq, args.prop_typed_variants, args.pos_int, r, logger)
    
    #variants_orig are used to simulate phenotypes. They need to be consistent with original tree and the typed status that might have been defined earlier with a variants file. 
    #The causal mutation should not be affected by a freq filter
    variants_orig = tvar.TVariantsFiltered(trees_orig, samp_ids, 0, 1, 1, args.pos_int, r, logger, args.variants_file)
    
    
    #--------------------------------
    # create phenotypes
    #--------------------------------    

    pheno = pt.Phenotypes(args.name, variants_orig, N, logger)
    pheno.simulateEnvNoise(args.pty_sd_envNoise, r)
    logger.info("- Simulating random noise with sd " + str(args.pty_sd_envNoise))
    
    if args.pty_sim_method is None:
        raise ValueError("Must provide a phenotype simulation method with --pty_sim_method")

    if args.pty_sim_method == 'uniform':
        logger.info("- Simulating phenotypes based on uniformly chosen variants with prop_causal_mutations: " + str(args.pty_prop_causal_mutations) + " and sd_beta_causal_mutations: " + str(args.pty_sd_beta_causal_mutations)) 
        pheno.simulateUniform(variants_orig, prop_causal_mutations=args.pty_prop_causal_mutations, sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r, logfile=logger)
        pheno.write_to_file(variants_orig, args.out, logger)
 
    elif args.pty_sim_method == 'fixed':
        if args.pty_fixed_betas == None:
            raise ValueError("Must provide beta values provided for 'fixed' phenotype using '--pty_fixed_betas'")

        logger.info("- Simulating phenotypes based on the following indeces: " + str(args.pty_fixed_variant_indeces) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants_orig, args.pty_fixed_variant_indeces, args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)

    elif args.pty_sim_method == 'singleTyped':
        if args.pty_fixed_betas == None:
            raise ValueError("Must provide beta values for 'singleTyped' phenotype using '--pty_fixed_betas'")
        if args.single_variant_af == None:
            raise ValueError("Must provide allele freq values for 'singleTyped' phenotype using '--single_variant_af'")
            
        fig, ax = plt.subplots(1,figsize=(30,30))            
        var_index, pos = variants_orig.findVariant(typed=True, freq = args.single_variant_af, interval = args.single_variant_interval, out = args.out, subplot = ax, random = r, logfile = logger)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

        logger.info("- Simulating a phenotypes based on the following typed variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        pheno.simulateFixed(variants_orig, [var_index], args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)
        
    elif args.pty_sim_method == 'singleUntyped':
        if args.pty_fixed_betas == None:
            raise ValueError("Must provide beta values using '--pty_fixed_betas' for phenotype 'singleUntyped'")
        
        fig, ax = plt.subplots(1,figsize=(30,30))            
        var_index, pos = variants_orig.findVariant(typed=False, freq = args.single_variant_af, interval = args.single_variant_interval, out = args.out, subplot = ax, random = r, logfile = logger)   
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')
 
        logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
        #to know which variants are untyped you need variants from simulated tree, not estimated tree
        if args.variants_file is None:
            raise ValueError("Must provide file with untyped variants to simulate phenotype with 'singleUntyped' model")
        pheno.simulateFixed(variants_orig, [var_index], args.pty_fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)
        
    elif args.pty_sim_method == 'oneTree':
        causal_tree = trees.at(args.causal_tree_pos)
        logger.info("- Simulating phenotypes based on all variants of the tree covering postion " + str(args.causal_tree_pos))
        if args.pty_sd_beta_causal_mutations is None:
            raise ValueError("pty_sd_beta_causal_mutations must be set to simulate phenotype with method 'oneTree'")
        pheno.simulateCausalRegion(variants_orig, left_bound = causal_tree.interval.left, right_bound = causal_tree.interval.right, sd_beta_causal_mutations = args.pty_sd_beta_causal_mutations, random = r, logfile = logger)
        pheno.write_to_file(variants_orig, args.out, logger)

    elif args.pty_sim_method == 'allelicHetero':
        if args.allelic_hetero_file == None:
            raise ValueError("No instruction file provided for allelic heterogeneity simulation")
        if args.single_variant_af != None or args.pty_fixed_betas != None:
            raise ValueError("Provided allele frequency or beta value as well as instruction file for allelic heterogeneity. Can accept only one type of instructions.")
        
        ah_info = pd.read_csv(args.allelic_hetero_file, delimiter = "\t")
        variant_indeces = []
        fixed_betas = []
        sum_betas = 0
        
        logger.info("- Searching for loci with requested allele frequencies.")
        logger.add()
        
        #start plot 
        # TODO: this plotting should not be done here
        fig, ax = plt.subplots(ah_info.shape[0],figsize=(30,30))        
        for index, row in ah_info.iterrows():
            #get allele freq
            f = -1
            if row['freq'] > 0.5:
                f = 1 - row['freq']
                logger.warning("- Allele frequencies above 0.5 are not allowed. Transformed " + str(row['freq'])  + " to " + str(f) + ".")
            else:
                f = row['freq']
            
            #get beta
            if not np.isnan(row['beta']) and not np.isnan(row['power']):
                raise ValueError("Cannot fix power and beta value. One value must be set to 'NA'. Beta is " + str(row['beta']) + " and power is " + str(row['power']))
            if np.isnan(row['beta']):
                beta = np.sqrt(row['power'] / (f * (1 - f)))
                #some betas should be negative
                r_num = r.random.uniform(0,1,1)
                if r_num < 0.5:
                    beta = -beta
            else:
                beta = row['beta']


            var_index, pos = variants_orig.findVariant(typed=False, freq = f, interval = [row["interval_start"], row["interval_end"]], out = args.out, subplot = ax[index], random = r, logfile = logger)   
            variant_indeces.append(var_index)
            fixed_betas.append(beta)
            sum_betas += beta

            fig.tight_layout()
            fig.set_size_inches(30, 30)
            fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

        logger.sub()

        logger.info("- Simulating allelic heterogeneous phenotype with total beta " + str(sum_betas))
        
        logger.info("- Simulating phenotypes:")
        logger.add()
        pheno.simulateFixed(variants_orig, variant_indeces, fixed_betas, logger)
        pheno.write_to_file(variants_orig, args.out, logger)
        logger.sub()

    #--------------------------------
    # run association tests and plot
    #--------------------------------
    
    if args.ass_method == "both":
        logger.info("- Running both GWAS and AIM for associating")
    else:
        logger.info("- Running " + args.ass_method + " for associating")
    
    if args.ass_method == "GWAS" or args.ass_method == "both":
        
        logger.info("- GWAS:")
        logger.add()
        pGWAS = gwas.TpGWAS(phenotypes = pheno, num_typed_variants = variants.number_typed)
        pGWAS.OLS(variants, logger)
        pGWAS.writeToFile(variants, args.out, logger)
        
        fig, ax = plt.subplots(1,figsize=(10,10))
        pGWAS.manhattan_plot(variants.info['position'], ax, logger)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(plots_dir + 'OLS_GWAS.png', bbox_inches='tight')# 
        
        logger.sub()

    if args.ass_method == "AIM" or args.ass_method == "both":
        
        logger.info("- AIM:")
        logger.add()
        
        if args.AIM_method is None:
            raise ValueError("ERROR: No method for tree association provided. Use '--AIM_method' to set method.")
            
        logger.info("- Reading tree estimations for tree-based association from " + args.tree_file)
        
        pheno.findCausalTrees(trees)        
            
        for m in args.AIM_method:         
                        
            logger.add()

            if m == "HE":
                logger.info("- Using GCTA Haseman-Elston to test for association between trees and phenotypes")
                
                logger.add()

                tGWAS = gwas.HE_tGWAS(trees, pheno)
                           
                if args.test_only_tree_at is None:    
                    tGWAS.run_association(trees, N, args.out, logger, args.covariance_scaled)
                else:
                    pheno.write_to_file_gcta(args.out, logger)        
                    tree = trees.at(args.test_only_tree_at)
                    tGWAS.run_association_one_tree(tree, N, args.out, logger, args.covariance_scaled)

                tGWAS.write_to_file(trees, args.out, logger)
                       
                # # TODO: move plotting function to tGWAS, should accept p values as argument
                # fig, ax = plt.subplots(5,figsize=(30,30))
                # tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_Jackknife, subplot=ax[1], logfile=logger, title_supplement = "HECP_Jackknife")
                # tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HECP_OLS, subplot=ax[2],  logfile=logger, title_supplement = "HECP_OLS")
                # tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_Jackknife, subplot=ax[3], logfile=logger, title_supplement = "HESD_Jackknife")
                # tGWAS.manhattan_plot_special_pvalues(range(trees.num_trees), tGWAS.p_values_HESD_OLS, ax[4], logfile=logger, title_supplement = "HESD_OLS")
                
                # fig.tight_layout()
                # fig.set_size_inches(30, 30)
                # fig.show()
                # fig.savefig(plots_dir + 'HE_AIM.png', bbox_inches='tight')#    
                
                logger.sub()

                
            if m == "REML":
                logger.info("- Using GCTA REML to test for association between trees and phenotypes")
                
                logger.add()

                tGWAS = gwas.REML_tGWAS(trees, pheno)
                
                if args.test_only_tree_at is None:    
                    tGWAS.run_association(trees, N, args.out, logger, args.covariance_scaled)
                else:
                    pheno.write_to_file_gcta(args.out, logger)        
                    tree = trees.at(args.test_only_tree_at)
                    tGWAS.run_association_one_tree(tree, N, args.out, logger, args.covariance_scaled)
                tGWAS.write_to_file(trees, args.out, logger)         
                
                logger.sub()
            
            logger.sub()
     
        logger.sub()

