#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 17:19:55 2022

@author: linkv
"""
import argparse
import datetime


class TParameters:
    def __init__(self):
        pass

    def initialize(self):
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
        parser.add_argument('--ploidy', type=int, choices=[1,2],
                            help = "Ploidy of individuals. Haplotypes will be assigned to individuals in increasing order")
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
        # TODO: ass method aim method and covariance type should be given in one parameter, maybe a file. Some combinations should be made impossible
        assoc.add_argument('--ass_method', choices = ["GWAS", "AIM", "both"], 
                           help = "Either run only GWAS, AIM or both")
        assoc.add_argument('--AIM_method', nargs='+', #choices = ["HE", "REML"],
                           help = "Use either Haseman-Elston or REML to test trees for association")
        assoc.add_argument('--covariance_type', type=str, choices=["scaled", "eGRM", "GRM"],
                           help = "Use scaled variance-covariance matrix calculated as the covariance scaled by N/trace, or use the eGRM calculated by egrm (Fan et al. 2022)")
        assoc.add_argument('--test_only_tree_at', type=float, #choices = ["HE", "REML"],
                           help = "Only test tree that is overlapping the given position for association")

        
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
        
        return(args)