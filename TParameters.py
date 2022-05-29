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

    @staticmethod
    def initialize():
        parser = argparse.ArgumentParser(description='Running association tests on variants and trees.')

        # general arguments
        parser.add_argument('--task', required=True,
                            choices=['simulate', 'impute', 'associate', 'downsampleVariantsWriteShapeit', 'ARGStatistics',
                                     'getTreeAtPosition', 'simulateMoreMutations'],
                            help='The task to be executed (simulate or associate)')
        parser.add_argument('--out', required=True, type=str,
                            help='Prefix of all output files')
        parser.add_argument('--seed', type=int,
                            default=datetime.datetime.now().hour * 10000 + datetime.datetime.now().minute * 100
                                    + datetime.datetime.now().second,
                            help='Set seed of random generator. Default is time stamp.')
        # parser.add_argument('--verbose', dest="verbose", 
        #                     help="Write output to screen")
        parser.add_argument('--tree_file', type=str,
                            help="File of trees to be used for association tests in tskit format")
        parser.add_argument('--variants_file', dest="variants_file",
                            help="File of variant info. Useful if they have been downsampled and you want to use the "
                                 "same typed variants")
        parser.add_argument('--tree_file_simulated',
                            help="File of simulated trees to be used for phenotype simulation in tskit format")
        parser.add_argument('--N_sample_pop', type=int, default=0,
                            help="Number of samples in sample population. Must be used on tree simulated with "
                                 "egrm simulate (reference individuals come after sample individuals")
        parser.add_argument('--N_ref_pop', type=int, default=0,
                            help="Number of samples in reference population. Must be used on tree simulated with "
                                 "egrm simulate (reference individuals come after sample individuals")
        parser.add_argument('--ploidy', type=int, choices=[1, 2],
                            help="Ploidy of individuals. Haplotypes will be assigned to individuals in increasing order")
        parser.add_argument('--ploidy_ref', type=int, choices=[1, 2], default=2,
                            help="Ploidy of individuals. Haplotypes will be assigned to individuals in increasing order")

        # limit data
        parser.add_argument('--min_allele_freq', type=float, default=0.01,
                            help="Minimum frequency an allele needs to have to be typed")
        parser.add_argument('--max_allele_freq', type=float, default=1,
                            help="Maximum frequency an allele needs to have to be typed")
        parser.add_argument('--prop_typed_variants', type=float, default=1,
                            help="Proportion of variants that are typed (out of the ones that pass the frequency "
                                 "filter).")
        parser.add_argument('--trees_interval', nargs='+', type=int,  # default=[49e6, 50e6],
                            help="Only test the trees and variants in this interval for association")
        # simulating trees
        sim = parser.add_argument_group('simulating trees')
        sim.add_argument('--sim_tree_simulator', dest="sim_tree_simulator", default="stdPopsim",
                         choices=["stdPopsim", "msprime"],
                         help="Method used for simulating. stdPopsim is real human chromosome")
        sim.add_argument('--N', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim")
        sim.add_argument('--N_African', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim")
        sim.add_argument('--N_European', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim")
        sim.add_argument('--mu', type=float,
                         help="Mutation rate for simulating with msprime (does not work for stdPopsim")
        sim.add_argument('--AH_tree_pos', type=int,
                         help="Add mutations only to the local tree covering this genomic position")
        sim.add_argument('--recomb_rate', type=float,
                         help="Recombination rate for simulating with msprime (does not work for stdPopsim")
        sim.add_argument('--sequence_length', type=float,
                         help="Sequence length for simulating with msprime (does not work for stdPopsim")
        sim.add_argument('--pos_float', action='store_true',
                         help="Should the positions of the variants be transformed into integers. Msprime simulates a "
                              "continuous genome, so if pos_int is true, the simulated positions are rounded and if "
                              "one position overlaps the previous, it is moved to the next position in the genome.")

        # simulating phenotypes
        pty = parser.add_argument_group('simulating phenotypes')
        pty.add_argument('--pty_sd_envNoise', type=float,
                         help="Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
        pty.add_argument('--pty_h_squared', type=float,
                         help="Heritability (h^2) of simulated phenotype. The environmental noise will be simulated in "
                              "order to satisfy the request.")
        pty.add_argument('--pty_sim_method',
                         choices=['null', 'uniform', 'fixed', 'singleTyped', 'singleUntyped', "allelicHetero",
                                  "oneTree"],
                         help="Phenotype simulations method")
        pty.add_argument('--pty_prop_causal_mutations', type=float, default=0,
                         help="Proportion of causal mutations to simulate at uniformly distributed positions if "
                              "pt.sim_method is set to 'uniform'. If set to 0, there will be no causal mutations "
                              "simulated randomly")
        pty.add_argument('--pty_sd_beta_causal_mutations', type=float,
                         help="Std. dev. for betas of causal mutations if pty_sim_method is set to 'uniform'.")
        pty.add_argument('--pty_fixed_betas', nargs='+', type=float,
                         help="Fixed betas of causal mutations if pt.sim_method is set to 'fixed_variants'.")
        pty.add_argument('--pty_fixed_variant_indeces', nargs='+', type=int,
                         help="Indeces of variants that should be simulated as causal if pty_sim_method is set to "
                              "'fixed'.")
        pty.add_argument('--single_variant_af', type=float,
                         help="Simulate a single, central causal variant that is typed. If there is no such variant "
                              "in the given range, will search for one with an allele frequency that is close.")
        pty.add_argument('--single_variant_interval', nargs='+', type=int, default=[49461796.0, 49602827.0],
                         help="Simulate a single, causal variant that is typed within the here given range. If there "
                              "is no such variant in the given range, will search for one with an allele frequency "
                              "that is close.")
        pty.add_argument('--allelic_hetero_file',
                         help="txt file with columns 'freq', 'typed', 'beta'")
        pty.add_argument('--causal_tree_pos', type=int,
                         help="Simulate phenotype with allelic heterogeneity by making all mutations of the local "
                              "tree covering this genomic position causal")

        # run associations
        assoc = parser.add_argument_group('running association tests')
        # TODO: ass method aim method and covariance type should be given in one parameter, maybe a file. Some
        #  combinations should be made impossible
        assoc.add_argument('--ass_method', choices=["GWAS", "AIM", "both"],
                           help="Either run only GWAS, AIM or both")
        assoc.add_argument('--AIM_method', nargs='+',
                           help="Use either Haseman-Elston or REML to test trees for association")
        assoc.add_argument('--covariance_type', type=str, choices=["scaled", "eGRM", "GRM"],
                           help="Use scaled variance-covariance matrix calculated as the covariance scaled by "
                                "N/trace, or use the eGRM calculated by egrm (Fan et al. 2022)")
        assoc.add_argument('--test_only_tree_at', type=float,
                           help="Only test tree that is overlapping the given position for association")
        assoc.add_argument('--skip_first_tree', action='store_true',
                           help='Do not run association test on first tree')
        assoc.add_argument('--do_imputation', action='store_true',
                           help="Imputation takes a long time to run. Use this parameter if imputation was already run "
                                "and you just want to read in files")
        assoc.add_argument('--imputation_ref_panel_tree_file', type=str,
                           help="tskit object of reference panel. If provided, imputation will be run before GWAS.")
        assoc.add_argument('--imputed_gen_file', type=str,
                           help="If imputation has already been run previously, set do_imputation to False and provide "
                                "the .gen file produced by impute2 with this parameter")
        assoc.add_argument('--genetic_map_file', type=str,
                           help="genetic map with columns for Position [bp], Rate [cM/Mb] and Map [cM]")

        args = parser.parse_args()

        return args
