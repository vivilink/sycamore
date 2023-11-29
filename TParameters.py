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
        parser = argparse.ArgumentParser(description='')

        # general arguments
        parser.add_argument('--task', required=True,
                            choices=['simulate', 'impute', 'associate', 'downsampleVariantsWriteShapeit',
                                     'ARGStatistics', 'getTreeAtPosition', 'simulateMoreMutations',
                                     'covarianceCorrelations', 'simulatePhenotypes', 'makeTreeChunks'],
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
        parser.add_argument('--min_allele_freq', type=float, default=0.0,
                            help="Minimum frequency an allele needs to have to be typed")
        parser.add_argument('--max_allele_freq', type=float, default=1,
                            help="Maximum frequency an allele needs to have to be typed")
        parser.add_argument('--prop_typed_variants', type=float, default=1,
                            help="Proportion of variants that are typed (out of the ones that pass the frequency "
                                 "filter).")
        parser.add_argument('--trees_interval', nargs='+', type=int,
                            help="Only test the trees and variants in this interval for association")
        parser.add_argument('--trees_interval_start', type=int,
                            help="Only test the trees and variants starting with trees_interval_start and the orginial "
                                 "end for association")
        parser.add_argument('--trees_interval_end', type=int,
                            help="Only test the trees and variants starting with trees_interval_end and the orginial "
                                 "start for association")
        parser.add_argument('--chunk_size', type=int, help="Chunk size for task 'makeTreeChunks'")

        # simulating trees
        sim = parser.add_argument_group('simulating trees')
        sim.add_argument('--sim_tree_simulator', dest="sim_tree_simulator", default="stdPopsim",
                         choices=["stdPopsim", "msprime"],
                         help="Method used for simulating. stdPopsim is real human chromosome")
        sim.add_argument('--N', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim, or number of haploid "
                              "individuals per population when simulating with msprime")
        sim.add_argument('--N_African', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim")
        sim.add_argument('--N_European', type=int,
                         help="Number of haploid African individuals to simulate with stdPopsim")
        sim.add_argument('--mu', type=float,
                         help="Mutation rate for simulating with msprime (does not work for stdPopsim)")
        sim.add_argument('--population_size', type=int,
                         help="Ancestral population size for simulating with msprime", default=1000)
        sim.add_argument("--sim_two_populations", action="store_true",
                         help="Simulate a population split with msprime")
        sim.add_argument("--split_time", type=int,
                         help="Split time between two populations")
        sim.add_argument('--AH_tree_pos', type=int,
                         help="Add mutations only to the local tree covering this genomic position")
        sim.add_argument('--recomb_rate',
                         help="Recombination rate for simulating with msprime. Can be float or rate map file "
                              "(does not work for stdPopsim)")
        sim.add_argument('--recomb_map_start_random', action='store_true',
                         help="Use recombination rates starting at random position (does not work for stdPopsim)")
        sim.add_argument('--mut_window_size', type=int,
                         help="window size for mutation rate map for simulating with msprime (does not work for "
                              "stdPopsim)")
        sim.add_argument('--mut_beta_shape1', type=float, default=10,
                         help="shape1 parameter for window-based mutation rate beta distribution")
        sim.add_argument('--mut_beta_shape2', type=float, default=1000000000,
                         help="shape2 parameter for window-based mutation rate beta distribution")
        sim.add_argument('--sequence_length', type=float,
                         help="Sequence length for simulating with msprime (does not work for stdPopsim)")
        sim.add_argument('--pos_float', action='store_true',
                         help="Should the positions of the variants be transformed into integers. Msprime simulates a "
                              "continuous genome, so if pos_int is true, the simulated positions are rounded and if "
                              "one position overlaps the previous, it is moved to the next position in the genome.")

        # phenotypes
        pty = parser.add_argument_group('phenotypes')
        pty.add_argument('--pheno_file', type=str,
                         help="File with phenotypes")
        pty.add_argument('--pheno_file_BMI', type=str,
                         help="File with BMI phenotypes")
        pty.add_argument('--simulate_phenotypes', action='store_true',
                         help="If set to true, phenotypes will be simulated. If set to false, they will be read from "
                              "file.")
        pty.add_argument('--pty_sd_envNoise', type=float,
                         help="Std. dev. for environmental noise. If set to 0, no noise will be simulated.")
        pty.add_argument('--pty_mean_envNoise', type=float, default=0,
                         help="Mean for environmental noise. Should normally be zero")
        pty.add_argument('--pty_h_squared', type=float,
                         help="Heritability (h^2) of simulated phenotype. The environmental noise will be simulated in "
                              "order to satisfy the request.")
        pty.add_argument('--pty_sim_method',
                         choices=['null', 'uniform', 'fixed', 'singleTyped', 'singleUntyped', "allelicHetero",
                                  "oneTree", "oneRegion"],
                         help="Phenotype simulations method")
        pty.add_argument('--allow_typed_causal_variants', action='store_true',
                         help="By default only non-typed variants will be simulated to be causal with "
                              "pty_sim_method='oneRegion' or 'oneTree'")
        pty.add_argument('--pty_prop_causal_mutations', type=float, default=0,
                         help="Proportion of causal mutations to simulate at uniformly distributed positions if "
                              "pty_sim_method is set to 'uniform'. If set to 0, there will be no causal mutations.")
        pty.add_argument('--pty_sd_beta_causal_mutations', type=str,
                         help="Std. dev. for betas of causal mutations if pty_sim_method is set to 'uniform', "
                              "'oneTree' or 'oneWindow'. If it can be converted to a float, betas will sampled from "
                              "N(0, pty_sd_beta_causal_mutations). If set to 'freq_dependant', betas will be sampled "
                              "from N(0, [2 * f * (1 - f)]^{-0.5} * h2g / p), where h2g is the heritability of the "
                              "trait and p is the number of causal SNPs.")
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
        pty.add_argument('--causal_region_coordinates', type=float, nargs='+',
                         help="Coordinates of causal region (must be list of length 2). End is not included.")
        pty.add_argument('--min_allele_freq_causal', type=float, default=0.0,
                         help="Simulate phenotype with all variants in a region with this min allele freq to be causal")
        pty.add_argument('--max_allele_freq_causal', type=float, default=1.0,
                         help="Simulate phenotype with all variants in a region with this max allele freq to be causal")
        pty.add_argument('--add_1_to_half_of_inds', action='store_true',
                         help="Add 1 to half of all individuals (only makes sense when simulating with "
                              "sim_two_populations")
        pty.add_argument('--population_disease_prevalence', type=float,
                         help="The disease prevalence in the population. Specifying this parameter turns the "
                              "quantitative phenotype into a liability score and considers the top x percent of "
                              "individuals to have the disease")

        # run associations
        assoc = parser.add_argument_group('running association tests')
        assoc.add_argument('--ass_method', nargs='+', type=str,
                           help="Provide association method [AIM, GWAS] and for AIM covariance type [eGRM, GRM, "
                                "scaled] with the following format: 'method:covariance'")
        assoc.add_argument("--AIM_method", type=str, nargs='+',
                           help="Provide association algorithm to be used by AIM [GCTA_HE, GCTA_REML, glimix_REML] "
                                "or any combination of these to test the local GRM for association")
        assoc.add_argument('--no_clean_up', action='store_true',
                           help="don't clean up temporary files created during association testing")
        assoc.add_argument('--test_only_tree_at', type=float,
                           help="Only test tree that is overlapping the given position for association")
        assoc.add_argument('--skip_first_tree', action='store_true',
                           help='Do not run association test on first tree. This is useful when tskit file was '
                                'extracted from a larger set of trees (the first tree in that case will be empty)')
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
        assoc.add_argument('--ass_window_size', type=int,
                           help="window size for region-based association tests")
        # assoc.add_argument('--write_covariance_picklefiles', action='store_true',
        #                    help="write covariances from association tests to picklefile")
        assoc.add_argument('--covariance_picklefiles', type=str, nargs='+',
                           help="calculate correlation between covariances along the genome (user must verify that "
                                "windows have same coordinates)")
        assoc.add_argument('--population_structure_matrix', type=str,
                           help="Prefix of covariance matrix in GCTA format used to correct for population structure")
        assoc.add_argument('--population_structure_pca_num_eigenvectors', type=int,
                           help="Correct for population structure with this number of PCA eigenvectors computed with "
                                "gcta")
        assoc.add_argument('--PC_keep', type=int, nargs='+',
                           help="List of PCs to keep (1-based)")
        assoc.add_argument('--global_GRM_and_PCs_model', action='store_true',
                           help="Correct for stratification with global GRM provided with '--population_strucure' and "
                                "also with PCs")
        assoc.add_argument('--coreGREML_model', action='store_true',
                           help="Correct for stratification with global GRM provided with '--population_structure' and "
                                "also the correlation between the local and global GRMs (this is the coreGREML model of "
                                "Zhou et al. 2023)")
        assoc.add_argument('--GCTA', type=str, default="gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1",
                           help="path to GCTA executable")
        assoc.add_argument("--num_gcta_threads", type=int, default=2,
                           help="number of threads on which to run GCTA. Going from 2 to 8 roughly doubles the speed")
        assoc.add_argument("--additional_gcta_params", type=str, nargs='+',
                           help="additional arguments to be passed to GCTA")
        assoc.add_argument("--additional_mtg2_params", type=str, nargs='+',
                           help="additional arguments to be passed to mtg2")
        assoc.add_argument("--limit_association_tests", type=int, default=1000000000,
                           help="Limit number of association tests to be run. Applies to trees and windows (tree-based tests)")

        # associations with real data
        assoc.add_argument('--relate_sample_names', type=str,
                           help="file with sample names used to run Relate")

        args = parser.parse_args()

        return args
