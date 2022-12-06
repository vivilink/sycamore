#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
import tskit
import TParameters as params
import TRandomGenerator as rg
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import TTree as tt
import TImputation as impute
import association_functions as af
from python_log_indenter import IndentedLoggerAdapter
import logging
import os
import sys

os.chdir(os.path.dirname(sys.argv[0]))

# -----------------------------
# initialize arguments
# -----------------------------

TParams = params.TParameters()
args = TParams.initialize()

# -----------------------------
# initialize logfile
# -----------------------------

logger = logging.getLogger()
file_handler = logging.FileHandler(args.out + ".log")
logger.addHandler(logging.StreamHandler())
logger.addHandler(file_handler)

logging.basicConfig(
    format='%(asctime)s %(levelname)s %(message)s',
    filename=args.out + ".log",
    filemode='w'  # this makes the log file be not append I think
)
logger = IndentedLoggerAdapter(logger)
logger.setLevel(logging.INFO)

# -----------------------------
# Header
# -----------------------------

logger.info("---------------------")
logger.info("AIM 0.1")
logger.info("---------------------")

# print arguments to logfile
logger.info("- The following parameters were passed: " + str(args))
logger.info("- Writing output files with prefix '" + str(args.out) + "'")
logger.info("- Adding plots to the following directory '" + str(args.out) + "_plots'")

# -----------------------------
# initialize random generator
# -----------------------------

r = rg.TRandomGenerator(args.seed)

logger.info("- randomGenerator seed is set to " + str(r.random.get_state()[1][0]))

# -----------------------
# Simulate
# -----------------------

if args.task == "simulate":
    logger.info("- TASK: simulate")

    if args.N is None:
        raise ValueError("Must provide sample size with argument '--N'")

    logger.info("- Simulating " + str(args.N) + " individuals with ploidy " + str(args.ploidy))

    if args.sim_tree_simulator == "stdPopsim":
        simulator = tsim.TSimulatorStdPopsim()
    elif args.sim_tree_simulator == "msprime":
        simulator = tsim.TSimulatorMSPrime()
    else:
        logger.error("use of any simulator besides stdPopSim not tested")
        raise ValueError("use of any simulator besides stdPopSim not tested")
    # TODO: trees should be put in a wrapper class, only this class should import tskit
    trees = simulator.run_simulation(arguments=args, randomGenerator=r, logfile=logger)

    sample_ids = trees.samples()
    N = len(sample_ids)
    if args.N != N:
        logger.warning("WARNING: Number of samples in tree (" + str(N) + ") does not match number of samples in "
                                                                         "arguments (" + str(args.N) + ")")
    inds = tind.Individuals(args.ploidy, N)
    variants = tvar.TVariants(ts_object=trees, samp_ids=sample_ids)
    variants.fill_info(ts_object=trees, samp_ids=sample_ids, pos_float=args.pos_float, logfile=logger)
    variants.write_variant_info(out=args.out, logfile=logger)

    tt.TTrees.writeStats(ts_object=trees, out=args.out, logfile=logger)

if args.task == "simulateMoreMutations":
    logger.info("- TASK: simulateMoreMutations")
    tsim.TSimulator.simulate_more_mutations(arguments=args, logfile=logger)

# -----------------------
# ARG statistics
# -----------------------
if args.task == "ARGStatistics":
    logger.info("- TASK: ARGStatistics")
    trees, args.trees_interval = tt.read_trees(tree_file=args.tree_file,
                                               trees_interval=args.trees_interval,
                                               trees_interval_start=args.trees_interval_start,
                                               trees_interval_end=args.trees_interval_end,
                                               logfile=logger)

    trees_class = tt.TTrees(ts_object=trees)
    trees_class.writeStats(ts_object=trees, out=args.out, logfile=logger)

# -----------------------
# Output single tree
# -----------------------
if args.task == "getTreeAtPosition":
    logger.info("- TASK: getTreeAtPosition")
    trees, args.trees_interval = tt.read_trees(tree_file=args.tree_file,
                                               trees_interval=args.trees_interval,
                                               trees_interval_start=args.trees_interval_start,
                                               trees_interval_end=args.trees_interval_end,
                                               logfile=logger)
    trees_class = tt.TTrees(trees)
    trees_class.extract_single_tree(trees, args.out, logger, position=args.test_only_tree_at)

# -----------------------
# Downsample variants
# -----------------------
if args.task == "downsampleVariantsWriteShapeit":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariantsWriteShapeit'")
    logger.info("- TASK: Downsampling variants")
    trees, args.trees_interval = tt.read_trees(tree_file=args.tree_file,
                                               trees_interval=args.trees_interval,
                                               trees_interval_start=args.trees_interval_start,
                                               trees_interval_end=args.trees_interval_end,
                                               logfile=logger)
    sample_ids = trees.samples()
    N = len(sample_ids)

    # --------------------------------
    # create diploids and variants
    # --------------------------------
    inds = tind.Individuals(args.ploidy, N)
    inds.write_shapeit2(args.out, logger)
    variants = tvar.TVariantsFiltered(trees, sample_ids, args.min_allele_freq, args.max_allele_freq,
                                      args.prop_typed_variants, args.pos_float, r, logger)
    # variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, 0.5, r)
    variants.write_variant_info(args.out, logger)
    variants.write_shapeit2(args.out, inds, logger)

    variants.write_genetic_map(out=args.out, logfile=logger)

# -----------------------
# Impute
# -----------------------
if args.task == "impute":
    logger.info("- TASK: Impute")
    if args.imputation_ref_panel_tree_file is None:
        raise ValueError("Most provide reference panel tree file using --imputation_ref_panel_tree_file")

    logger.info("-  Imputing genotypes from " + args.tree_file
                + " with reference panel trees from "
                + args.imputation_ref_panel_tree_file)

    # TODO: this needs to be standardized and used in all tasks
    trees = tskit.load(args.tree_file)
    sample_ids = trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(args.ploidy, N)
    trees = tt.TTrees.remove_monomorphic(trees)
    variants = tvar.TVariantsFiltered(ts_object=trees, samp_ids=sample_ids, min_allele_freq=args.min_allele_freq,
                                      max_allele_freq=args.max_allele_freq,
                                      prop_typed_variants=args.prop_typed_variants, pos_float=args.pos_float, random=r,
                                      logfile=logger, filtered_variants_file=None)

    # impute
    imputation_obj = impute.TImpute()
    name_imputation_output = imputation_obj.run_impute(trees_sample=trees,
                                                       variants_sample=variants,
                                                       inds=inds,
                                                       imputation_ref_panel_tree_file=args.imputation_ref_panel_tree_file,
                                                       ploidy_ref=args.ploidy_ref,
                                                       genetic_map_file=args.genetic_map_file,
                                                       out=args.out, logfile=logger)

# ----------------------------------------------------------------
# Compare covariance matrices
# ----------------------------------------------------------------

if args.task == "covarianceCorrelations":
    af.run_covariance_correlation(args=args, logfile=logger)

# ----------------------------------------------------------------
# Read simulation to simulate phenotypes and perform association
# ----------------------------------------------------------------

if args.task == "associate":
    af.run_association_testing(args=args, random=r, logfile=logger)

# -----------------------
# Write variants plink files
# -----------------------

if args.task == "writeToPlink":
    raise ValueError("task 'writeToPlink' is planned but not yet implemented")

# if __name__ == "__main__":
#     sys.exit(main(sys.argv[1:]))
