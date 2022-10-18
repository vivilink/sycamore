#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
import tskit
import TParameters as params
import TRandomGenerator as rg
import TPhenotypes as pt
import TAssociationTesting as at
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import TTree as tt
import TImputation as impute
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

plots_dir = args.out + "_plots/"
if not os.path.exists(plots_dir):
    os.mkdir(plots_dir)

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


# ---------------------
# read tree file
# ---------------------
def read_trees(tree_file, trees_interval, trees_interval_start, trees_interval_end):
    logger.info("- Reading tree from " + tree_file)
    trees_full = tskit.load(tree_file)
    # if args.trees_interval is not None and (args.trees_interval_start is not None or args.trees_interval_end is not None):
    #     raise ValueError("Cannot set 'trees_interval' and 'trees_interval_start' or 'trees_interval_end' at the same "
    #                      "time.")
    if trees_interval is None:
        trees_interval = [0, trees_full.sequence_length]
    if trees_interval_start:
        trees_interval[0] = trees_interval_start
    if trees_interval_end:
        trees_interval[1] = trees_interval_end

    logger.info(
        "- Running association only on the trees overlapping the following interval: " + str(trees_interval))
    trees_extract = trees_full.keep_intervals([trees_interval], simplify=True)
    args.trees_interval = trees_interval

    return trees_extract, trees_interval


# -----------------------
# ARG statistics
# -----------------------
if args.task == "ARGStatistics":
    logger.info("- TASK: ARGStatistics")
    trees, args.trees_interval = read_trees(tree_file=args.tree_file,
                                            trees_interval=args.trees_interval,
                                            trees_interval_start=args.trees_interval_start,
                                            trees_interval_end=args.trees_interval_end)

    trees_class = tt.TTrees(ts_object=trees)
    trees_class.writeStats(ts_object=trees, out=args.out, logfile=logger)

# -----------------------
# Output single tree
# -----------------------
if args.task == "getTreeAtPosition":
    logger.info("- TASK: getTreeAtPosition")
    trees, args.trees_interval = read_trees(tree_file=args.tree_file,
                                            trees_interval=args.trees_interval,
                                            trees_interval_start=args.trees_interval_start,
                                            trees_interval_end=args.trees_interval_end)
    trees_class = tt.TTrees(trees)
    trees_class.extract_single_tree(trees, args.out, logger, position=args.test_only_tree_at)

# -----------------------
# Downsample variants
# -----------------------
if args.task == "downsampleVariantsWriteShapeit":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariantsWriteShapeit'")
    logger.info("- TASK: Downsampling variants")
    trees, args.trees_interval = read_trees(tree_file=args.tree_file,
                                            trees_interval=args.trees_interval,
                                            trees_interval_start=args.trees_interval_start,
                                            trees_interval_end=args.trees_interval_end)
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

    logger.info(
        "-  Imputing genotypes from " + args.tree_file + " with reference panel trees from " + args.imputation_ref_panel_tree_file)

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

# -----------------------
# Write variants plink files
# -----------------------

if args.task == "writeToPlink":
    raise ValueError("task 'writeToPlink' is planned but not yet implemented")

# ----------------------------------------------------------------
# Read simulation to simulate phenotypes and perform association
# ----------------------------------------------------------------

if args.task == "associate":

    if args.tree_file_simulated is None or args.tree_file is None:
        raise ValueError(
            "Both the simulated and estimated trees need to be provided with 'tree_file_simulated' and 'tree_file'.")

    logger.info("- TASK: Associate")
    logger.info("- Reading simulated tree used for simulating phenotypes from " + args.tree_file_simulated)
    trees_orig = tskit.load(args.tree_file_simulated)

    logger.info(
        "- Reading tree used for Aim's association testing, and for defining variants to be tested by GWAS, from " + args.tree_file)
    trees, args.trees_interval = read_trees(tree_file=args.tree_file,
                                            trees_interval=args.trees_interval,
                                            trees_interval_start=args.trees_interval_start,
                                            trees_interval_end=args.trees_interval_end)

    if trees_orig.num_samples != trees.num_samples:
        raise ValueError(
            "The trees provided with params 'tree_file_simulated' and 'tree_file' must have the same number of samples")

    # --------------------------------
    # create diploids and variants
    # --------------------------------
    sample_ids = trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(args.ploidy, N)
    trees = tt.TTrees.remove_monomorphic(trees)

    # TODO: trees_orig and variants_orig should be initialized at the same time, e.g. together in one function. We
    #  should not have 2 tree files and 2 variant files just floating around separately. Maybe TTrees class could
    #  contain the tskit sequence and variants

    # TODO: find way to save variants in their tskit format without needing to read the original tree. I only need
    #  original tree in association task for this. It would be nice if the only tree that needs to be read would be
    #  estimated tree

    #  do not provide variant file here but have it estimated from tree, otherwise variants and tree
    #  won't match (tree only contains typed variants). The variant file is only useful for simulating phenotypes to
    #  be able to keep track of untyped variants (i.e. for variants_orig)
    variants = tvar.TVariantsFiltered(ts_object=trees, samp_ids=sample_ids, min_allele_freq=args.min_allele_freq,
                                      max_allele_freq=args.max_allele_freq,
                                      prop_typed_variants=args.prop_typed_variants, pos_float=args.pos_float, random=r,
                                      logfile=logger, filtered_variants_file=None)
    variants.write_variant_info(out=args.out + "_sample", logfile=logger)

    # variants_orig are used to simulate phenotypes. They need to be consistent with original tree and the typed
    # status that might have been defined earlier with a variants file. The causal mutation should not be affected by
    # a freq filter
    variants_orig = tvar.TVariantsFiltered(ts_object=trees_orig, samp_ids=sample_ids, min_allele_freq=0,
                                           max_allele_freq=1, prop_typed_variants=1, pos_float=args.pos_float, random=r,
                                           logfile=logger, filtered_variants_file=args.variants_file)

    # --------------------------------
    # create phenotypes
    # --------------------------------
    logger.info("- Phenotypes:")
    logger.add()

    pheno = pt.Phenotypes(variants=variants_orig, inds=inds, logfile=logger)

    pheno.simulate(args=args, r=r, logfile=logger, variants_orig=variants_orig, trees=trees, inds=inds,
                   plots_dir=plots_dir)
    logger.sub()

    # --------------------------------
    # run association tests and plot
    # --------------------------------
    covariance_types = []
    for m in args.ass_method:
        method = m.split(':')[0]
        logger.info("- Running " + m + " for associating")

        logger.add()

        if method == "GWAS":
            at.run_association_GWAS(trees=trees, inds=inds, variants=variants, pheno=pheno, args=args, impute=impute,
                                    logfile=logger)

        elif method == "AIM":
            at.run_association_AIM(trees=trees,
                                   inds=inds,
                                   variants=variants,
                                   pheno=pheno,
                                   args=args,
                                   ass_method=m,
                                   window_size=args.ass_window_size,
                                   logfile=logger)

        else:
            raise ValueError("Unknown association method '" + m + "'. Must be 'AIM' or 'GWAS'.")

        logger.sub()

    logger.info("- Done running association tests")

# if __name__ == "__main__":
#     sys.exit(main(sys.argv[1:]))
