#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 16:57:18 2021

@author: linkv
"""
import tskit
import TParameters as params
import TPhenotypes as pt
import TAssociationTesting as gwas
import TVariants as tvar
import TIndividuals as tind
import TSimulator as tsim
import TTree as tt
from python_log_indenter import IndentedLoggerAdapter
import logging
import os
import sys
from numpy.random import RandomState

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
class TRandomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)


r = TRandomGenerator(args.seed)

logger.info("- randomGenerator seed is set to " + str(r.random.get_state()[1][0]))

# -----------------------
# Simulate
# -----------------------

if args.task == "simulate":
    logger.info("- TASK: simulate")

    if args.N is None:
        raise ValueError("Must provide sample size with argument '--N'")

    if args.sim_tree_simulator == "stdPopsim":
        simulator = tsim.TSimulatorStdPopsim()
    elif args.sim_tree_simulator == "msprime":
        simulator = tsim.TSimulatorMSPrime()
    else:
        logger.error("use of any simulator besides stdPopSim not tested")
        raise ValueError("use of any simulator besides stdPopSim not tested")

    trees = simulator.run_simulation(arguments=args, randomGenerator=r, logfile=logger)

    samp_ids = trees.samples()
    N = len(samp_ids)
    if args.N != N:
        logger.warning("Number of samples in tree does not match number of samples in arguments")
    inds = tind.Individuals(args.ploidy, N)
    variants = tvar.TVariants(ts_object=trees, samp_ids=samp_ids, pos_int=args.pos_int)
    variants.write_variant_info(ts_object=trees, samp_ids=samp_ids, out=args.out, logfile=logger)

    tt.TTrees.writeStats(ts_object=trees, out=args.out, logfile=logger)

if args.task == "simulateMoreMutations":
    logger.info("- TASK: simulateMoreMutations")
    tsim.TSimulator.simulate_more_mutations(arguments=args, logfile=logger)

# -----------------------
# ARG statistics
# -----------------------
if args.task == "ARGStatistics":
    logger.info("- TASK: ARGStatistics")
    logger.info("- Reading tree from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    trees_class = tt.TTrees(ts_object=trees)
    trees_class.writeStats(ts_object=trees, out=args.out, logfile=logger)

# -----------------------
# Output single tree
# -----------------------
if args.task == "getTreeAtPosition":
    logger.info("- TASK: getTreeAtPosition")
    logger.info("- Reading tree from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    trees_class = tt.TTrees(trees)
    trees_class.extract_single_tree(trees, args.out, logger, position=args.test_only_tree_at)  # 49027865

# -----------------------
# Downsample variants
# -----------------------
if args.task == "downsampleVariants":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariant'")
    logger.info("- TASK: Downsampling variants")
    logger.info("- Reading tree simulations from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    samp_ids = trees.samples()

    # --------------------------------
    # create diploids and variants
    # --------------------------------
    inds = tind.Individuals(args.ploidy, args.N)
    inds.write_shapeit2(args.out, logger)
    variants = tvar.TVariantsFiltered(trees, samp_ids, args.min_allele_freq, args.max_allele_freq,
                                      args.prop_typed_variants, args.pos_int, r, logger)
    # variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, 0.5, r)
    variants.write_variant_info(args.out, logger)
    variants.write_shapeit2(args.out, inds, logger)

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
    trees = tskit.load(args.tree_file)

    if trees_orig.num_samples != trees.num_samples:
        raise ValueError(
            "The trees provided with params 'tree_file_simulated' and 'tree_file' must have the same number of samples")

    samp_ids = trees.samples()
    N = len(samp_ids)

    # --------------------------------
    # create diploids and variants
    # --------------------------------

    inds = tind.Individuals(args.ploidy, N)

    # TODO: trees_orig and variants_orig should be initialized at the same time, e.g. together in one function. We should not have 2 tree files and 2 variant files just floating around separately

    #  TODO: find way to save variants in their tskit format without needing to read the original tree. I only need original tree in
    #   association task for this. It would be nice if the only tree that needs to be read would be estimated tree do
    #   not provide variant file here but have it estimated from tree, otherwise variants and tree won't match (tree
    #   only contains typed variants). The variant file is only useful for simulating phenotypes to be able to keep
    #   track of untyped variants
    variants = tvar.TVariantsFiltered(ts_object=trees, samp_ids=samp_ids, min_allele_freq=args.min_allele_freq,
                                      max_allele_freq=args.max_allele_freq,
                                      prop_typed_variants=args.prop_typed_variants, pos_int=args.pos_int, random=r,
                                      logfile=logger)

    # variants_orig are used to simulate phenotypes. They need to be consistent with original tree and the typed
    # status that might have been defined earlier with a variants file. The causal mutation should not be affected by
    # a freq filter
    variants_orig = tvar.TVariantsFiltered(ts_object=trees_orig, samp_ids=samp_ids, min_allele_freq=0,
                                           max_allele_freq=1, prop_typed_variants=1, pos_int=args.pos_int, random=r,
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

    if args.ass_method == "both":
        logger.info("- Running both GWAS and AIM for associating")
    else:
        logger.info("- Running " + args.ass_method + " for associating")

    if args.ass_method == "GWAS" or args.ass_method == "both":
        logger.info("- GWAS:")
        logger.add()
        GWAS = gwas.TAssociationTesting_GWAS(phenotypes=pheno, num_typed_variants=variants.number_typed)

        GWAS.OLS(variants, inds, logger)
        GWAS.write_to_file(variants, args.out, logger)
        GWAS.manhattan_plot(variant_positions=variants.info['position'], plots_dir=plots_dir)

        logger.sub()

    if args.ass_method == "AIM" or args.ass_method == "both":

        logger.info("- AIM:")
        logger.add()

        if args.AIM_method is None:
            raise ValueError("ERROR: No method for tree association provided. Use '--AIM_method' to set method.")
        if args.covariance_type is None:
            raise ValueError(
                "ERROR: No method for covariance calculation provided. Use '--covariance_type' to set method.")

        logger.info("- Reading tree estimations for tree-based association from " + args.tree_file)

        pheno.find_causal_trees(trees)

        for m in args.AIM_method:

            logger.add()

            if m == "HE":
                treeWAS = gwas.TAssociationTesting_trees_gcta_HE(trees, pheno)

                # write phenotypes in gcta format
                if args.covariance_type == "eGRM" or args.covariance_type == "GRM":
                    pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logger)
                else:
                    pheno.write_to_file_gcta_scaled(out=args.out, logfile=logger)

                    # run association
                if args.test_only_tree_at is None:
                    logger.info("- Running associations test using GCTA Haseman-Elston for a sequence of trees")
                    logger.add()
                    treeWAS.run_association(ts_object=trees, variants=variants, inds=inds, out=args.out, logfile=logger,
                                            covariance_type=args.covariance_type, skip_first_tree=args.skip_first_tree)
                    logger.sub()
                else:
                    logger.info("- Running associations test using GCTA Haseman-Elston for a single tree")
                    logger.add()
                    tree = trees.at(args.test_only_tree_at)
                    tree_obj = tt.TTree(tree)
                    treeWAS.run_association_one_tree(ts_object=trees, variants=variants, tree_obj=tree_obj, inds=inds,
                                                     out=args.out, logfile=logger, covariance_type=args.covariance_type,
                                                     skip_first_tree=args.skip_first_tree)
                    logger.sub()

                treeWAS.write_to_file(trees, args.out, logger)

            if m == "REML":
                treeWAS = gwas.TAssociationTesting_trees_gcta_REML(trees, pheno)

                # write phenotypes in gcta format
                if args.covariance_type == "eGRM" or args.covariance_type == "GRM":
                    pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logger)
                else:
                    pheno.write_to_file_gcta_scaled(out=args.out, logfile=logger)

                # run association
                if args.test_only_tree_at is None:
                    logger.info("- Running associations test using GCTA REML for a sequence of trees")
                    logger.add()
                    treeWAS.run_association(ts_object=trees, variants=variants, inds=inds, out=args.out, logfile=logger,
                                            covariance_type=args.covariance_type, skip_first_tree=args.skip_first_tree)
                    logger.sub()
                else:
                    logger.info("- Running associations test using GCTA REML for a single tree")
                    logger.add()
                    tree = trees.at(args.test_only_tree_at)
                    tree_obj = tt.TTree(tree)
                    treeWAS.run_association_one_tree(ts_object=trees, variants=variants, tree_obj=tree_obj, inds=inds,
                                                     out=args.out, logfile=logger, covariance_type=args.covariance_type,
                                                     skip_first_tree=args.skip_first_tree)
                    logger.sub()

                treeWAS.write_to_file(trees, args.out, logger)

                logger.sub()

            logger.sub()

        logger.sub()

        logger.info("- Done running AIM")
