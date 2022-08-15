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
import TAssociationTesting as gwas
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
def read_trees():
    logger.info("- Reading tree from " + args.tree_file)
    trees = tskit.load(args.tree_file)
    if args.trees_interval is not None:
        logger.info(
            "- Running association only on the trees overlapping the following interval: " + str(args.trees_interval))
        trees = trees.keep_intervals([args.trees_interval], simplify=True)
    else:
        args.trees_interval = [0, trees.sequence_length]

    return trees

# -----------------------
# ARG statistics
# -----------------------
if args.task == "ARGStatistics":
    logger.info("- TASK: ARGStatistics")
    trees = read_trees()
    trees_class = tt.TTrees(ts_object=trees)
    trees_class.writeStats(ts_object=trees, out=args.out, logfile=logger)

# -----------------------
# Output single tree
# -----------------------
if args.task == "getTreeAtPosition":
    logger.info("- TASK: getTreeAtPosition")
    trees = read_trees()
    trees_class = tt.TTrees(trees)
    trees_class.extract_single_tree(trees, args.out, logger, position=args.test_only_tree_at)

# -----------------------
# Downsample variants
# -----------------------
if args.task == "downsampleVariantsWriteShapeit":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariantsWriteShapeit'")
    logger.info("- TASK: Downsampling variants")
    trees = read_trees()
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
    trees = read_trees()

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

    if args.ass_method == "both":
        logger.info("- Running both GWAS and AIM for associating")
    else:
        logger.info("- Running " + args.ass_method + " for associating")

    if args.ass_method == "GWAS" or args.ass_method == "both":
        logger.info("- GWAS:")
        logger.add()

        if args.imputation_ref_panel_tree_file is not None:
            logger.info("- Using genotypes imputed with impute2 for GWAS:")
            logger.add()

            if args.do_imputation:
                # impute
                imputation_obj = impute.TImpute()
                name_imputation_output = imputation_obj.run_impute(trees_sample=trees,
                                                                   variants_sample=variants,
                                                                   inds=inds,
                                                                   imputation_ref_panel_tree_file=args.imputation_ref_panel_tree_file,
                                                                   ploidy_ref=args.ploidy_ref,
                                                                   genetic_map_file=args.genetic_map_file,
                                                                   out=args.out, logfile=logger)
            else:
                if args.imputed_gen_file is None:
                    raise ValueError(
                        "When --do_imputation is set to False, the imputed genotypes must be provided "
                        "with --imputed_gen_file parameter")

                logger.info(
                    "- Assuming imputation was already run, reading imputed genotypes from " + args.imputed_gen_file)
                name_imputation_output = args.imputed_gen_file

            # read imputed genotypes
            gt_matrix_imputed, pos = impute.TImpute.read_imputed_gt(name_imputation_output=name_imputation_output,
                                                                    variants_sample=variants,
                                                                    trees_interval=args.trees_interval,
                                                                    logfile=logger)
            logger.sub()

            # run association tests
            GWAS = gwas.TAssociationTesting_GWAS(phenotypes=pheno, num_typed_variants=gt_matrix_imputed.shape[1])
            GWAS.test_with_positions_from_X_matrix(X=gt_matrix_imputed, positions=pos,
                                                   variants_sample=variants,
                                                   logfile=logger)
            GWAS.write_to_file_with_X_matrix(positions=pos, name=args.out, logfile=logger)

        else:
            logger.info("- Using genotypes from tree file for GWAS:")
            # run association tests
            GWAS = gwas.TAssociationTesting_GWAS(phenotypes=pheno, num_typed_variants=variants.num_typed)
            GWAS.test_with_variants_object(variants, inds, logger)
            GWAS.write_to_file(variants, args.out, logger)
            # GWAS.manhattan_plot(variant_positions=variants.info['position'], plots_dir=plots_dir)

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

            if m == "HE":

                if args.test_only_tree_at is None:
                    logger.info("- Running associations test using GCTA Haseman-Elston for a sequence of trees")
                else:
                    logger.info("- Running associations test using GCTA Haseman-Elston for a single tree")
                logger.add()
                treeWAS = gwas.TAssociationTesting_trees_gcta_HE(trees, pheno)

                # write phenotypes in gcta format
                if args.covariance_type == "eGRM" or args.covariance_type == "GRM":
                    pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logger)
                else:
                    pheno.write_to_file_gcta_scaled(out=args.out, logfile=logger)

                # run association
                if args.test_only_tree_at is None:
                    treeWAS.run_association(ts_object=trees, variants=variants, inds=inds, out=args.out, logfile=logger,
                                            covariance_type=args.covariance_type, skip_first_tree=args.skip_first_tree)
                else:
                    tree = trees.at(args.test_only_tree_at)
                    tree_obj = tt.TTree(tree)
                    treeWAS.run_association_one_tree(ts_object=trees, variants=variants, tree_obj=tree_obj, inds=inds,
                                                     out=args.out, logfile=logger, covariance_type=args.covariance_type,
                                                     skip_first_tree=args.skip_first_tree)

                treeWAS.write_to_file(trees, args.out, logger)
                logger.sub()

            if m == "REML":

                if args.test_only_tree_at is None:
                    logger.info("- Running associations test using GCTA REML for a sequence of trees")
                else:
                    logger.info("- Running associations test using GCTA REML for a single tree")
                logger.add()
                treeWAS = gwas.TAssociationTesting_trees_gcta_REML(trees, pheno)

                # write phenotypes in gcta format
                if args.covariance_type == "eGRM" or args.covariance_type == "GRM":
                    pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logger)
                else:
                    pheno.write_to_file_gcta_scaled(out=args.out, logfile=logger)

                # run association
                if args.test_only_tree_at is None:
                    treeWAS.run_association(ts_object=trees, variants=variants, inds=inds, out=args.out, logfile=logger,
                                            covariance_type=args.covariance_type, skip_first_tree=args.skip_first_tree)
                else:
                    tree = trees.at(args.test_only_tree_at)
                    tree_obj = tt.TTree(tree)
                    treeWAS.run_association_one_tree(ts_object=trees, variants=variants, tree_obj=tree_obj, inds=inds,
                                                     out=args.out, logfile=logger, covariance_type=args.covariance_type,
                                                     skip_first_tree=args.skip_first_tree)

                treeWAS.write_to_file(trees, args.out, logger)

                logger.sub()

        logger.sub()

        logger.info("- Done running association tests")


# if __name__ == "__main__":
#     sys.exit(main(sys.argv[1:]))