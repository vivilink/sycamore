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
import TPhenotypes as pt
import association_functions as af
from python_log_indenter import IndentedLoggerAdapter
import logging
import os
import sys
import pandas as pd

os.chdir(os.path.dirname(sys.argv[0]))

# -----------------------------
# initialize arguments
# -----------------------------

TParams = params.TParameters()
args = TParams.initialize()

# -----------------------------
# initialize logfile
# -----------------------------

logger = logging.getLogger("indented")
file_handler = logging.FileHandler(args.out + ".log")
# logger.addHandler(logging.StreamHandler())
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

logger.info("--------------------->")
logger.info("sycamore 1.0")
logger.info("--------------------->")

# print arguments to logfile
# logger.info("- The following parameters were passed: " + str(args))
logger.info("- Writing output files with prefix '" + str(args.out) + "'")
# logger.info("- Adding plots to the following directory '" + str(args.out) + "_plots'")

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
    if args.sim_two_populations is True:
        trees = simulator.run_simulation_two_populations(arguments=args, randomGenerator=r, logfile=logger)
    else:
        trees = simulator.run_simulation(arguments=args, randomGenerator=r, logfile=logger)

    sample_ids = trees.samples()
    N = len(sample_ids)
    if args.N != N:
        logger.warning("WARNING: Number of samples in tree (" + str(N) + ") does not match number of samples in "
                                                                         "arguments (" + str(args.N) + ")")
    inds = tind.Individuals(ploidy=args.ploidy,
                            num_haplotypes=N,
                            relate_sample_names_file=args.relate_sample_names,
                            logfile=logger)
    variants = tvar.TVariants(tskit_object=trees, samp_ids=sample_ids)
    variants.fill_info(tskit_object=trees, samp_ids=sample_ids, pos_float=args.pos_float, logfile=logger)
    variants.write_variant_info(out=args.out, logfile=logger)
    variants.write_genetic_map_relate(out=args.out, logfile=logger)
    variants.write_genetic_map_argNeedle(out=args.out, logfile=logger, chrom=args.chromosome)

    tt.TTrees.writeStats(ts_object=trees, out=args.out, logfile=logger)

if args.task == "simulateMoreMutations":
    logger.info("- TASK: simulateMoreMutations")
    tsim.TSimulator.simulate_more_mutations(arguments=args, logfile=logger)

# -----------------------
# ARG statistics
# -----------------------
if args.task == "ARGStatistics":
    logger.info("- TASK: ARGStatistics")
    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)
    trees_object.writeStats(ts_object=trees_object.trees, out=args.out, logfile=logger)

# -----------------------
# Output single tree
# -----------------------
if args.task == "getTreeAtPosition":
    logger.info("- TASK: getTreeAtPosition")

    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)

    trees_object.extract_single_tree(trees_object.trees, args.out, logger, position=args.test_only_tree_at)

# -----------------------
# Output tree chunks
# -----------------------

if args.task == "makeTreeChunks":
    """
    Cut an ARG into chunks. Output tskit tree files with size 'chunk_size' 
    """
    logger.info("- TASK: makeTreeChunks")
    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)
    print("trees_object.actual_trees_interval in ARGWAS.py", trees_object.actual_trees_interval)

    if args.chunk_size is None:
        raise ValueError("Must provide window size")
    chunk_starts, chunk_ends = af.get_window_starts_and_ends(window_size=args.chunk_size,
                                                             trees_interval=trees_object.actual_trees_interval)
    num_windows = len(chunk_ends)

    for w in range(num_windows):
        trees_extract = trees_object.trees.keep_intervals([[chunk_starts[w], chunk_ends[w]]], simplify=True)
        print("trees_extract.num_trees", trees_extract.num_trees)
        print("coords of trees_extract", trees_extract.sequence_length)
        outname = args.out + "_chunk" + str(w) + ".trees"
        trees_extract.dump(outname)
        logger.info("- Wrote trees with coordinates [" + str(chunk_starts[w]) + "," + str(chunk_ends[w]) + ")"
                    + " to " + outname)

# -----------------------
# Downsample variants
# -----------------------
if args.task == "downsampleVariantsWriteShapeit":
    if args.prop_typed_variants is None:
        raise ValueError("Must provide downsampling probability to task 'downsampleVariantsWriteShapeit'")
    logger.info("- TASK: Downsampling variants")
    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)
    sample_ids = trees_object.trees.samples()
    N = len(sample_ids)

    # --------------------------------
    # create diploids and variants
    # --------------------------------
    inds = tind.Individuals(ploidy=args.ploidy,
                            num_haplotypes=N,
                            relate_sample_names_file=args.relate_sample_names,
                            logfile=logger)
    inds.write_shapeit2_relate(args.out, logger)
    inds.write_sample_argNeedle(args.out, logger)

    variants = tvar.TVariantsFiltered(trees_object.trees, sample_ids, args.min_allele_freq, args.max_allele_freq,
                                      args.prop_typed_variants, args.pos_float, r, logger)
    # variants = tvar.TVariantsFiltered(trees, samp_ids, 0.01, 1, 0.5, r)
    variants.write_variant_info(args.out, logger)
    variants.write_haps(args.out, inds, args.chromosome, logger)

    variants.write_genetic_map_relate(out=args.out, logfile=logger)
    variants.write_genetic_map_argNeedle(out=args.out, logfile=logger, chrom=args.chromosome)

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

    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)

    sample_ids = trees_object.trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(ploidy=args.ploidy,
                            num_haplotypes=N,
                            relate_sample_names_file=args.relate_sample_names,
                            logfile=logger)
    variants = tvar.TVariantsFiltered(tskit_object=trees_object, samp_ids=sample_ids,
                                      min_allele_freq=args.min_allele_freq,
                                      max_allele_freq=args.max_allele_freq,
                                      prop_typed_variants=args.prop_typed_variants, pos_float=args.pos_float, random=r,
                                      logfile=logger, filtered_variants_file=None)

    # impute
    imputation_obj = impute.TImpute()
    name_imputation_output = imputation_obj.run_impute(trees_sample=trees_object.trees,
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
# only simulate phenotypes
# ----------------------------------------------------------------

if args.task == "simulatePhenotypes":
    if args.tree_file is None:
        raise ValueError("The estimated trees need to be provided with 'tree_file'.")

    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logger)
    sample_ids = trees_object.trees.samples()

    plots_dir = args.out + "_plots/"
    # if not os.path.exists(plots_dir):
    #     os.mkdir(plots_dir)
    N = len(sample_ids)
    inds = tind.Individuals(ploidy=args.ploidy, num_haplotypes=N, relate_sample_names_file=args.relate_sample_names,
                            logfile=logger)

    pheno = pt.make_phenotypes(args=args,
                               trees=trees_object.trees,
                               sample_ids=sample_ids,
                               inds=inds,
                               plots_dir=plots_dir,
                               random=r,
                               logfile=logger)

    if args.population_disease_prevalence:
        logger.add()
        logger.info(
            "- Adding binary disease status with a population prevalence of " + str(args.population_disease_prevalence))
        pheno.add_disease_status(prevalence=args.population_disease_prevalence, logfile=logger)
        pheno.write_to_file_gcta_eGRM_disease_status(inds=inds, out=args.out, logfile=logger)

    pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logger)
    logger.sub()


# ----------------------------------------------------------------
# Read simulation to simulate phenotypes and perform association
# ----------------------------------------------------------------

if args.task == "associate":
    af.run_association_testing(args=args, random=r, logfile=logger)

# ----------------------------------------------------------------
# Create case control sample
# ----------------------------------------------------------------

if args.task == "transformToBinaryAndAscertain":
    pt.transformAndAscertain(args=args, random=r, logger=logger)

    # print(pheno_file['3'].iloc[pheno_file['3'] == 0.0])

# # ----------------------------------------------------------------
# # Create case control sample
# # ----------------------------------------------------------------
#
# if args.task == "ascertainCaseControlSample":
#     if args.disease_status_file is None:
#         raise ValueError("Provide disease status file with 'disease_status_file'")
#     if args.sample_size is None:
#         raise ValueError("Provide expected sample size with 'sample_size'")
#
#     # read file and find population size
#     pheno_file = pd.read_csv(args.disease_status_file, delim_whitespace=True, header=None)
#     pheno_file.columns = ['1', '2', '3']
#     population_size = pheno_file.shape[0]
#     if population_size < args.sample_size:
#         raise ValueError("Sample size cannot be larger than number of phenotypes in phen file (" + str(population_size) + ")")
#
#     # find cases
#     cases = pheno_file.loc[pheno_file['3'] == 1.0, :]
#     num_cases = cases.shape[0]
#
#     # find controls
#     num_controls = args.sample_size - num_cases
#     if num_controls <= 0:
#         raise ValueError("Sample size needs to be larger than number of cases")
#     logger.info("- Keeping " + str(num_cases) + " cases and " + str(num_controls) + " controls from a total of " +
#                 str(population_size) + " individuals")
#     controls = pheno_file.loc[pheno_file['3'] == 0.0, :].sample(n=num_controls)
#
#     # concatenate and write
#     pheno_file_new = pd.concat([cases, controls])
#     logger.info("- Writing ascertained sample's disease status to '" + args.out + "_disease_status_ascertained.phen'")
#     pheno_file_new.to_csv(args.out + "_disease_status_ascertained.phen", sep=' ', index=False, header=False)
#
#     # print(pheno_file['3'].iloc[pheno_file['3'] == 0.0])

# ----------------------------------------------------------------
# Write variants plink files
# ----------------------------------------------------------------

if args.task == "writeToPlink":
    raise ValueError("task 'writeToPlink' is planned but not yet implemented")

# ----------------------------------------------------------------
# Write variants plink files
# ----------------------------------------------------------------
logger.info("<---------------------")
logger.info("sycamore done")
logger.info("<---------------------")

# if __name__ == "__main__":
#     sys.exit(main(sys.argv[1:]))
