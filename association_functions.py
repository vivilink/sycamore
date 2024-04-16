#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 17:37:12 2022
@author: linkv
"""

import numpy as np
import statsmodels.api as sm
import pickle
import os
import tskit
import TTree as tt
import time
import TVariants as tvar
import TCovariance as cov
import TIndividuals as tind
import TPhenotypes as pt
import TAssociationTesting as at
import TImputation as impute
import glob
import TParameters
from typing import IO
from python_log_indenter import IndentedLoggerAdapter


def run_association_testing(args, random, logfile):
    """
    Simulate phenotypes and perform association

    @param args: TParameters
    @param imputation_obj: TImputation
    @param random: TRandom
    @param logfile:
    @return:
    """

    logfile.info("- TASK: Associate")
    # logfile.add()

    # --------------------------------
    # read ARG
    # --------------------------------
    if args.tree_file is None:
        raise ValueError("The estimated trees need to be provided with 'tree_file'.")

    trees_object = tt.TTrees(tree_file=args.tree_file,
                             trees_interval=args.trees_interval,
                             trees_interval_start=args.trees_interval_start,
                             trees_interval_end=args.trees_interval_end,
                             skip_first_tree=args.skip_first_tree,
                             logfile=logfile)

    trees = trees_object.trees

    # interval = np.array([173584350,175000000])
    # trees_tmp = trees.keep_intervals([interval], simplify=True)
    # trees_tmp.dump("/data/ARGWAS/hawaiians/chr5.part-05_1mb.trees")

    # plots_dir = args.out + "_plots/"
    # if not os.path.exists(plots_dir):
    #     os.mkdir(plots_dir)

    sample_ids = trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(ploidy=args.ploidy, num_haplotypes=N, relate_sample_names_file=args.relate_sample_names,
                            logfile=logfile)

    # -----------------------------------------
    # create variants object from tskit trees
    # -----------------------------------------

    #  do not provide variant file here but have it estimated from tree, otherwise variants and tree
    #  won't match (tree only contains typed variants). The variant file is only useful for simulating phenotypes to
    #  be able to keep track of untyped variants (i.e. for variants_orig)
    logfile.info("- Compiling variant info from trees file")
    variants = tvar.TVariantsFiltered(tskit_object=trees,
                                      samp_ids=sample_ids,
                                      min_allele_freq=args.min_allele_freq,
                                      max_allele_freq=args.max_allele_freq,
                                      prop_typed_variants=args.prop_typed_variants,
                                      pos_float=args.pos_float,
                                      random=random,
                                      logfile=logfile,
                                      filtered_variants_file=None)
    variants.write_variant_info(out=args.out + "_sample", logfile=logfile)

    # --------------------------------
    # create phenotypes
    # --------------------------------

    pheno = pt.make_phenotypes(args=args,
                               trees=trees,
                               sample_ids=sample_ids,
                               inds=inds,
                               plots_dir="",
                               random=random,
                               logfile=logfile)

    # --------------------------------
    # run association tests and plot
    # --------------------------------
    if args.ass_method is None:
        raise ValueError("Must provide association method using 'ass_method'")
    for m in args.ass_method:
        method = m.split(':')[0]
        logfile.info("- Running " + m + " for associating")

        logfile.add()

        if method == "GWAS":
            run_association_GWAS(trees=trees,
                                 samp_ids=sample_ids,
                                 inds=inds,
                                 variants=variants,
                                 pheno=pheno,
                                 args=args,
                                 logfile=logfile)

        elif method == "AIM":
            run_association_AIM(trees=trees,
                                inds=inds,
                                variants=variants,
                                pheno=pheno,
                                args=args,
                                ass_method=m,
                                window_size=args.ass_window_size,
                                trees_interval=trees_object.actual_trees_interval,
                                logfile=logfile)

        else:
            raise ValueError("Unknown association method '" + m + "'. Must be 'AIM' or 'GWAS'.")

        logfile.sub()

    logfile.info("- Done running association tests")

    # --------------------------------
    # clean up
    # --------------------------------
    if not args.no_clean_up:
        remove_files_with_pattern(args.out + '*HE-CP_result.txt')
        remove_files_with_pattern(args.out + '*HE-SD_result.txt')
        remove_files_with_pattern(args.out + '*.HEreg')
        remove_files_with_pattern(args.out + '*tmp.out')
        # remove_files_with_pattern(args.out + '*phenotypes.phen')
        remove_files_with_pattern(args.out + '*REML_result.txt')
        remove_files_with_pattern(args.out + '*REML.hsq')
        remove_files_with_pattern(args.out + '*REML.log')
        remove_files_with_pattern(args.out + '*.grm.id')
        remove_files_with_pattern(args.out + '*.grm.bin')
        remove_files_with_pattern(args.out + '*.grm.N.bin')

    logfile.sub()


def OLS(genotypes, phenotypes):
    # add intercept
    # print("!!!!!len genotypes", len(genotypes))
    # print("!!!!!len phenotypes", len(phenotypes))
    # print("!!!!!len phenotypes nan", len(np.isnan(phenotypes)))
    # print("!!!!!len genotypes nan", len(np.isnan(genotypes)))
    genotypes = genotypes[~np.isnan(genotypes)]
    # raise ValueError("here")
    #    tmp = np.isnan(phenotypes)
    # np.delete(genotypes, tmp)

    genotypes_test = sm.tools.add_constant(genotypes)
    # phenotypes = np.delete(phenotypes, np.isnan(phenotypes))
    PVALUE = sm.OLS(phenotypes, genotypes_test, missing='drop').fit().pvalues[1]
    return PVALUE


def run_association_GWAS(trees, samp_ids, inds, variants: tvar, pheno: pt, args: TParameters, logfile):
    outname = args.out + "_GWAS"

    if args.imputation_ref_panel_tree_file is not None:
        logfile.info("- Using genotypes imputed with impute2 for GWAS:")
        logfile.add()

        if args.do_imputation:
            # impute
            imputation_obj = impute.TImpute()
            name_imputation_output = imputation_obj.run_impute(trees_sample=trees,
                                                               variants_sample=variants,
                                                               inds=inds,
                                                               imputation_ref_panel_tree_file=args.imputation_ref_panel_tree_file,
                                                               ploidy_ref=args.ploidy_ref,
                                                               genetic_map_file=args.genetic_map_file,
                                                               out=outname, logfile=logfile)
        else:
            if args.imputed_gen_file is None:
                raise ValueError(
                    "When --do_imputation is set to False, the imputed genotypes must be provided "
                    "with --imputed_gen_file parameter")

            logfile.info(
                "- Assuming imputation was already run, reading imputed genotypes from " + args.imputed_gen_file)
            name_imputation_output = args.imputed_gen_file

        # read imputed genotypes
        gt_matrix_imputed, pos = impute.TImpute.read_imputed_gt(name_imputation_output=name_imputation_output,
                                                                variants_sample=variants,
                                                                trees_interval=args.trees_interval,
                                                                logfile=logfile)
        logfile.sub()

        # run association tests
        GWAS = at.TAssociationTestingGWAS(phenotypes=pheno, num_typed_variants=gt_matrix_imputed.shape[1])
        GWAS.test_with_positions_from_X_matrix(X=gt_matrix_imputed, positions=pos,
                                               variants_sample=variants,
                                               phenotypes=pheno,
                                               logfile=logfile)
        GWAS.write_to_file_with_X_matrix(positions=pos, name=outname, logfile=logfile)

    else:
        logfile.info("- Using genotypes from tree file for GWAS")
        # run association tests
        GWAS = at.TAssociationTestingGWAS(phenotypes=pheno, num_typed_variants=variants.num_typed)
        GWAS.test_with_variants_object(trees=trees, samp_ids=samp_ids, variants=variants, phenotypes=pheno, inds=inds, logfile=logfile)
        GWAS.write_to_file(variants, outname, logfile)
        # GWAS.manhattan_plot(variant_positions=variants.info['position'], plots_dir=plots_dir)


def get_AIM_test_object(test_name: str, phenotypes, pheno_file, num_associations, outname: str, logfile,
                        args: TParameters):
    """
    :param test_name: str describing the association algorithm to be used
    :param phenotypes:
    :param pheno_file:
    :param num_associations:
    :param outname:
    :param logfile:
    :param args:
    :return: TAssociationTestingRegions object
    """
    # rename deprecated test names
    if test_name == "HE":
        logfile.info("WARNING: AIM method 'HE' is deprecated. Use 'GCTA_HE'")
        test_name = "GCTA_HE"
    if test_name == "REML":
        logfile.info("WARNING: AIM method 'REML' is deprecated. Use 'GCTA_REML'")
        test_name = "GCTA_REML"

    # create object
    if test_name == "GCTA_HE":
        test_obj = at.TAssociationTestingRegionsGCTA_HE(phenotypes, num_associations, test_name=test_name,
                                                        pheno_file=pheno_file, outname=outname, args=args,
                                                        logfile=logfile)
    elif test_name == "GCTA_REML":
        test_obj = at.TAssociationTestingRegionsGCTA_REML(phenotypes, num_associations, test_name=test_name,
                                                          pheno_file=pheno_file, outname=outname, args=args,
                                                          logfile=logfile)
    elif test_name == "glimix_REML":
        test_obj = at.TAssociationTestingRegionsGlimix(phenotypes, num_associations)
    elif test_name == "mtg2":
        raise ValueError("mtg2 is not implemented yet!")
        # test_obj = at.TAssociationTestingRegionsMtg2(phenotypes, num_associations)
    else:
        raise ValueError("Did not recognize " + str(test_name) + " as a association test type")

    return test_obj


def get_window_starts_and_ends(window_size: int, trees_interval: list):
    """
    Get coordinates of windows of non-overlapping windows of size window_size
    """

    # get window ends
    window_ends = []
    if window_size >= (trees_interval[1] - trees_interval[0]):
        window_ends.append(trees_interval[1])
    else:
        num_tests = (trees_interval[1] - trees_interval[0]) / window_size
        for w in range(int(np.floor(num_tests))):
            window_ends.append(trees_interval[0] + (w + 1) * window_size)
        # add last bit (smaller window)
        if window_ends[-1] < trees_interval[1]:
            window_ends.append(trees_interval[1])

    # get window starts
    window_starts = [x - window_size for x in window_ends[:-1]]
    window_starts.append(window_starts[-1] + window_size)  # the last start should not be last end - window size

    return window_starts, window_ends


def get_proportion_of_tree_within_window(window_start: int, window_end: int, tree_start: float, tree_end: float):
    """
    Return proportion of tree that is within a genomic window. I guess the window_end and tree_end are not included.
    @param window_start:
    @param window_end:
    @param tree_start:
    @param tree_end:
    @return float: proportion
    """
    tree_length = tree_end - tree_start
    if (tree_start >= window_start and tree_end < window_end) \
            or (tree_start <= window_start and tree_end >= window_end):
        # tree is completely in window
        return 1.0
    elif tree_start < window_end <= tree_end:
        proportion = (window_end - tree_start) / tree_length
        # tree overlaps with window end
        if 0 > proportion or proportion > 1:
            raise ValueError("Proportion " + str(proportion) + " out of bounds")
        return proportion
    elif window_start < tree_end < window_end:
        # tree overlaps with window start
        proportion = (tree_end - window_start) / tree_length
        if 0 > proportion or proportion > 1:
            raise ValueError("Proportion " + str(proportion) + " out of bounds")
        return (tree_end - window_start) / tree_length
    else:
        # there is no overlap
        return 0.0


def write_matrices_for_testing(cholesky_global_GRM_for_cor: cov, covariance_obj: cov, inds: tind, outname: str,
                               covariances_picklefile: IO, index: int, logfile: IndentedLoggerAdapter):

    if covariance_obj.write(out=outname, inds=inds, covariances_picklefile=covariances_picklefile):
        if cholesky_global_GRM_for_cor:
            # calculate cholesky of local GRM --> get A
            covariance_obj.calculate_cholesky_decomposition(ploidy=inds.ploidy, logfile=logfile)
            A = covariance_obj.get_cholesky_decomposition(ploidy=inds.ploidy)

            # add A to global_GRM
            Q = A + cholesky_global_GRM_for_cor.get_cholesky_decomposition(ploidy=inds.ploidy)

            # write sum (Q) to GCTA
            cholesky_global_GRM_for_cor.write_gcta_format(covariance_matrix=Q,
                                                          mu=cholesky_global_GRM_for_cor.mu,
                                                          inds=inds, out=outname + "_cov")
        return True

    else:
        print("did not run association because covariance objects could not be written at index", index)


def test_window_for_association(covariance_obj: cov, inds: tind, AIM_methods: list, outname: str, window_index: int,
                                phenotypes_obj: pt, cholesky_global_GRM_for_cor: cov, covariances_picklefile: IO,
                                logfile: IndentedLoggerAdapter):
    """
    Finalize covariance object and run association test. If global GRM is given and the testing method is GCTA,
    first also calculate the correlation between the local and global matrix

    :param covariance_obj:
    :param inds:
    :param AIM_methods:
    :param outname:
    :param window_index:
    :param phenotypes_obj:
    :param global_GRM_for_cor: TCovariance, if provided the correlation between the local and global matrix are written to file
    :param covariances_picklefile: bool, should local covariance matrix be written to a pickle file
    :return:
    """

    covariance_obj.finalize(inds=inds)

    for m in AIM_methods:
        if m.name == "regions_glimix":
            m.test(index=window_index,
                   out=outname,
                   inds=inds,
                   phenotypes_object=phenotypes_obj,
                   covariance_object=covariance_obj,
                   covar=None,
                   covariances_picklefile=None,
                   )
        elif m.name == "regions_mtg2":
            m.test(index=window_index,
                   out=outname,
                   inds=inds,
                   phenotypes_object=phenotypes_obj,
                   covariance_object=covariance_obj,
                   covar=None,
                   covariances_picklefile=None,
                   )
        elif m.name == "regions_GCTA_HE" or m.name == "regions_GCTA_REML":

            if write_matrices_for_testing(cholesky_global_GRM_for_cor=cholesky_global_GRM_for_cor,
                                          covariance_obj=covariance_obj,
                                          inds=inds,
                                          outname=outname,
                                          index=window_index,
                                          covariances_picklefile=covariances_picklefile,
                                          logfile=logfile):

                m.test(index=window_index,
                       out=outname,
                       inds=inds,
                       phenotypes_object=None,
                       covariance_object=covariance_obj,
                       covar=None,
                       covariances_picklefile=covariances_picklefile)

        else:
            raise ValueError("There is no association test implemented for '" + m.name + "'")

    covariance_obj.clear()


def loop_windows_variant_based_covariance_testing(covariance_obj: cov, AIM_methods: list, variants: tvar,
                                                  window_ends: list, window_starts: list, num_tests: int, inds: tind,
                                                  covariances_picklefile: IO, pheno: pt,
                                                  cholesky_global_GRM_for_cor: cov, logfile, outname: str):
    """
    Write covariance calculated based on variants within a window (can be one tree) to file and test it for association with
    phenotypes. Currently, the only covariance type based on variants is GRM.

    :param covariance_obj:
    :param AIM_methods:
    :param variants:
    :param window_ends:
    :param window_starts:
    :param num_tests:
    :param inds:
    :param covariances_picklefile:
    :param pheno:
    :param cholesky_global_GRM_for_cor:
    :param logfile:
    :param outname:
    :return:
    """
    window_ends_copy = window_ends.copy()
    window_starts_copy = window_starts.copy()

    # log progress
    start = time.time()

    for w in range(num_tests):  #
        covariance_obj.calculate_GRM(window_beginning=window_starts[w], window_end=window_ends[w],
                                     variants=variants, inds=inds)
        tmpCov = covariance_obj.get_covariance_matrix(inds.ploidy)

        if tmpCov is not None:
            for m in AIM_methods:
                # covariance_obj.write(out=outname, inds=inds, covariances_picklefile=covariances_picklefile) removed
                # because already in run association function
                m.test(index=w,
                       out=outname,
                       covariance_object=covariance_obj,
                       phenotypes_object=pheno,
                       inds=inds,
                       covar=None,
                       covariances_picklefile=covariances_picklefile)
            covariance_obj.clear()

        # log progress
        if w % 1000 == 0:
            end = time.time()
            logfile.info("- Ran AIM for " + str(w) + " windows in " + str(round(end - start)) + " s")

    for m in AIM_methods:
        m.write_association_results_to_file(window_starts=window_starts_copy,
                                            window_ends=window_ends_copy,
                                            phenotypes=pheno,
                                            logfile=logfile,
                                            out=outname)


def loop_windows_tree_based_covariance_testing(trees, covariance_obj: cov, AIM_methods: list,
                                               window_ends: list,
                                               window_starts: list, window_size: int, skip_first_tree: bool,
                                               inds: tind, pheno: pt, covariances_picklefile: IO,
                                               cholesky_global_GRM_for_cor: cov, logfile, outname: str,
                                               limit_association_tests: int):
    """
    Loop over windows, write necessary files and run association tests
    :param trees:
    :param covariance_obj:
    :param AIM_methods:
    :param window_ends:
    :param window_starts:
    :param window_size: if not specified, trees will be used
    :param skip_first_tree:
    :param inds: TInds
    :param pheno:
    :param covariances_picklefile: when true picklefiles with all covariances are written to a pickle file
    :param cholesky_global_GRM_for_cor:
    :param logfile:
    :param outname:
    :param limit_association_tests:
    :return:
    """

    window_ends_copy = window_ends.copy()
    window_starts_copy = window_starts.copy()
    window_index = 0

    # log progress
    start = time.time()

    # windows are given by trees
    if window_size is None:
        tree = tskit.Tree(trees)
        while tree.next() and window_index < limit_association_tests:
            tree_obj = tt.TTree(tree)

            # print("----------------")
            # print("tree_obj.start", tree_obj.start, ", tree_obj.end", tree_obj.end, ", window_index", window_index)

            if tree_obj.is_testable(skip_first_tree):
                if window_size is None:
                    covariance_obj.add_tree(tree_obj=tree_obj, inds=inds, proportion=1.0)
                    test_window_for_association(covariance_obj=covariance_obj,
                                                inds=inds,
                                                AIM_methods=AIM_methods,
                                                phenotypes_obj=pheno,
                                                outname=outname,
                                                window_index=window_index,
                                                cholesky_global_GRM_for_cor=cholesky_global_GRM_for_cor,
                                                covariances_picklefile=covariances_picklefile,
                                                logfile=logfile)

            window_index += 1
            # log progress
            if tree.index % 100 == 0:
                end = time.time()
                logfile.info("- Ran AIM for " + str(tree.index) + " trees in " + str(round(end - start)) + " s")

    # there is a window size
    else:
        tree = tskit.Tree(trees)
        while tree.next() and window_index < limit_association_tests:
            tree_obj = tt.TTree(tree)

            # print("----------------")
            # print("tree_obj.start", tree_obj.start, ", tree_obj.end", tree_obj.end, ", window_index", window_index,
            #       ", window_starts[0]", window_starts[0], ", window_ends[0]", window_ends[0])

            if tree_obj.is_testable(skip_first_tree):
                if tree_obj.start >= window_ends[0]:
                    logfile.info("WARNING: Tree start (" + str(tree_obj.start) + ") is >= window end ("
                                 + str(window_ends[0]) + "). This should only happen at beginning of ARG if "
                                                         "trees_interval parameter was not set to coincide with part of ARG tested. "
                                                         "Increasing window index.")
                    while tree_obj.start >= window_ends[0]:
                        window_ends.pop(0)
                        window_starts.pop(0)
                        window_index += 1

                # tree start is definitely within window based on previous while loop
                proportion = get_proportion_of_tree_within_window(window_start=window_starts[0],
                                                                  window_end=window_ends[0],
                                                                  tree_start=tree_obj.start,
                                                                  tree_end=tree_obj.end)
                # print("tree is in window with index", window_index, "with proportion", proportion, ". window has coordinates", window_starts[0], window_ends[0])

                if 0.0 < proportion <= 1.0:
                    # part of tree is in this window, part in next --> needs to be added to both windows, and first
                    # window needs to be tested
                    covariance_obj.add_tree(tree_obj=tree_obj, inds=inds, proportion=proportion)
                    # print("added first part of tree with proportion", proportion, "to window with index",
                    #       window_index)

                    # add rest of tree to next window. Tree might span multiple windows, so need while loop
                    while (proportion < 1.0 and tree_obj.end >= window_ends[0]) \
                            or (proportion == 1.0 and tree_obj.end >= window_ends[
                        0]) and window_index < limit_association_tests:
                        test_window_for_association(covariance_obj=covariance_obj,
                                                    inds=inds,
                                                    AIM_methods=AIM_methods,
                                                    outname=outname,
                                                    window_index=window_index,
                                                    phenotypes_obj=pheno,
                                                    cholesky_global_GRM_for_cor=cholesky_global_GRM_for_cor,
                                                    covariances_picklefile=covariances_picklefile,
                                                    logfile=logfile)

                        if len(window_ends) == 1:  # that was the last window
                            break

                        # move to next window
                        window_ends.pop(0)
                        window_starts.pop(0)
                        window_index += 1

                        proportion = get_proportion_of_tree_within_window(window_start=window_starts[0],
                                                                          window_end=window_ends[0],
                                                                          tree_start=tree_obj.start,
                                                                          tree_end=tree_obj.end)
                        covariance_obj.add_tree(tree_obj=tree_obj, inds=inds, proportion=proportion)

                        # print("added next part of tree to window with index", window_index, "coordinates",
                        #       window_starts[0], window_ends[0], "and proportion", proportion)

                # log progress
                if tree.index % 100 == 0:
                    end = time.time()
                    logfile.info("- Ran AIM for " + str(tree.index) + " trees in " + str(round(end - start)) + " s")
            # else:
            #     print("tree was not usable")

    # write association test results to file
    for m in AIM_methods:
        m.write_association_results_to_file(window_starts=window_starts_copy,
                                            window_ends=window_ends_copy,
                                            out=outname,
                                            phenotypes=pheno,
                                            logfile=logfile)


def run_association_AIM(trees, inds, variants, pheno, args, ass_method, window_size, trees_interval,
                        logfile):
    # ----------------
    # initialize
    # ----------------

    if args.AIM_method is None:
        raise ValueError("ERROR: No method for tree association provided. Use '--AIM_method' to set method.")

    logfile.info("- Reading tree estimations for tree-based association from " + args.tree_file)

    # determine covariance type
    covariance = ass_method.split(':')[1]
    if covariance not in ["scaled", "eGRM", "GRM"]:
        raise ValueError("Unknown covariance method '" + covariance + "'. Must be one of 'scaled', 'eGRM', 'GRM'.")
    logfile.info("- Writing output files with suffix '_" + covariance + "'")
    outname = args.out + "_" + covariance
    logfile.info("- Running associations tests using covariance type " + covariance + " for a sequence of trees")
    logfile.add()

    # define number and coordinates of windows
    window_starts = trees.breakpoints(as_array=True)[0:trees.num_trees]
    window_ends = trees.breakpoints(as_array=True)[1:]
    num_tests = trees.num_trees

    if window_size is not None:
        window_starts, window_ends = get_window_starts_and_ends(window_size=window_size, trees_interval=trees_interval)
        num_tests = len(window_ends)

    # initialize and write phenotypes
    if covariance == "eGRM" or covariance == "GRM":
        pheno.set_missing_phenotype_status(inds)
        pheno.write_to_file_gcta_eGRM(inds=inds, out=args.out, logfile=logfile)
    else:
        pheno.write_to_file_gcta_scaled(out=outname, logfile=logfile)
    if args.simulate_phenotypes is True:
        pheno.find_causal_trees(trees)
        pheno.find_causal_windows(window_ends=window_ends, window_starts=window_starts)

    pheno_file = args.out + "_phenotypes.phen"
    # if args.simulate_phenotypes:
    #     pheno_file = args.out + "_phenotypes.phen"

    # create association method objects
    logfile.info("- Using association test methods " + str(args.AIM_method) + " for a sequence of trees")
    AIM_methods = []
    for m in args.AIM_method:
        test_obj = get_AIM_test_object(test_name=m,
                                       phenotypes=pheno,
                                       pheno_file=pheno_file,
                                       num_associations=num_tests,
                                       outname=outname,
                                       logfile=logfile,
                                       args=args)
        AIM_methods.append(test_obj)

    # create covariance type object
    covariance_obj = cov.get_covariance_object(covariance)

    # ----------------
    # run association tests
    # ----------------
    logfile.info("- Running association tests")

    # write covariances to picklefile for later comparison?
    covariance_picklefiles = None
    if args.covariance_picklefiles:
        covariances_picklefiles = open(outname + "_matrices_" + covariance_obj.covariance_type + ".pickle", "wb")

    cholesky_global_GRM_for_cor = None
    if args.coreGREML_model:
        global_GRM = cov.get_covariance_object("base")
        global_GRM.read_gcta_format(args.population_structure_matrix, inds.ploidy)
        global_GRM.calculate_cholesky_decomposition(ploidy=inds.ploidy, logfile=logfile)
        cholesky_global_GRM_for_cor = global_GRM
        cholesky_global_GRM_for_cor.forget_original_matrix()

    # variant based covariance
    if covariance == "GRM":
        loop_windows_variant_based_covariance_testing(covariance_obj=covariance_obj,
                                                      AIM_methods=AIM_methods,
                                                      variants=variants,
                                                      window_ends=window_ends,
                                                      window_starts=window_starts,
                                                      num_tests=num_tests,
                                                      inds=inds,
                                                      covariances_picklefile=covariance_picklefiles,
                                                      cholesky_global_GRM_for_cor=cholesky_global_GRM_for_cor,
                                                      pheno=pheno,
                                                      logfile=logfile,
                                                      outname=outname)

    # window based covariance (need to loop over trees)
    else:
        loop_windows_tree_based_covariance_testing(trees=trees,
                                                   covariance_obj=covariance_obj,
                                                   AIM_methods=AIM_methods,
                                                   window_ends=window_ends,
                                                   window_starts=window_starts,
                                                   window_size=window_size,
                                                   inds=inds,
                                                   skip_first_tree=args.skip_first_tree,
                                                   covariances_picklefile=covariance_picklefiles,
                                                   cholesky_global_GRM_for_cor=cholesky_global_GRM_for_cor,
                                                   pheno=pheno,
                                                   logfile=logfile,
                                                   outname=outname,
                                                   limit_association_tests=args.limit_association_tests)

    # close covariance picklefiles
    if args.covariance_picklefiles:
        covariance_picklefiles.close()

    logfile.sub()


def remove_files_with_pattern(pattern):
    """

    @param pattern: str, can contain '*'
    @return:
    """
    fileList = glob.glob(pattern)
    for file in fileList:
        if os.path.exists(file):
            os.remove(file)


def run_covariance_correlation(args, logfile):
    """

    @param args: TParameters
    @param logfile:
    @return:
    """
    logfile.info("- TASK: covarianceCorrelations")
    logfile.info("- Calculating correlations for all windows from the following covariance picklefiles " +
                 str(args.covariance_picklefiles))

    if args.covariance_picklefiles is None or len(args.covariance_picklefiles) != 2:
        raise ValueError("Must provide two files with covariance matrices for same windows with "
                         "'--covariance_picklefiles'")

    # open picklefiles
    f1 = open(args.covariance_picklefiles[0], 'rb')
    f2 = open(args.covariance_picklefiles[1], 'rb')

    # open output file
    output = open(args.out + "_correlations.csv", 'w')

    while True:
        try:
            L1 = pickle.load(f1)
            L2 = pickle.load(f2)

            diags = np.diag_indices(L1.shape[0])
            non_diags = np.where(~np.eye(L1.shape[0], dtype=bool))
            V1 = L1[non_diags].flatten()
            V2 = L2[non_diags].flatten()

            corr = np.corrcoef(V1, V2)
            output.write(str(corr[0][1]) + "\n")

        except EOFError:
            logfile.info("- Completed reading covariance picklefile")
            break

    f1.close()
    f2.close()
    output.close()
