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
import stat
import pandas as pd


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
    logfile.add()

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

    plots_dir = args.out + "_plots/"
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    sample_ids = trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(ploidy=args.ploidy, num_haplotypes=N, relate_sample_names_file=args.relate_sample_names,
                            logfile=logfile)

    # --------------------------------
    # create diploids and variants
    # --------------------------------

    # TODO: trees_orig and variants_orig should be initialized at the same time, e.g. together in one function. We
    #  should not have 2 tree files and 2 variant files just floating around separately. Maybe TTrees class could
    #  contain the tskit sequence and variants

    # TODO: find way to save variants in their tskit format without needing to read the original tree. I only need
    #  original tree in association task for this. It would be nice if the only tree that needs to be read would be
    #  estimated tree

    #  do not provide variant file here but have it estimated from tree, otherwise variants and tree
    #  won't match (tree only contains typed variants). The variant file is only useful for simulating phenotypes to
    #  be able to keep track of untyped variants (i.e. for variants_orig)
    logfile.info("- Compiling variant info from trees file")
    variants = tvar.TVariantsFiltered(ts_object=trees,
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
                               plots_dir=plots_dir,
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


def run_association_GWAS(trees, inds, variants, pheno, args, logfile):
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
        logfile.info("- Using genotypes from tree file for GWAS:")
        # run association tests
        GWAS = at.TAssociationTestingGWAS(phenotypes=pheno, num_typed_variants=variants.num_typed)
        GWAS.test_with_variants_object(variants=variants, phenotypes=pheno, inds=inds, logfile=logfile)
        GWAS.write_to_file(variants, outname, logfile)
        # GWAS.manhattan_plot(variant_positions=variants.info['position'], plots_dir=plots_dir)


def write_GCTA_command_script(test_name, pheno_file, outname, logfile, args):
    with open(outname + "_run_" + test_name + ".sh", 'w') as f:
        f.write("#!/bin/bash\n")
        if args.population_structure and args.population_structure_pca_num_eigenvectors is None:
            logfile.info("- Writing gcta command file to test a model containing the local GRM and a global GRM as "
                         "random effects")
            write_GCTA_command_file_mgrm(testing_method=test_name,
                                         outname=outname,
                                         pheno_file=pheno_file,
                                         outfile=f,
                                         GCTA=args.GCTA,
                                         num_GCTA_threads=args.num_gcta_threads,
                                         additional_gcta_params=args.additional_gcta_params)

        elif args.population_structure and args.population_structure_pca_num_eigenvectors \
                and args.do_all_stratification_correction:
            logfile.info("- Writing gcta command file to run a PCA on the population structure GRM, and then test a "
                         "model containing the local GRM and a global GRM as a random effects, and the PCs as "
                         "fixed effects")
            write_GCTA_command_file_mgrm_pca(testing_method=test_name,
                                             outname=outname,
                                             pheno_file=pheno_file,
                                             num_eigenvectors=args.population_structure_pca_num_eigenvectors,
                                             population_structure_matrix=args.population_structure,
                                             outfile=f,
                                             GCTA=args.GCTA,
                                             num_GCTA_threads=args.num_gcta_threads)

        elif args.population_structure and args.population_structure_pca_num_eigenvectors \
                and not args.do_all_stratification_correction:
            logfile.info("- Writing gcta command file to run a PCA on the population structure GRM, and then test a "
                         "model containing the local GRM as a random effect, and the PCs as fixed effects")
            write_GCTA_command_file_grm_pca(testing_method=test_name,
                                            outname=outname,
                                            pheno_file=pheno_file,
                                            num_eigenvectors=args.population_structure_pca_num_eigenvectors,
                                            population_structure_matrix=args.population_structure,
                                            outfile=f,
                                            GCTA=args.GCTA,
                                            num_GCTA_threads=args.num_gcta_threads)
        else:
            logfile.info("- Writing gcta command file to test a model containing the local GRM as a random effect")
            write_GCTA_command_file_grm(testing_method=test_name,
                                        outname=outname,
                                        pheno_file=pheno_file,
                                        outfile=f,
                                        GCTA=args.GCTA,
                                        num_GCTA_threads=args.num_gcta_threads)

    st = os.stat(outname + "_run_" + test_name + ".sh")
    os.chmod(outname + "_run_" + test_name + ".sh", st.st_mode | stat.S_IEXEC)


def get_AIM_test_object(test_name, phenotypes, pheno_file, num_associations, outname, logfile, args):
    # rename deprecated test names
    if test_name == "HE":
        logfile.info("WARNING: AIM method 'HE' is deprecated. Use 'GCTA_HE'")
        test_name = "GCTA_HE"
    if test_name == "REML":
        logfile.info("WARNING: AIM method 'REML' is deprecated. Use 'GCTA_REML'")
        test_name = "GCTA_REML"

    # create object
    if test_name == "GCTA_HE":
        write_GCTA_command_script(test_name=test_name, pheno_file=pheno_file, outname=outname, args=args, logfile=logfile)
        test_obj = at.TAssociationTestingRegionsGCTA_HE(phenotypes, num_associations)
    elif test_name == "GCTA_REML":
        write_GCTA_command_script(test_name=test_name, pheno_file=pheno_file, outname=outname, args=args, logfile=logfile)
        test_obj = at.TAssociationTestingRegionsGCTA_REML(phenotypes, num_associations)
    elif test_name == "glimix_REML":
        test_obj = at.TAssociationTestingRegionsGlimix(phenotypes, num_associations)
    else:
        raise ValueError("Did not recognize " + str(test_name) + " as a association test type")

    return test_obj


def get_window_starts_and_ends(window_size, trees_interval):
    """
    @param skip_first_tree: bool
    @param trees_object: TTrees
    @param trees_interval: int genomic region covered by ARG
    @param window_size: int
    @return list: list of window ends

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


def get_proportion_of_tree_within_window(window_start, window_end, tree_start, tree_end):
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


def test_window_for_association(covariance_obj, inds, AIM_methods, outname, window_index, phenotypes_obj,
                                covariances_picklefile):
    covariance_obj.finalize(inds=inds)

    for m in AIM_methods:
        if m.name == "regions_glimix":
            m.run_association(index=window_index,
                              out=outname,
                              inds=inds,
                              phenotypes_object=phenotypes_obj,
                              covariance_object=covariance_obj,
                              covar=None,
                              covariances_picklefile=None,
                              )
        else:
            m.run_association(index=window_index,
                              out=outname,
                              inds=inds,
                              phenotypes_object=None,
                              covariance_object=covariance_obj,
                              covar=None,
                              covariances_picklefile=covariances_picklefile)

    covariance_obj.clear()


def run_variant_based_covariance_testing(covariance_obj, AIM_methods, variants, window_ends, window_starts, num_tests,
                                         inds, covariances_picklefile, pheno, logfile, outname, population_structure):
    """
    Write covariance calculated based on variants within a window (can be one tree) to file and test it for association with
    phenotypes. Currently, the only covariance type based on variants is GRM.
    @param covariances_picklefile: write covariance matrix to a pickle file
    @param covariance_obj: TCovariance
    @param AIM_methods: list
    @param variants: TVariants
    @param window_ends: list
    @param window_starts: list
    @param num_tests: int
    @param inds: TInds
    @param logfile: Tlogger
    @param outname: str
    @return: None
    :param pheno: TPhenotype
    """
    window_ends_copy = window_ends.copy()
    window_starts_copy = window_starts.copy()

    # log progress
    start = time.time()

    for w in range(num_tests):  #
        covariance_obj.calculate_GRM(window_beginning=window_starts[w], window_end=window_ends[w],
                                     variants=variants, inds=inds)
        tmpCov = covariance_obj.covariance_matrix

        if tmpCov is not None:
            for m in AIM_methods:
                # covariance_obj.write(out=outname, inds=inds, covariances_picklefile=covariances_picklefile) removed
                # because already in run association function
                m.run_association(index=w,
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


def run_tree_based_covariance_testing(trees, covariance_obj, AIM_methods, window_ends, window_starts,
                                      window_size, skip_first_tree, inds, pheno, covariances_picklefile,
                                      logfile, outname, limit_association_tests):
    """

    @param population_structure: str prefix of covariance matrix files used for correcting for population structure (can be None)
    @param covariances_picklefile:
    @param trees:
    @param covariance_obj:
    @param AIM_methods:
    @param window_ends:
    @param window_starts:
    @param window_size: if not specified, trees will be used
    @param skip_first_tree: bool
    @param inds: TInds
    @param write_covariance_picklefiles: should picklefiles with all covariances be written to a pickle file
    @param logfile:
    @param outname:
    @return:
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
                                                covariances_picklefile=covariances_picklefile)

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
                                                    covariances_picklefile=covariances_picklefile)

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


def write_GCTA_command_file_mgrm(testing_method, outname, pheno_file, outfile, GCTA, num_GCTA_threads, additional_gcta_params):
    """
    Write executable bash script for running association test with multiple random effects using GCTA

    @param testing_method:
    @param outname:
    @param pheno_file:
    @param outfile:
    @param GCTA:
    @param num_GCTA_threads:
    @return:
    """

    if testing_method == "GCTA_REML":
        gcta_string = GCTA + " --reml --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out " \
                      + outname + "_REML --reml-lrt 1 --threads " + str(num_GCTA_threads) + " --reml-maxit 500 "
        for p in additional_gcta_params:
            gcta_string += " --" + p
        outfile.write(gcta_string + " > " + outname + "_tmp.out\n")
    elif testing_method == "GCTA_HE":
        gcta_string = GCTA + " --HEreg --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out "\
                      + outname + "_HE --reml-lrt 1 --threads " + str(num_GCTA_threads) + " --reml-maxit 500 "
        for p in additional_gcta_params:
            gcta_string += " --" + p
        outfile.write(gcta_string + " > " + outname + "_tmp.out\n")

        # grep results
        outfile.write("sed -n '2,6p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\t\' > "
                      + outname + "_HE-CP_result.txt\n")
        outfile.write("sed -n '9,13p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\t\' > "
                      + outname + "_HE-SD_result.txt\n")


def write_GCTA_command_file_mgrm_pca(testing_method, outname, pheno_file, outfile, num_eigenvectors,
                                     population_structure_matrix, GCTA, num_GCTA_threads):
    """
    Write executable bash script for running association test with multiple random effects and fixed effects using GCTA

    @param testing_method:
    @param outname:
    @param pheno_file:
    @param outfile:
    @param GCTA:
    @param num_GCTA_threads:
    @return:
    """

    outfile.write(GCTA + " --grm " + population_structure_matrix + " --pca " + str(num_eigenvectors) + " --out "
                  + outname + "> " + outname + "_tmp2.out\n\n")

    if testing_method == "GCTA_REML":
        outfile.write(GCTA + " --reml --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out "
                      + outname + "_REML --reml-lrt 1 " + " --qcovar " + outname + ".eigenvec --threads " + str(
            num_GCTA_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
    elif testing_method == "GCTA_HE":
        outfile.write(
            GCTA + " --HEreg --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out "
            + outname + "_HE --reml-lrt 1 " + " --qcovar " + outname + ".eigenvec --threads " + str(
                num_GCTA_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
        # grep results
        outfile.write("sed -n '2,6p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\t\' > "
                      + outname + "_HE-CP_result.txt\n")
        outfile.write("sed -n '9,13p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\t\' > "
                      + outname + "_HE-SD_result.txt\n")


def write_GCTA_command_file_grm(testing_method, outname, pheno_file, outfile, GCTA, num_GCTA_threads):
    """
    Write executable bash script for running association test with only the local eGRM as random effects using GCTA

    @param testing_method:
    @param outname:
    @param pheno_file:
    @param outfile:
    @param GCTA:
    @param num_GCTA_threads:
    @return:
    """
    if testing_method == "GCTA_REML":
        outfile.write(GCTA + " --reml --grm " + outname + " --pheno " + pheno_file + " --out " + outname
                      + "_REML --threads " + str(
            num_GCTA_threads) + " --reml-maxit 500  > " + outname + "_tmp.out\n")
    elif testing_method == "GCTA_HE":
        outfile.write(GCTA + " --HEreg --grm " + outname + " --pheno " + pheno_file + " --out " + outname
                      + "_HE --threads " + str(
            num_GCTA_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
        # grep results
        outfile.write("sed -n '2,4p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                      + outname + "_HE-CP_result.txt\n")
        outfile.write("sed -n '7,9p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                      + outname + "_HE-SD_result.txt\n")
    else:
        raise ValueError("Unrecognized AIM testing method")


def write_GCTA_command_file_grm_pca(testing_method, outname, pheno_file, outfile, num_eigenvectors,
                                    population_structure_matrix, GCTA, num_GCTA_threads):
    """
    Write executable bash script for running association test with local eGRM as random effects and PCA of global
    population structure matrix using GCTA

    @param num_eigenvectors:
    @param population_structure_matrix:
    @param testing_method:
    @param outname:
    @param pheno_file:
    @param outfile:
    @param GCTA:
    @param num_GCTA_threads:
    @return:
    """
    outfile.write(GCTA + " --grm " + population_structure_matrix + " --pca " + str(num_eigenvectors) + " --out "
                  + outname + "> " + outname + "_tmp2.out\n\n")

    if testing_method == "GCTA_REML":
        outfile.write(
            GCTA + " --reml --grm " + outname + " --pheno " + pheno_file + " --out " + outname + "_REML" + " --qcovar " + outname + ".eigenvec --threads "
            + str(num_GCTA_threads) + " --reml-maxit 500  > " + outname + "_tmp.out\n")
    elif testing_method == "GCTA_HE":
        outfile.write(
            GCTA + " --HEreg --grm " + outname + " --pheno " + pheno_file + " --out " + outname + "_HE --qcovar " + outname + ".eigenvec "
            + " --threads " + str(num_GCTA_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
        # grep results
        outfile.write("sed -n '2,4p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                      + outname + "_HE-CP_result.txt\n")
        outfile.write("sed -n '7,9p' " + outname + "_" + testing_method + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                      + outname + "_HE-SD_result.txt\n")
    else:
        raise ValueError("Unknown testing method '" + str(testing_method) + "'")


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

    # write GCTA files and scripts
    if args.population_structure:
        with open(outname + '_multi_grm.txt', 'w') as f:
            f.write(outname + '\n')
            f.write(args.population_structure + '\n')

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
    covariances_picklefile = None
    if args.write_covariance_picklefiles:
        covariances_picklefile = open(outname + "_matrices_" + covariance_obj.covariance_type + ".pickle", "wb")

    # variant based covariance
    if covariance == "GRM":
        run_variant_based_covariance_testing(covariance_obj=covariance_obj,
                                             AIM_methods=AIM_methods,
                                             variants=variants,
                                             window_ends=window_ends,
                                             window_starts=window_starts,
                                             num_tests=num_tests,
                                             inds=inds,
                                             covariances_picklefile=covariances_picklefile,
                                             pheno=pheno,
                                             logfile=logfile,
                                             outname=outname,
                                             population_structure=args.population_structure)

    # window based covariance (need to loop over trees)
    else:
        run_tree_based_covariance_testing(trees=trees,
                                          covariance_obj=covariance_obj,
                                          AIM_methods=AIM_methods,
                                          window_ends=window_ends,
                                          window_starts=window_starts,
                                          window_size=window_size,
                                          inds=inds,
                                          skip_first_tree=args.skip_first_tree,
                                          covariances_picklefile=covariances_picklefile,
                                          pheno=pheno,
                                          logfile=logfile,
                                          outname=outname,
                                          limit_association_tests=args.limit_association_tests)

    # close covariance picklefiles
    if args.write_covariance_picklefiles:
        covariances_picklefile.close()

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
