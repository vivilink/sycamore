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

    trees, args.trees_interval = tt.read_trees(tree_file=args.tree_file,
                                               trees_interval=args.trees_interval,
                                               trees_interval_start=args.trees_interval_start,
                                               trees_interval_end=args.trees_interval_end,
                                               logfile=logfile)

    # interval = np.array([173584350,175000000])
    # trees_tmp = trees.keep_intervals([interval], simplify=True)
    # trees_tmp.dump("/data/ARGWAS/hawaiians/chr5.part-05_1mb.trees")

    plots_dir = args.out + "_plots/"
    if not os.path.exists(plots_dir):
        os.mkdir(plots_dir)

    # --------------------------------
    # create diploids and variants
    # --------------------------------
    sample_ids = trees.samples()
    N = len(sample_ids)
    inds = tind.Individuals(ploidy=args.ploidy, num_haplotypes=N, relate_sample_names_file=args.relate_sample_names, logfile=logfile)
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

    if args.simulate_phenotypes:
        if args.tree_file_simulated is None:
            raise ValueError("To simulate phenotypes based on untyped variants the simulated trees need to be "
                             "provided with 'tree_file_simulated'.")

        logfile.info("- Reading simulated tree used for simulating phenotypes from " + args.tree_file_simulated)
        trees_orig = tskit.load(args.tree_file_simulated)

        if trees_orig.num_samples != trees.num_samples:
            raise ValueError("The trees provided with params 'tree_file_simulated' and 'tree_file' must have the same "
                             "number of samples")

        # variants_orig are used to simulate phenotypes. They need to be consistent with original tree and the typed
        # status that might have been defined earlier with a variants file. The causal mutation should not be affected
        # by a freq filter
        variants_orig = tvar.TVariantsFiltered(ts_object=trees_orig,
                                               samp_ids=sample_ids,
                                               min_allele_freq=0,
                                               max_allele_freq=1,
                                               prop_typed_variants=1,
                                               pos_float=args.pos_float,
                                               random=random,
                                               logfile=logfile,
                                               filtered_variants_file=args.variants_file)

        logfile.info("- Phenotypes:")
        logfile.add()

        pheno = pt.PhenotypesSimulated(variants=variants_orig, num_inds=inds.num_inds)

        pheno.simulate(args=args, r=random, logfile=logfile, variants_orig=variants_orig, trees=trees, inds=inds,
                       plots_dir=plots_dir)
        logfile.sub()

    else:
        if args.pheno_file_BMI:
            pheno = pt.PhenotypesBMI()
            pheno.initialize_from_file(filename=args.pheno_file_BMI, inds=inds, out=args.out, logfile=logfile)
            # TODO: maybe restrict tree to inds for which we have phenotypes here
        else:
            pheno = pt.Phenotypes()
            pheno.initialize_from_file(filename=args.pheno_file, inds=inds, out=args.out, logfile=logfile)

    # --------------------------------
    # run association tests and plot
    # --------------------------------
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


def get_AIM_test_object(test_name, phenotypes, num_associations):
    if test_name == "HE":
        test_obj = at.TAssociationTestingRegionsGCTA_HE(phenotypes, num_associations)
    elif test_name == "REML":
        test_obj = at.TAssociationTestingRegionsGCTA_REML(phenotypes, num_associations)
    else:
        raise ValueError("Did not recognize " + str(test_name) + " as a association test type")

    return test_obj


def get_window_ends(window_size, trees_interval):
    """
    @param window_size int: window size
    @param trees_interval int: genomic region covered by ARG
    @return list: list of window ends
    """
    window_ends = []
    num_tests = (trees_interval[1] - trees_interval[0]) / window_size

    for w in range(int(np.floor(num_tests))):
        window_ends.append(trees_interval[0] + (w + 1) * window_size)
    # add last bit (smaller window)
    if window_ends[-1] < trees_interval[1]:
        window_ends.append(trees_interval[1])

    return window_ends


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


def write_and_test_window_for_association(covariance_obj, inds, AIM_methods, outname, window_index,
                                          covariances_picklefile):
    if covariance_obj.covariance_matrix_haploid is None:
        raise ValueError("trying to test empty covariance matrix for association at window index " + str(window_index))
    covariance_obj.finalize(inds=inds)
    covariance_obj.write(out=outname, inds=inds, covariances_picklefile=covariances_picklefile)
    for m in AIM_methods:
        m.run_association(index=window_index, out=outname)
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
    """
    window_ends_copy = window_ends.copy()
    window_starts_copy = window_starts.copy()

    # log progress
    start = time.time()

    for w in range(num_tests):  #
        tmpCov, tmpMu = covariance_obj.get_GRM(window_beginning=window_starts[w], window_end=window_ends[w],
                                               variants=variants, inds=inds)
        if tmpCov is not None:
            covariance_obj.write(out=outname, inds=inds, covariances_picklefile=covariances_picklefile, logfile=logfile)
            for m in AIM_methods:
                m.run_association(index=w, out=outname)
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
                                      logfile, outname, population_structure):
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

    # log progress
    start = time.time()

    # windows are given by trees
    if window_size is None:
        window_index = 0
        for tree in trees.trees():
            tree_obj = tt.TTree(tree)

            # print("----------------")
            # print("tree_obj.start", tree_obj.start, ", tree_obj.end", tree_obj.end, ", window_index", window_index)

            if tree_obj.is_testable(skip_first_tree):
                if window_size is None:
                    covariance_obj.add_tree(tree_obj=tree_obj, inds=inds, proportion=1.0)
                    write_and_test_window_for_association(covariance_obj=covariance_obj,
                                                          inds=inds,
                                                          AIM_methods=AIM_methods,
                                                          outname=outname,
                                                          window_index=window_index,
                                                          covariances_picklefile=covariances_picklefile)

            window_index += 1
            # log progress
            if tree.index % 1 == 0:
                end = time.time()
                logfile.info("- Ran AIM for " + str(tree.index) + " trees in " + str(round(end - start)) + " s")

    # there is a window size
    else:
        window_index = 0
        for tree in trees.trees():
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
                            or (proportion == 1.0 and tree_obj.end >= window_ends[0]):
                        write_and_test_window_for_association(covariance_obj=covariance_obj,
                                                              inds=inds,
                                                              AIM_methods=AIM_methods,
                                                              outname=outname,
                                                              window_index=window_index,
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
                if tree.index % 1 == 0:
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


def run_association_AIM(trees, inds, variants, pheno, args, ass_method, window_size,
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
    num_tests = trees.num_trees

    # remove the next start that is included, i think this tree is removed due to incompleteness when taking tree subset
    window_starts = trees.breakpoints(as_array=True)[0:trees.num_trees]
    window_ends = trees.breakpoints(as_array=True)[1:]

    if window_size is not None:
        window_ends = get_window_ends(window_size=window_size, trees_interval=args.trees_interval)
        window_starts = [x - window_size for x in window_ends]
        num_tests = len(window_ends)

    # initialize and write phenotypes
    if covariance == "eGRM" or covariance == "GRM":
        pheno.set_missing_phenotype_status(inds)
        pheno.write_to_file_gcta_eGRM(inds=inds, out=outname, logfile=logfile)
    else:
        pheno.write_to_file_gcta_scaled(out=outname, logfile=logfile)
    if args.simulate_phenotypes is True:
        pheno.find_causal_trees(trees)
        pheno.find_causal_windows(window_ends=window_ends, window_starts=window_starts)

    # write GCTA files and scripts
    if args.population_structure is not None:
        with open(outname + '_multi_grm.txt', 'w') as f:
            f.write(outname + '\n')
            f.write(args.population_structure + '\n')

    pheno_file = outname + "_phenotypes.phen"
    # if args.simulate_phenotypes:
    #     pheno_file = args.out + "_phenotypes.phen"

    # create association method objects
    logfile.info("- Running associations tests using test methods " + str(
        args.AIM_method) + " for a sequence of trees")
    AIM_methods = []
    for m in args.AIM_method:
        with open(outname + "_run_gcta_" + m + ".sh", 'w') as f:
            f.write("#!/bin/bash\n")
            if args.population_structure:
                if m == "REML":
                    f.write(args.GCTA + " --reml --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out "
                            + outname + "_REML --reml-lrt 1 --threads " + str(args.num_gcta_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
                elif m == "HE":
                    f.write(
                        args.GCTA + " --HEreg --mgrm " + outname + "_multi_grm.txt --pheno " + pheno_file + " --out "
                        + outname + "_HE --reml-lrt 1 --threads " + str(args.num_gcta_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
                    # grep results
                    f.write("sed -n '2,6p' " + outname + "_" + m + ".HEreg | unexpand -a | tr -s \'\t\' > "
                            + outname + "_HE-CP_result.txt\n")
                    f.write("sed -n '9,13p' " + outname + "_" + m + ".HEreg | unexpand -a | tr -s \'\t\' > "
                            + outname + "_HE-SD_result.txt\n")

            else:
                if m == "REML":
                    f.write(args.GCTA + " --reml --grm " + outname + " --pheno " + pheno_file + " --out " + outname
                            + "_REML --threads " + str(args.num_gcta_threads) + " --reml-maxit 500  > " + outname + "_tmp.out\n")
                elif m == "HE":
                    f.write(args.GCTA + " --HEreg --grm " + outname + " --pheno " + pheno_file + " --out " + outname
                            + "_HE --threads " + str(args.num_gcta_threads) + " --reml-maxit 500 > " + outname + "_tmp.out\n")
                    # grep results
                    f.write("sed -n '2,4p' " + outname  + "_" + m + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                            + outname + "_HE-CP_result.txt\n")
                    f.write("sed -n '7,9p' " + outname + "_" + m + ".HEreg | unexpand -a | tr -s \'\\t\' > "
                            + outname + "_HE-SD_result.txt\n")

        st = os.stat(outname + "_run_gcta_" + m + ".sh")
        os.chmod(outname + "_run_gcta_" + m + ".sh", st.st_mode | stat.S_IEXEC)

        test_obj = get_AIM_test_object(m, phenotypes=pheno, num_associations=num_tests)
        AIM_methods.append(test_obj)

    # create covariance type object
    covariance_obj = cov.get_covariance_object(covariance)

    # ----------------
    # run association tests
    # ----------------

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

    # tree based covariance
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
                                          population_structure=args.population_structure)

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
