#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:37:12 2021

@author: linkv
"""
import numpy as np
import statsmodels.api as sm
import scipy
import TTree as tt
# from limix_lmm.lmm_core import LMMCore
import utils as ut
import time
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
import os
import sys
import TCovariance as cov


def OLS(genotypes, phenotypes):
    # add intercept
    genotypes_test = sm.tools.add_constant(genotypes)
    PVALUE = sm.OLS(phenotypes.y, genotypes_test).fit().pvalues[1]
    return PVALUE


def run_association_GWAS(trees, inds, variants, pheno, args, impute, logfile):
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
        GWAS = TAssociationTestingGWAS(phenotypes=pheno, num_typed_variants=gt_matrix_imputed.shape[1])
        GWAS.test_with_positions_from_X_matrix(X=gt_matrix_imputed, positions=pos,
                                               variants_sample=variants,
                                               logfile=logfile)
        GWAS.write_to_file_with_X_matrix(positions=pos, name=outname, logfile=logfile)

    else:
        logfile.info("- Using genotypes from tree file for GWAS:")
        # run association tests
        GWAS = TAssociationTestingGWAS(phenotypes=pheno, num_typed_variants=variants.num_typed)
        GWAS.test_with_variants_object(variants, inds, logfile)
        GWAS.write_to_file(variants, outname, logfile)
        # GWAS.manhattan_plot(variant_positions=variants.info['position'], plots_dir=plots_dir)


def get_AIM_test_object(test_name, phenotypes, num_associations):
    if test_name == "HE":
        test_obj = TAssociationTestingRegionsGCTA_HE(phenotypes, num_associations)
    elif test_name == "REML":
        test_obj = TAssociationTestingRegionsGCTA_REML(phenotypes, num_associations)
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


def write_and_test_window_for_association(covariance_obj, inds, AIM_methods, outname, window_index):
    if covariance_obj.covariance_matrix_haploid is None:
        raise ValueError("trying to test empty covariance matrix for association at window index " + str(window_index))
    covariance_obj.finalize(inds=inds)
    written = covariance_obj.write(out=outname, inds=inds)
    if written:
        for m in AIM_methods:
            m.run_association(index=window_index, out=outname)
        covariance_obj.clear()
    else:
        raise ValueError(
            "covariance matrix " + covariance_obj.covariance_type + " was not written to file")


def run_variant_based_covariance_testing(covariance_obj, AIM_methods, variants, window_ends, window_starts, num_tests,
                                         inds, logfile, outname):
    """
    Write covariance calculated based on variants within a window (can be one tree) to file and test it for association with
    phenotypes. Currently, the only covariance type based on variants is GRM.

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

    for w in range(num_tests):
        tmpCov, tmpMu = covariance_obj.get_GRM(window_beginning=window_starts[w], window_end=window_ends[w], variants=variants,
                               inds=inds)
        if tmpCov is not None:
            written = covariance_obj.write(out=outname, inds=inds, logfile=logfile)
            if written:
                for m in AIM_methods:
                    m.run_association(index=w, out=outname)
            else:
                raise ValueError("Covariance matrix " + covariance_obj.covariance_type + " was not written to file")
            covariance_obj.clear()

        # log progress
        if w % 10 == 0:
            end = time.time()
            logfile.info("- Ran AIM for " + str(w) + " windows in " + str(round(end - start)) + " s")

    for m in AIM_methods:
        m.write_to_file(window_starts=window_starts_copy,
                        window_ends=window_ends_copy,
                        out=outname,
                        logfile=logfile)


def run_tree_based_covariance_testing(trees, covariance_obj, AIM_methods, window_ends, window_starts,
                                      window_size, skip_first_tree, inds, logfile, outname):

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
                                                          window_index=window_index)

            window_index += 1

    # there is a window size
    else:
        window_index = 0
        for tree in trees.trees():
            tree_obj = tt.TTree(tree)

            # print("----------------")
            # print("tree_obj.start", tree_obj.start, ", tree_obj.end", tree_obj.end, ", window_index", window_index,
            #       ", window_starts[0]", window_starts[0], ", window_ends[0]", window_ends[0])

            if tree_obj.is_testable(skip_first_tree):
                # calculate one covariance matrix per tree
                if window_size is None:
                    covariance_obj.add_tree(tree_obj=tree_obj, inds=inds, proportion=1.0)
                    write_and_test_window_for_association(covariance_obj=covariance_obj,
                                                          inds=inds,
                                                          AIM_methods=AIM_methods,
                                                          outname=outname,
                                                          window_index=window_index)
                    window_index += 1

                # calculate one covariance matrix per region
                else:
                    # increase window index until we get to the one the new tree starts in
                    if tree_obj.start >= window_ends[0]:
                        raise ValueError("Tree start is >= window end. This should not happen")

                    # tree start is definitely within window based on previous while loop
                    proportion = get_proportion_of_tree_within_window(window_start=window_starts[0],
                                                                      window_end=window_ends[0],
                                                                      tree_start=tree_obj.start,
                                                                      tree_end=tree_obj.end)
                    # print("tree is in window with index", window_index, "with proportion", proportion)

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
                                                                  window_index=window_index)

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
                if tree.index % 10 == 0:
                    end = time.time()
                    logfile.info("- Ran AIM for " + str(tree.index) + " trees in " + str(round(end - start)) + " s")
            # else:
            #     print("tree was not usable")

    # write association test results to file
    for m in AIM_methods:
        m.write_to_file(window_starts=window_starts_copy,
                        window_ends=window_ends_copy,
                        out=outname,
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
    pheno.find_causal_trees(trees)
    pheno.find_causal_windows(window_ends=window_ends, window_starts=window_starts)
    if covariance == "eGRM" or covariance == "GRM":
        pheno.write_to_file_gcta_eGRM(inds=inds, out=outname, logfile=logfile)
    else:
        pheno.write_to_file_gcta_scaled(out=outname, logfile=logfile)

    # create association method objects
    logfile.info("- Running associations tests using test methods " + str(
        args.AIM_method) + " for a sequence of trees")
    AIM_methods = []
    for m in args.AIM_method:
        test_obj = get_AIM_test_object(m, phenotypes=pheno, num_associations=num_tests)
        AIM_methods.append(test_obj)

    # create covariance type object
    covariance_obj = cov.get_covariance_object(covariance)

    # ----------------
    # run association tests
    # ----------------

    # variant based covariance
    if covariance == "GRM":
        run_variant_based_covariance_testing(covariance_obj=covariance_obj,
                                             AIM_methods=AIM_methods,
                                             variants=variants,
                                             window_ends=window_ends,
                                             window_starts=window_starts,
                                             num_tests=num_tests,
                                             inds=inds,
                                             logfile=logfile,
                                             outname=outname)

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
                                          logfile=logfile,
                                          outname=outname)

    logfile.sub()


class TAssociationTesting:
    def __init__(self, phenotypes):
        self.phenotypes = phenotypes
        self.num_associations = -1
        self.p_values = np.empty(0)
        self.p_values.fill(np.nan)
        # self.p_values_init = False

    # def _check_compatibility(self, ts_object, phenotypes):
    # #check if ts_object is compatible with phenotype
    # if phenotypes.num_variants != self.num_associations:
    #   raise ValueError("Phenotype object must contain same number
    #                      of sampled variants as ts_object to create a GWAS object. Phenotype: " +
    #                      str(phenotypes.num_variants) + " and ts_object: " + str(self.num_associations))
    # if phenotypes.N != ts_object.num_samples:
    #   raise ValueError( "Phenotype object must contain same number of samples as ts_object to create a GWAS object.
    #                       Phenotype: " + str(phenotypes.N) + " and ts_object: " + str(ts_object.num_samples))

    def p_value_dist(self, subplot, num_bins):
        subplot.hist(self.p_values, bins=num_bins)
        subplot.set(xlabel='p-values', ylabel='density', title=self.phenotypes.name)

    def chi_squared(self, num_bins):
        h, bin_edges = np.histogram(self.p_values, num_bins)
        tmp = scipy.stats.chisquare(f_obs=h, f_exp=((len(self.p_values)) / (num_bins)))
        return tmp


class TAssociationTestingGWAS(TAssociationTesting):
    """
    SNP based association testing using Ordinary Least Squares regression
    """

    def __init__(self, phenotypes, num_typed_variants):

        super().__init__(phenotypes)

        self.num_typed_variants = num_typed_variants
        self.num_associations = self.num_typed_variants
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        self.imputed_status = np.repeat(False, self.num_associations)

    def test_with_variants_object(self, variants, inds, logfile):
        # counter respective to typed variants
        i = 0
        for v, variant in enumerate(variants.variants):
            if variants.info.iloc[v]['typed']:
                if inds.ploidy == 2:
                    genotypes = inds.get_diploid_genotypes(variant.genotypes)
                else:
                    genotypes = variant.genotypes

                PVALUE = OLS(genotypes=genotypes, phenotypes=self.phenotypes)
                self.p_values[i] = PVALUE
                i += 1
        logfile.info("- Ran OLS for " + str(variants.num_typed) + " variants")

    def test_with_positions_from_X_matrix(self, X, positions, variants_sample, logfile):
        # counter respective to typed variants
        if len(positions) != X.shape[1]:
            raise ValueError("X genotype matrix does not have same number of columns (" + str(X.shape[1])
                             + ") as positions (" + str(len(positions)) + ")")
        positions_sample = set(variants_sample.info['position'].values)
        for v in range(X.shape[1]):
            genotypes = X[:, v]
            PVALUE = OLS(genotypes=genotypes, phenotypes=self.phenotypes)
            self.p_values[v] = PVALUE
            if not (positions[v] in positions_sample):
                self.imputed_status[v] = True
        logfile.info("- Ran OLS for " + str(X.shape[1]) + " variants")

    def write_to_file_with_X_matrix(self, positions, name, logfile):
        # results for each variant
        table = pd.DataFrame()
        table['start'] = positions
        table['end'] = positions
        # table['typed'] = variants.info['typed']
        table['p_value'] = self.p_values
        table['imputed_status'] = self.imputed_status
        # table['causal'] = np.repeat("FALSE", self.num_associations)
        # table.loc[self.phenotypes.causal_variant_indeces, 'causal'] = "TRUE"
        # table['betas'] = self.phenotypes.betas
        logfile.info("- Writing results from OLS to '" + name + "_variants_results.csv'")
        table.to_csv(name + "_variants_results.csv", index=False, header=True)

        # summary statistics
        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        logfile.info("- Writing stats from OLS to '" + name + "_variants_stats.csv'")
        stats.to_csv(name + "_variants_stats.csv", index=False, header=True)

    def write_to_file(self, variants, name, logfile):
        # results for each variant
        info_typed = variants.info.loc[variants.info['typed'] == True]
        info_typed['index'] = range(variants.num_typed)
        info_typed.set_index(info_typed['index'], drop=True, inplace=True)

        table = pd.DataFrame()
        table['start'] = info_typed['position']
        table['end'] = info_typed['position']
        # table['typed'] = variants.info['typed']
        table['p_value'] = self.p_values
        table['imputed_status'] = self.imputed_status
        # table['causal'] = np.repeat("FALSE", self.num_associations)
        # table.loc[self.phenotypes.causal_variant_indeces, 'causal'] = "TRUE"
        # table['betas'] = self.phenotypes.betas 
        logfile.info("- Writing results from OLS to '" + name + "_variants_results.csv'")
        table.to_csv(name + "_variants_results.csv", index=False, header=True)

        # summary statistics
        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        logfile.info("- Writing stats from OLS to '" + name + "_variants_stats.csv'")
        stats.to_csv(name + "_variants_stats.csv", index=False, header=True)

    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, size=1, n_snps_lowess=0,
                              *args):
        """
        Parameters
        ----------
        variant_positions : float?
            Positions of the variants in the genome.
        subplot : matplotlib.axes._subplots.AxesSubplot
            Which subplot is used for plotting.
        index_min : int
            Min index of variants that should be plotted. Refers to index within TVariants object.
        index_max : int
            Max index of variants that should be plotted. Refers to index within TVariants object.
        size : int, optional
            Size of the points in scatter plot. The default is 1.
        n_snps_lowess : int, optional
            Number of SNPs that are used to calculate lowess smoothing. The default is 0, meaning that no lowess curve is plotted
        *args : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            The variant index within the variant class (so after filtering) must be between zero and the max.

        Returns
        -------
        None.

        """
        if index_min < 0 or index_max < 0:
            raise ValueError("data subset indeces for manhattan plotting must be positive")
        if index_max > len(self.p_values) or index_min > len(self.p_values):
            raise ValueError("data subset index cannot be larger than number of p-values")

        p_values = self.p_values
        p_values[(np.where(p_values == 0.0))] = np.nextafter(0, 1)
        q_values = -np.log10(p_values)

        subplot.scatter(variant_positions[index_min:index_max], q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='variant position', ylabel='q-value', title="GWAS")
        for v, var in enumerate(self.phenotypes.causal_variants):
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            colscale = self.phenotypes.causal_power[v]
            subplot.axvline(x=var.site.position, color=str(0.3), alpha=0.5, lw=colscale * 100)
            subplot.axvline(x=var.site.position, color="black", lw=0.5)

        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess / len(q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=q_values[index_min:index_max],
                                          exog=variant_positions[index_min:index_max], frac=fraction,
                                          return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")

        subplot.axhline(y=8, color="red", lw=0.5)

    def manhattan_plot(self, variant_positions, plots_dir, *args):
        fig, ax = plt.subplots(1, figsize=(10, 10))
        self.manhattan_plot_subset(variant_positions=variant_positions, subplot=ax, index_min=0,
                                   index_max=len(self.p_values), *args)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(plots_dir + 'OLS_GWAS.png', bbox_inches='tight')


class TAssociationTestingRegions(TAssociationTesting):
    """
    base class for all tree-based association tests
    """

    def __init__(self, phenotypes, num_associations):

        super().__init__(phenotypes)

        self.num_associations = num_associations
        # self._check_compatibility(ts_object, phenotypes)

    def manhattan_plot(self, variant_positions, subplot, logfile, *args):
        self.manhattan_plot_subset(variant_positions=variant_positions, subplot=subplot, index_min=0,
                                   index_max=self.num_associations, p_values=self.p_values)

    def manhattan_plot_special_pvalues(self, variant_positions, p_values, subplot, logfile, title_supplement="", *args):
        logfile.info("Plotting " + str(self.num_associations) + " associations")
        self.manhattan_plot_subset(variant_positions=variant_positions, subplot=subplot, index_min=0,
                                   index_max=self.num_associations, p_values=p_values,
                                   title_supplement=title_supplement)

    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, p_values,
                              title_supplement="", size=1, n_snps_lowess=0, *args):
        """
        Parameters
        ----------
        variant_positions : vector of floats
            Positions of the variants in the genome.
        subplot : matplotlib.axes._subplots.AxesSubplot
            Which subplot is used for plotting.
        index_min : int
            Min index of variants that should be plotted. Refers to index within TVariants object.
        index_max : int
            Max index of variants that should be plotted. Refers to index within TVariants object.
        p_values: np.array
            p-values to be plotted
        title_supplement: str
            Supplementary text to print in plot title
        size : int, optional
            Size of the points in scatter plot. The default is 1.
        n_snps_lowess : int, optional
            Number of SNPs that are used to calculate lowess smoothing. The default is 0, meaning that no lowess curve is plotted
        *args : TYPE
            DESCRIPTION.

        Raises
        ------
        ValueError
            The variant index within the variant class (so after filtering) must be between zero and the max.

        Returns
        -------
        None.

        """

        if index_min < 0 or index_max < 0:
            raise ValueError("data subset indeces for manhattan plotting must be positive")
        if index_max > self.num_associations or index_min > self.num_associations:
            raise ValueError("data subset index cannot be larger than number of p-values")

        p_values[(np.where(p_values == 0.0))] = np.nextafter(0, 1)
        q_values = -np.log10(p_values)

        subplot.scatter(variant_positions[index_min:index_max], q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='tree index', ylabel='q-value', title=self.phenotypes.name + title_supplement)
        for t in self.phenotypes.causal_tree_indeces:
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            # colscale = self.phenotypes.causal_power[v] 
            # subplot.axvline(x=t, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=t, color="black", lw=0.5)

        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess / len(q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=q_values[index_min:index_max],
                                          exog=variant_positions[index_min:index_max], frac=fraction,
                                          return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")

        subplot.axhline(y=8, color="red", lw=0.5)


class TAssociationTestingRegionsGCTA(TAssociationTestingRegions):
    """
    tree-based asssociation testing using GCTA 
    """

    def __init__(self, phenotypes, num_associations):
        super().__init__(phenotypes, num_associations)

    def run_association(self, index, out):
        self.run_association_one_window_gcta(index=index, out=out)

    def run_association_one_window_gcta(self, index, out):
        raise ValueError("function run_association_one_window_gcta not implemented in base class")


class TAssociationTestingRegionsGCTA_HE(TAssociationTestingRegionsGCTA):
    """
    tree-based association testing using GCTA Haseman-Elston algorithm
    """

    def __init__(self, phenotypes, num_associations):

        super().__init__(phenotypes, num_associations)

        # p-value containers
        self.p_values_HECP_OLS = np.empty(self.num_associations)
        self.p_values_HECP_OLS.fill(np.nan)
        self.p_values_HECP_Jackknife = np.empty(self.num_associations)
        self.p_values_HECP_Jackknife.fill(np.nan)
        self.p_values_HESD_OLS = np.empty(self.num_associations)
        self.p_values_HESD_OLS.fill(np.nan)
        self.p_values_HESD_Jackknife = np.empty(self.num_associations)
        self.p_values_HESD_Jackknife.fill(np.nan)

        # other statistics
        self.V_G_over_Vp_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_HECP.fill(np.nan)
        self.V_G_over_Vp_HESD = np.empty(self.num_associations)
        self.V_G_over_Vp_HESD.fill(np.nan)
        self.V_G_over_Vp_SE_OLS_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_OLS_HECP.fill(np.nan)
        self.V_G_over_Vp_SE_OLS_HESD = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_OLS_HESD.fill(np.nan)
        self.V_G_over_Vp_SE_Jackknife_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_Jackknife_HECP.fill(np.nan)
        self.V_G_over_Vp_SE_Jackknife_HESD = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_Jackknife_HESD.fill(np.nan)

    def run_association(self, index, out):
        self.run_association_one_window_gcta(index=index, out=out)

    def run_association_one_window_gcta(self, index, out):
        # create gcta input files, run gcta and parse output
        exit_code = subprocess.call([os.path.dirname(sys.argv[0]) + "/run_gcta_HE.sh", out])
        # exit_code = subprocess.call([os.getcwd() + "/run_gcta_HE.sh", out])

        # read results
        HE_CP = pd.read_table(out + "_HE-CP_result.txt")
        HE_SD = pd.read_table(out + "_HE-SD_result.txt")

        # p-values
        self.p_values_HECP_OLS[index] = HE_CP["P_OLS"][1]
        if HE_CP["P_OLS"][1] < 0:
            raise ValueError("window index", index, "produced negative p-value for CP OLS")

        self.p_values_HECP_Jackknife[index] = HE_CP["P_Jackknife"][1]
        if HE_CP["P_Jackknife"][1] < 0:
            raise ValueError("window index", index, "produced negative p-value for CP Jackknife")

        self.p_values_HESD_OLS[index] = HE_SD["P_OLS"][1]
        if HE_SD["P_OLS"][1] < 0:
            raise ValueError("window index", index, "produced negative p-value for SD OLS")

        self.p_values_HESD_Jackknife[index] = HE_SD["P_Jackknife"][1]
        if HE_SD["P_Jackknife"][1] < 0:
            raise ValueError("window index", index, "produced negative p-value for SD Jackknife")

        # other statistics
        self.V_G_over_Vp_HECP[index] = HE_CP["Estimate"][1]
        self.V_G_over_Vp_HESD[index] = HE_SD["Estimate"][1]

        self.V_G_over_Vp_SE_OLS_HECP[index] = HE_CP["SE_OLS"][1]
        self.V_G_over_Vp_SE_OLS_HESD[index] = HE_SD["SE_OLS"][1]
        self.V_G_over_Vp_SE_Jackknife_HECP[index] = HE_CP["SE_Jackknife"][1]
        self.V_G_over_Vp_SE_Jackknife_HESD[index] = HE_SD["SE_Jackknife"][1]

    def write_to_file(self, window_starts, window_ends, out, logfile):
        table = pd.DataFrame()
        table['start'] = window_starts
        table['end'] = window_ends

        # p-values
        print("self.p_values_HECP_OLS", self.p_values_HECP_OLS)
        table['p_values_HECP_OLS'] = self.p_values_HECP_OLS
        table['p_values_HECP_Jackknife'] = self.p_values_HECP_Jackknife
        table['p_values_HESD_OLS'] = self.p_values_HESD_OLS
        table['p_values_HESD_Jackknife'] = self.p_values_HESD_Jackknife

        # other stats
        table['V_G_over_Vp_HECP'] = self.V_G_over_Vp_HECP
        table['V_G_over_Vp_HESD'] = self.V_G_over_Vp_HESD
        table['V_G_over_Vp_SE_OLS_HECP'] = self.V_G_over_Vp_SE_OLS_HECP
        table['V_G_over_Vp_SE_OLS_HESD'] = self.V_G_over_Vp_SE_OLS_HESD
        table['V_G_over_Vp_SE_Jackknife_HECP'] = self.V_G_over_Vp_SE_Jackknife_HECP
        table['V_G_over_Vp_SE_Jackknife_HESD'] = self.V_G_over_Vp_SE_Jackknife_HESD

        # causal or not
        table['causal'] = np.repeat("FALSE", self.num_associations)
        table.loc[self.phenotypes.causal_window_indeces, 'causal'] = "TRUE"
        # table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"

        table.to_csv(out + "_trees_HE_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + out + "_trees_HE_results.csv'")

        stats = pd.DataFrame({'min_p_value_HECP_OLS': [np.nanmin(self.p_values_HECP_OLS)],
                              'min_p_value_HECP_Jackknife': [np.nanmin(self.p_values_HECP_Jackknife)],
                              'min_p_value_HESD_OLS': [np.nanmin(self.p_values_HESD_OLS)],
                              'min_p_value_HESD_Jackknife': [np.nanmin(self.p_values_HESD_Jackknife)],
                              'max_p_value_HECP_OLS': [np.nanmax(self.p_values_HECP_OLS)],
                              'max_p_value_HECP_Jackknife': [np.nanmax(self.p_values_HECP_Jackknife)],
                              'max_p_value_HESD_OLS': [np.nanmax(self.p_values_HESD_OLS)],
                              'max_p_value_HESD_Jackknife': [np.nanmax(self.p_values_HESD_Jackknife)]
                              })
        stats.to_csv(out + "_trees_HE_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from HE to '" + out + "_trees_HE_stats.csv'")


class TAssociationTestingRegionsGCTA_REML(TAssociationTestingRegionsGCTA):
    """
    tree-based association testing using CGTA REML algorithm
    """

    def __init__(self, ts_object, phenotypes):

        super().__init__(ts_object, phenotypes)

        # results containers
        self.p_values = np.empty(self.num_associations)
        self.p_values.fill(np.nan)
        self.V_G = np.empty(self.num_associations)
        self.V_G.fill(np.nan)
        self.V_e = np.empty(self.num_associations)
        self.V_e.fill(np.nan)
        self.Vp = np.empty(self.num_associations)
        self.Vp.fill(np.nan)
        self.V_G_over_Vp = np.empty(self.num_associations)
        self.V_G_over_Vp.fill(np.nan)
        self.logL = np.empty(self.num_associations)
        self.logL.fill(np.nan)
        self.logL0 = np.empty(self.num_associations)
        self.logL0.fill(np.nan)
        self.LRT = np.empty(self.num_associations)
        self.LRT.fill(np.nan)

        self.V_G_SE = np.empty(self.num_associations)
        self.V_G_SE.fill(np.nan)
        self.V_e_SE = np.empty(self.num_associations)
        self.V_e_SE.fill(np.nan)
        self.Vp_SE = np.empty(self.num_associations)
        self.Vp_SE.fill(np.nan)
        self.V_G_over_Vp_SE = np.empty(self.num_associations)
        self.V_G_over_Vp_SE.fill(np.nan)

    def run_association(self, index, out):
        self.run_association_one_window_gcta(index=index, out=out)

    def run_association_one_window_gcta(self, index, out):
        # create gcta input files, run gcta and parse output
        exit_code = subprocess.call([os.path.dirname(sys.argv[0]) + "/run_gcta_REML.sh", out])

        # read results
        result = pd.read_table(out + "_REML.hsq")
        result_pvalue = float(result['Variance'][result['Source'] == 'Pval'])
        if result_pvalue < 0:
            raise ValueError("Negative p-value for window with index " + str(index))
        if result_pvalue > 1:
            raise ValueError("p-value larger than 1 for window with index " + str(index))

        self.p_values[index] = result_pvalue
        if result_pvalue < 0:
            raise ValueError("window index", index, "produced negative p-value with REML")

        self.V_G[index] = float(result['Variance'][result['Source'] == 'V(G)'])
        self.V_e[index] = float(result['Variance'][result['Source'] == 'V(e)'])
        self.Vp[index] = float(result['Variance'][result['Source'] == 'Vp'])
        self.V_G_over_Vp[index] = float(result['Variance'][result['Source'] == 'V(G)/Vp'])
        self.logL[index] = float(result['Variance'][result['Source'] == 'logL'])
        self.logL0[index] = float(result['Variance'][result['Source'] == 'logL0'])
        self.LRT[index] = float(result['Variance'][result['Source'] == 'LRT'])

        self.V_G_SE[index] = float(result['SE'][result['Source'] == 'V(G)'])
        self.V_e_SE[index] = float(result['SE'][result['Source'] == 'V(e)'])
        self.Vp_SE[index] = float(result['SE'][result['Source'] == 'Vp'])
        self.V_G_over_Vp_SE[index] = float(result['SE'][result['Source'] == 'V(G)/Vp'])

    def write_to_file(self, window_starts, window_ends, out, logfile):
        table = pd.DataFrame()
        table['start'] = window_starts
        table['end'] = window_ends
        table['p_values'] = self.p_values
        print("self.p_values", self.p_values)
        table['V_G'] = self.V_G
        table['V_e'] = self.V_e
        table['Vp'] = self.Vp
        table['V_G_over_Vp'] = self.V_G_over_Vp
        table['logL'] = self.logL
        table['logL0'] = self.logL0
        table['LRT'] = self.LRT
        table['V_G_SE'] = self.V_G_SE
        table['V_e_SE'] = self.V_e_SE
        table['Vp_SE'] = self.Vp_SE
        table['V_G_over_Vp_SE'] = self.V_G_over_Vp_SE

        table['causal'] = np.repeat("FALSE", self.num_associations)
        table.loc[self.phenotypes.causal_window_indeces, 'causal'] = "TRUE"

        table.to_csv(out + "_trees_REML_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + out + "_trees_REML_results.csv'")

        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        stats.to_csv(out + "_trees_REML_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from tree association tests to '" + out + "_trees_REML_stats.csv'")


class TTreeAssociationMantel(TAssociationTestingRegions):

    def __init__(self, ts_object, phenotypes):

        super().__init__(ts_object, phenotypes)

        # p-value containers
        self.p_values = np.empty(self.num_associations)

    def run_Mantel(self, ts_object, phenotypes, inds):
        # test for associations
        diffs = phenotypes.diffs()
        start = time.time()
        for tree in ts_object.trees():
            if tree.index % 100 == 0:
                end = time.time()
                print("Ran Mantel for", tree.index, "trees in ", round(end - start), "s")
            if tree.total_branch_length == 0:
                print("tree's total branch length is zero")
                continue
            tree_obj = tt.TTree(tree)
            tmrca = tree_obj.TMRCA(inds.num_haplotypes)
            # print("tmrca",tmrca)
            self.p_values[tree.index] = ut.mantel(tmrca, diffs)
            if self.p_values[tree.index] < 0:
                print(tmrca)
                raise ValueError("p-value is negative")

# def runLimix(self, ts_object, N, y, F, random):
# self.lrt = np.empty(self.num_associations)

#     raise ValueError("Limix not currently implemented")

#     G = np.zeros(N).reshape(N,1)
#     # G = np.random.binomial(1, 0.5, N).reshape(N,1)
#     # Inter = np.zeros(N).reshape(N,1)
#     Inter = None

#     for tree in ts_object.trees():
#         # if tree.index == ts_object.num_trees-1:
#         #     continue
#         # if tree.index % 1000 == 0:
#         print("tree index: ",tree.index)
#         tree_obj = tt.TTree(tree, N)
#         lmm = LMMCore(y, F, tree_obj.solving_function)
#         lmm.process(G, Inter) #this needs step param to produce only one p-value per tree. for this i need the number of sites per tree, or just use 1?
#         self.p_values[tree.index] = lmm.getPv()
#         print("p-value", self.p_values[tree.index])
#         # raise ValueError("printing covariance")
#         # beta = lmm.getBetaSNP()
#         # beta_ste = lmm.getBetaSNPste()
#         # self.lrt[tree.index] = lmm.getLRT() #likelihood ratio


# -----------------
# example_tree = trees.aslist()[10995]
# tree_obj = tt.TTree(example_tree, N)

# tmrca = np.zeros([N, N])
# height = 0
# for c in example_tree.nodes():
#     print("c",c)
#     descendants = list(example_tree.samples(c))
#     n = len(descendants)
#     if(n == 0 or n == N or example_tree.time(c) == 0): #The branch length for a node that has no parent (e.g., a root) is defined as zero.
#         continue
#     t = example_tree.time(example_tree.parent(c)) - example_tree.time(c)
#     tmrca[np.ix_(descendants, descendants)] -= t
#     height = max(height, example_tree.time(example_tree.parent(c))) #time returns the time of a node
# tmrca += height
# # covariance = (tmrca+tmrca.T)/2 #why does caoqi do this??
# np.fill_diagonal(covariance, 0)

# #test if matrix is positive semidefinite
# np.linalg.cholesky(covariance)
# np.exp(-covariance)
# inv = np.linalg.inv(covariance)
# tmp = np.dot(inv, array)
# # print("shape of my dot product",np.shape(tmp))
# tree_obj = tt.TTree(example_tree, N)
# F = sp.zeros(N)
# F.reshape(N,1)
# y = pheno_random.y.reshape(N,1)
# lmm = LMMCore(y, F.reshape(N,1), tree_obj.solving_function)
# Inter = sp.zeros(N).reshape(N,1)
# G = sp.zeros(N).reshape(N,1)
# lmm.process(G, Inter)
# lmm.getPv()
