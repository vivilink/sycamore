#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:37:12 2021
@author: linkv
"""
import numpy as np
import statsmodels.api as sm
import scipy
import utils as ut
import time
import pandas as pd
import matplotlib.pyplot as plt
import subprocess
# from limix_lmm.lmm_core import LMMCore
from glimix_core.lmm import LMM
from numpy_sugar.linalg import economic_qs
from scipy import stats
import association_functions as af
import TTree as tt


class TAssociationTesting:
    def __init__(self, phenotypes):
        self.name = "base"
        # self.phenotypes = phenotypes
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

    def p_value_dist(self, subplot, num_bins, phenotypes):
        subplot.hist(self.p_values, bins=num_bins)
        subplot.set(xlabel='p-values', ylabel='density', title=phenotypes.name)

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

        self.name = "GWAS"
        self.num_typed_variants = num_typed_variants
        self.num_associations = self.num_typed_variants
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        self.imputed_status = np.repeat(False, self.num_associations)

    def test_with_variants_object(self, variants, phenotypes, inds, logfile):
        # counter respective to typed variants
        i = 0
        for v, variant in enumerate(variants.variants):
            if variants.info.iloc[v]['typed']:
                if inds.ploidy == 2:
                    genotypes = inds.get_diploid_genotypes(variant.genotypes)
                else:
                    genotypes = variant.genotypes
                #
                # if len(genotypes) != len(phenotypes.y):
                #     # TODO: remove this after debugging on real data
                #     raise ValueError("Genotypes length (" + str(len(genotypes)) + ") is not same as phenotypes length ("
                #                      + str(len(phenotypes.y)) + ")")

                PVALUE = af.OLS(genotypes=genotypes, phenotypes=phenotypes.y)
                self.p_values[i] = PVALUE
                i += 1
        logfile.info("- Ran OLS for " + str(variants.num_typed) + " variants")

    def test_with_positions_from_X_matrix(self, X, positions, variants_sample, phenotypes, logfile):
        # counter respective to typed variants
        if len(positions) != X.shape[1]:
            raise ValueError("X genotype matrix does not have same number of columns (" + str(X.shape[1])
                             + ") as positions (" + str(len(positions)) + ")")
        positions_sample = set(variants_sample.info['position'].values)
        for v in range(X.shape[1]):
            genotypes = X[:, v]
            PVALUE = af.OLS(genotypes=genotypes, phenotypes=phenotypes)
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

    def manhattan_plot_subset(self, variant_positions, phenotypes, subplot, index_min, index_max, size=1,
                              n_snps_lowess=0,
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
        for v, var in enumerate(phenotypes.causal_variants):
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            colscale = phenotypes.causal_power[v]
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

        self.name = "regions_base"
        self.num_associations = num_associations
        # self._check_compatibility(ts_object, phenotypes)

    def manhattan_plot(self, variant_positions, subplot, phenotypes, *args):
        self.manhattan_plot_subset(variant_positions=variant_positions, subplot=subplot, index_min=0,
                                   index_max=self.num_associations, p_values=self.p_values, phenotypes=phenotypes)

    def manhattan_plot_special_pvalues(self, variant_positions, p_values, subplot, logfile, phenotypes,
                                       title_supplement="", *args):
        logfile.info("Plotting " + str(self.num_associations) + " associations")
        self.manhattan_plot_subset(variant_positions=variant_positions, subplot=subplot, index_min=0,
                                   index_max=self.num_associations, p_values=p_values,
                                   title_supplement=title_supplement, phenotypes=phenotypes)

    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, p_values, phenotypes,
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
        phenotypes : TPhenotypes
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
        subplot.set(xlabel='tree index', ylabel='q-value', title=phenotypes.name + title_supplement)
        for t in phenotypes.causal_tree_indeces:
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

    def run_association(self, covariance_object, phenotypes_object, inds, index, covar, covariances_picklefile, out):
        raise ValueError("Function not defined for base class")

class TAssociationTestingRegionsMtg2(TAssociationTestingRegions):
    """
    tree-based association testing using GCTA2 (https://datadryad.org/stash/dataset/doi:10.5061/dryad.bk3j9kd8c)
    """

    def __init__(self, phenotypes, num_associations):
        super().__init__(phenotypes, num_associations)
        self.name = "regions_mtg2"

        # p-value container
        self.p_values = np.empty(self.num_associations)
        self.p_values.fill(np.nan)

    def run_association(self, covariance_object, phenotypes_object, inds, index, covar, covariances_picklefile, out):
        phenotypes_object.write_to_file_fam(inds=inds, out=out)
        


class TAssociationTestingRegionsGlimix(TAssociationTestingRegions):
    """
    tree-based association testing using glimix lmm (https://glimix-core.readthedocs.io/en/latest/lmm.html)
    """

    def __init__(self, phenotypes, num_associations):
        super().__init__(phenotypes, num_associations)
        self.name = "regions_glimix"

        # p-value container
        self.p_values = np.empty(self.num_associations)
        self.p_values.fill(np.nan)

        # other statistics
        self.V_G = np.empty(self.num_associations)
        self.V_G.fill(np.nan)
        self.V_G_over_Vp = np.empty(self.num_associations)
        self.V_G_over_Vp.fill(np.nan)

    def run_association(self, covariance_object, phenotypes_object, inds, index, covar, covariances_picklefile, out):
        """
        Run LMM to test for association between a GRM and phenotypes using limix, copied from
        https://github.com/mancusolab/sushie/blob/main/sushie/utils.py

        :param covariance_object: TCovariance, has GRM
        :param n: number of individuals
        :param covar: :math:`n \\times m` matrix for covariates.
        :param y: phenotypes
        :param index: index of genomic window
        :param out: prefix for output files

        Returns:
            :py:obj:`Tuple[float, float, float, float, float]`: A tuple of
                #. genetic variance (:py:obj:`float`) of the complex trait,
                #. :math:`h_g^2` (:py:obj:`float`) from `limix <https://github.com/limix/limix>`_ definition,
                #. :math:`h_g^2` (:py:obj:`float`) from `gcta <https://yanglab.westlake.edu.cn/software/gcta/>`_ definition,
                #. LRT test statistics (:py:obj:`float`) for :math:`h_g^2`,
                #. LRT :math:`p` value (:py:obj:`float`) for :math:`h_g^2`.
        """
        QS = economic_qs(covariance_object.covariance_matrix)
        if covar is None:
            covar = np.ones(inds.num_inds)

        # create REML object
        method = LMM(phenotypes_object.y, covar, QS, restricted=True)

        # alternative model
        method.fit(verbose=False)
        genetic_variance = method.scale * (1 - method.delta)
        random_variance = method.scale * method.delta
        v = np.var(method.mean())  # mean of the prior
        h2g_w_v = genetic_variance / (v + genetic_variance + random_variance)  # limix definition of heritability
        h2g_wo_v = genetic_variance / (genetic_variance + random_variance)  # GCTA definition of heritability
        alt_lk = method.lml()

        # null model, delta = 1 removes genetic variance component
        method.delta = 1
        method.fix("delta")
        method.fit(verbose=False)
        null_lk = method.lml()

        # likelihood ratio test
        lrt_stats = -2 * (null_lk - alt_lk)
        p_value = stats.chi2.sf(lrt_stats, 1) / 2

        # store values
        self.p_values[index] = p_value
        self.V_G[index] = genetic_variance
        self.V_G_over_Vp[index] = h2g_wo_v
        # return genetic_variance, h2g_w_v, h2g_wo_v, lrt_stats, p_value

    def write_association_results_to_file(self, window_starts, window_ends, out, phenotypes, logfile):
        table = pd.DataFrame()
        table['start'] = window_starts
        table['end'] = window_ends

        # p-values
        table['p_values'] = self.p_values

        # other stats
        table['V_G'] = self.V_G
        table['V_G_over_Vp'] = self.V_G_over_Vp

        # causal or not
        table['causal'] = np.repeat("FALSE", self.num_associations)
        table.loc[phenotypes.causal_window_indeces, 'causal'] = "TRUE"
        # table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"

        table.to_csv(out + "_trees_glimix_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + out + "_trees_glimix_results.csv'")

        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        stats.to_csv(out + "_trees_glimix_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from glimix to '" + out + "_trees_HE_stats.csv'")


class TAssociationTestingRegionsGCTA(TAssociationTestingRegions):
    """
    Run LMM to test for association between a GRM and phenotypes using GCTA
    """

    def __init__(self, phenotypes, num_associations):
        super().__init__(phenotypes, num_associations)
        self.name = "regions_GCTA"

    def run_association(self, covariance_object, phenotypes_object, inds, index, covar, covariances_picklefile, out):
        if covariance_object.write(out=out, inds=inds, covariances_picklefile=covariances_picklefile):
            self.run_association_one_window_gcta(index=index, out=out)
        else:
            print("did not run association because covariance object was not written at index", index)

    def run_association_one_window_gcta(self, index, out):
        raise ValueError("function run_association_one_window_gcta not implemented in base class")


class TAssociationTestingRegionsGCTA_HE(TAssociationTestingRegionsGCTA):
    """
    tree-based association testing using GCTA Haseman-Elston algorithm
    """

    def __init__(self, phenotypes, num_associations):

        super().__init__(phenotypes, num_associations)
        self.name = "regions_GCTA_HE"

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

    # def run_association(self, covariance_object, covariances_picklefile, inds, y=None, covar=None, n=None, index=None, out=None):
    #     covariance_object.write(out=out, inds=inds, covariances_picklefile=covariances_picklefile)
    #     self.run_association_one_window_gcta(index=index, out=out)

    def run_association_one_window_gcta(self, index, out):
        # create gcta input files, run gcta and parse output

        exit_code = subprocess.call([out + "_run_GCTA_HE.sh"])

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

        # delete GCTA results file to make sure it's not used again
        af.remove_files_with_pattern(out + '*.HEreg')

    def write_association_results_to_file(self, window_starts, window_ends, out, phenotypes, logfile):
        table = pd.DataFrame()
        table['start'] = window_starts
        table['end'] = window_ends

        # p-values
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
        table.loc[phenotypes.causal_window_indeces, 'causal'] = "TRUE"
        # table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"

        table.to_csv(out + "_trees_GCTA_HE_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + out + "_trees_GCTA_HE_results.csv'")

        stats = pd.DataFrame({'min_p_value_HECP_OLS': [np.nanmin(self.p_values_HECP_OLS)],
                              'min_p_value_HECP_Jackknife': [np.nanmin(self.p_values_HECP_Jackknife)],
                              'min_p_value_HESD_OLS': [np.nanmin(self.p_values_HESD_OLS)],
                              'min_p_value_HESD_Jackknife': [np.nanmin(self.p_values_HESD_Jackknife)],
                              'max_p_value_HECP_OLS': [np.nanmax(self.p_values_HECP_OLS)],
                              'max_p_value_HECP_Jackknife': [np.nanmax(self.p_values_HECP_Jackknife)],
                              'max_p_value_HESD_OLS': [np.nanmax(self.p_values_HESD_OLS)],
                              'max_p_value_HESD_Jackknife': [np.nanmax(self.p_values_HESD_Jackknife)]
                              })
        stats.to_csv(out + "_trees_GCTA_HE_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from HE to '" + out + "_trees_GCTA_HE_stats.csv'")


class TAssociationTestingRegionsGCTA_REML(TAssociationTestingRegionsGCTA):
    """
    tree-based association testing using CGTA REML algorithm
    """

    def __init__(self, ts_object, phenotypes):

        super().__init__(ts_object, phenotypes)
        self.name = "regions_GCTA_REML"

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

    # def run_association(self, GRM=None, y=None, covar=None, n=None, index=None, out=None):
    #     self.run_association_one_window_gcta(index=index, out=out)

    def run_association_one_window_gcta(self, index, out):
        """
        create multi_grm.txt according to https://yanglab.westlake.edu.cn/software/gcta/#GREMLinWGSorimputeddata
        @param index: index of window or tree
        @param out: prefix of output files
        @param population_structure: prefix of files containing covariance matrix used to correct for population structure
        @return:
        """
        # create gcta input files, run gcta and parse output
        exit_code = subprocess.call([out + "_run_GCTA_REML.sh"])

        # read results
        result = pd.read_table(out + "_REML.hsq")

        # replace different names in case of model with local and global GRM
        result.replace('V(G1)', 'V(G)', inplace=True)
        result.replace('Sum of V(G)/Vp', 'V(G)/Vp', inplace=True)

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

        # delete GCTA results file to make sure it's not used again
        af.remove_files_with_pattern(out + '*REML.hsq')

    def write_association_results_to_file(self, window_starts, window_ends, out, phenotypes, logfile):
        table = pd.DataFrame()
        table['start'] = window_starts
        table['end'] = window_ends
        table['p_values'] = self.p_values
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
        table.loc[phenotypes.causal_window_indeces, 'causal'] = "TRUE"

        table.to_csv(out + "_trees_GCTA_REML_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + out + "_trees_GCTA_REML_results.csv'")

        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        stats.to_csv(out + "_trees_GCTA_REML_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from tree association tests to '" + out + "_trees_GCTA_REML_stats.csv'")


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
