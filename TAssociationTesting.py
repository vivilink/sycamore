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
import struct


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


class TAssociationTesting_GWAS(TAssociationTesting):
    """
    SNP based association testing
    """

    def __init__(self, phenotypes, num_typed_variants):

        super().__init__(phenotypes)

        self.num_typed_variants = num_typed_variants
        self.num_associations = self.num_typed_variants
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)

    def OLS(self, variants, inds, logfile):
        # counter respective to typed variants
        i = 0
        for v, variant in enumerate(variants.variants):
            if variants.info.iloc[v]['typed']:
                if inds.ploidy == 2:
                    genotypes = inds.get_diploid_genotypes(variant.genotypes)
                else:
                    genotypes = variant.genotypes

                # add intercept
                genotypes_test = sm.tools.add_constant(genotypes)
                PVALUE = sm.OLS(self.phenotypes.y, genotypes_test).fit().pvalues[1]
                self.p_values[i] = PVALUE
                i += 1
        logfile.info("- Ran OLS for " + str(variants.number_typed) + " variants")

    def write_to_file(self, variants, name, logfile):
        # results for each variant
        info_typed = variants.info.loc[variants.info['typed'] == True]
        info_typed['index'] = range(variants.number_typed)
        info_typed.set_index(info_typed['index'], drop=True, inplace=True)

        table = pd.DataFrame()
        table['start'] = info_typed['position']
        table['end'] = info_typed['position']
        # table['typed'] = variants.info['typed']
        table['p_value'] = self.p_values
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


class TAssociationTesting_trees(TAssociationTesting):
    """
    base class for all tree-based association tests
    """

    def __init__(self, ts_object, phenotypes):

        super().__init__(phenotypes)

        self.num_associations = ts_object.num_trees
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


class TAssociationTesting_trees_gcta(TAssociationTesting_trees):
    """
    tree-based asssociation testing using GCTA 
    """

    def __init__(self, ts_object, phenotypes):
        super().__init__(ts_object, phenotypes)

    def run_association(self, ts_object, variants, inds, out, logfile, covariance_type, skip_first_tree):
        # log progress
        start = time.time()

        for tree in ts_object.trees():
            tree_obj = tt.TTree(tree)

            self.run_association_one_tree(ts_object=ts_object, variants=variants, tree_obj=tree_obj, inds=inds, out=out,
                                          logfile=logfile, covariance_type=covariance_type,
                                          skip_first_tree=skip_first_tree)
            # log progress
            if tree.index % 100 == 0:
                end = time.time()
                logfile.info("- Ran AIM for " + str(tree.index) + " trees in " + str(round(end - start)) + " s")

        logfile.info("- Done running associations")

    def run_association_one_tree(self, ts_object, variants, tree_obj, inds, out, logfile, covariance_type,
                                 skip_first_tree):
        """
        :param ts_object: TreeSequence
        :param variants: TVariants
        :param tree_obj: TTree
        :param inds: TInds
        :param out: str
        :param logfile:
        :param covariance_type: str
        :param skip_first_tree: bool
        :return:
        """
        # logfile.info("starting association testing for tree with corrdinates: " + str(tree.interval.left) + ",
        # "  + str(tree.interval.right)) calculate covariance and write to file

        # TODO: this condition is because if you extract a region from tskit sequence, the first tree goes
        #  from zero to the first tree. This causes problems with eGRM. Needs to be investigated what the problem is
        #  and a better condition needs to be found!
        if tree_obj.height != -1 and not (skip_first_tree and tree_obj.index == 0):
            # TODO: calculating and writing should be separate functions, only write if tree is valid and matrix is
            #  not empty. Only possible to do this when eGRM functionality runs internally
            covariance = self.calculate_and_write_covariance_matrix_to_gcta_file(ts_object=ts_object, variants=variants,
                                                                                 tree_obj=tree_obj, inds=inds,
                                                                                 covariance_type=covariance_type,
                                                                                 out=out, logfile=logfile,
                                                                                 skip_first_tree=skip_first_tree)
            if covariance is not None:
                self.run_association_one_tree_gcta(tree_obj, out)

    @staticmethod
    def write_covariance_matrix_R(covariance, out):
        with open(out + '_GRM_covariance.txt', 'w') as f:
            np.savetxt(f, covariance)
        f.close()

        subprocess.call(["Rscript", os.path.dirname(sys.argv[0]) + "/create_gcta_GRM.R", out])

    @staticmethod
    def write_covariance_matrix_bin(covariance, mu, inds, out):
        """
        Write out eGRM in GCTA binary format.
        :param: covariance numpy.ndarray of expected relatedness
        :param: mu floating point number of expected mutations
        :param: numpy.ndarray/list of individual IDs
        :param: str of output
        :returns: None
        """
        # K = prefix_path.grm.bin; relatedness diagonal + lower diagonal
        # mu = prefix_path.grm.N.bin; number of shared mutations between individuals on diagonal + lower diagonal
        # samples = prefix_path.grm.id; 2 column text = family_id individual_id
        n, n = covariance.shape
        with open("{}.grm.bin".format(out), "wb") as grmfile:
            for idx in range(n):
                for jdx in range(idx + 1):
                    val = struct.pack("f", covariance[idx, jdx])
                    grmfile.write(val)

        with open("{}.grm.N.bin".format(out), "wb") as grmfile:
            for idx in range(n):
                for jdx in range(idx + 1):
                    val = struct.pack("f", mu)
                    grmfile.write(val)

        with open("{}.grm.id".format(out), "w") as grmfile:
            for idx in range(n):
                fid = 0
                iid = inds.names[idx]
                grmfile.write("\t".join([str(fid), str(iid)]) + os.linesep)

    # TODO: could be static?
    def calculate_and_write_covariance_matrix_to_gcta_file(self, ts_object, variants, tree_obj, inds, covariance_type,
                                                           out, skip_first_tree, logfile):
        """
        Writes covariance and other files necessary to run gcta. The program egrm does that automatically, the scaled

        Parameters
        ----------
        ts_object : tskit.treeSequence
        variants : TVariants
        tree_obj : TTree.py
        inds : TInds
        covariance_type : str
        out : str
        skip_first_tree: bool
        logfile : IndentedLoggerAdapter

        Raises
        ------
        ValueError
            If the covariance_type is not a recognized method.

        Returns
        -------
        Covariance: ndarray(inds.num_inds, inds.num_inds).

        """
        if covariance_type == "scaled":
            covariance = tree_obj.get_covariance_scaled(inds=inds)
            self.write_covariance_matrix_R(covariance=covariance, out=out)

        elif covariance_type == "eGRM":
            # trees = ts_object.keep_intervals(np.array([[tree_obj.start, tree_obj.end]]), simplify=True)
            covariance, mu = tree_obj.get_eGRM(tskit_obj=ts_object, tree_obj=tree_obj, inds=inds, out=out, logfile=logfile,
                                               skip_first_tree=skip_first_tree)
            self.write_covariance_matrix_bin(covariance=covariance, mu=mu, inds=inds, out=out)

            # if np.trace(covariance) != inds.num_inds:
            # raise ValueError("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
            # inds.num_inds) + " but obtained " + str(np.trace(covariance)))
            # logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
            #     inds.num_inds) + " but obtained " + str(np.trace(covariance)))

        elif covariance_type == "GRM":
            if inds.ploidy == 2:
                raise ValueError("GRM not implemented for diploids")
            covariance, mu = tree_obj.get_GRM(variants=variants, inds=inds, out=out, logfile=logfile)
            if covariance is None:
                return None
            if np.trace(covariance) != inds.num_inds:
                logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
                    inds.num_inds) + " but obtained " + str(np.trace(covariance)))
            self.write_covariance_matrix_bin(covariance=covariance, mu=mu, inds=inds, out=out)

        else:
            raise ValueError("Did not recognize " + str(covariance_type) + " as a covariance type")

        return covariance

    def run_association_one_tree_gcta(self, tree, out):
        raise ValueError("function run_association_one_tree_gcta not implemented in base class")


class TAssociationTesting_trees_gcta_HE(TAssociationTesting_trees_gcta):
    """
    tree-based association testing using GCTA Haseman-Elston algorithm
    """

    def __init__(self, ts_object, phenotypes):

        super().__init__(ts_object, phenotypes)

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

    def run_association_one_tree_gcta(self, tree, out):
        # create gcta input files, run gcta and parse output
        exit_code = subprocess.call([os.path.dirname(sys.argv[0]) + "/run_gcta_HE.sh", out])
        # exit_code = subprocess.call([os.getcwd() + "/run_gcta_HE.sh", out])

        # read results
        HE_CP = pd.read_table(out + "_HE-CP_result.txt")
        HE_SD = pd.read_table(out + "_HE-SD_result.txt")

        # p-values
        self.p_values_HECP_OLS[tree.index] = HE_CP["P_OLS"][1]
        if HE_CP["P_OLS"][1] < 0:
            raise ValueError("tree index", tree.index, "produced negative p-value for CP OLS")

        self.p_values_HECP_Jackknife[tree.index] = HE_CP["P_Jackknife"][1]
        if HE_CP["P_Jackknife"][1] < 0:
            raise ValueError("tree index", tree.index, "produced negative p-value for CP Jackknife")

        self.p_values_HESD_OLS[tree.index] = HE_SD["P_OLS"][1]
        if HE_SD["P_OLS"][1] < 0:
            raise ValueError("tree index", tree.index, "produced negative p-value for SD OLS")

        self.p_values_HESD_Jackknife[tree.index] = HE_SD["P_Jackknife"][1]
        if HE_SD["P_Jackknife"][1] < 0:
            raise ValueError("tree index", tree.index, "produced negative p-value for SD Jackknife")

        # other statistics
        self.V_G_over_Vp_HECP[tree.index] = HE_CP["Estimate"][1]
        self.V_G_over_Vp_HESD[tree.index] = HE_SD["Estimate"][1]

        self.V_G_over_Vp_SE_OLS_HECP[tree.index] = HE_CP["SE_OLS"][1]
        self.V_G_over_Vp_SE_OLS_HESD[tree.index] = HE_SD["SE_OLS"][1]
        self.V_G_over_Vp_SE_Jackknife_HECP[tree.index] = HE_CP["SE_Jackknife"][1]
        self.V_G_over_Vp_SE_Jackknife_HESD[tree.index] = HE_SD["SE_Jackknife"][1]

    def write_to_file(self, ts_object, out, logfile):
        table = pd.DataFrame()
        table['start'] = ts_object.breakpoints(as_array=True)[
                         0:self.num_associations]  # otherwise the next start is included, i think this tree is removed due to incompleteness when taking tree subset
        table['end'] = table['start']

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
        table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"

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


class TAssociationTesting_trees_gcta_REML(TAssociationTesting_trees_gcta):
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

    def run_association_one_tree_gcta(self, tree, out):
        # create gcta input files, run gcta and parse output
        exit_code = subprocess.call([os.path.dirname(sys.argv[0]) + "/run_gcta_REML.sh", out])

        # read results
        result = pd.read_table(out + "_REML.hsq")
        result_pvalue = float(result['Variance'][result['Source'] == 'Pval'])
        if result_pvalue < 0:
            raise ValueError("Negative p-value for tree starting at " + str(tree.interval.left))
        if result_pvalue > 1:
            raise ValueError("p-value larger than 1 for tree starting at " + str(tree.interval.left))

        self.p_values[tree.index] = result_pvalue
        if (result_pvalue < 0):
            raise ValueError("tree index", tree.index, "produced negative p-value with REML")

        self.V_G[tree.index] = float(result['Variance'][result['Source'] == 'V(G)'])
        self.V_e[tree.index] = float(result['Variance'][result['Source'] == 'V(e)'])
        self.Vp[tree.index] = float(result['Variance'][result['Source'] == 'Vp'])
        self.V_G_over_Vp[tree.index] = float(result['Variance'][result['Source'] == 'V(G)/Vp'])
        self.logL[tree.index] = float(result['Variance'][result['Source'] == 'logL'])
        self.logL0[tree.index] = float(result['Variance'][result['Source'] == 'logL0'])
        self.LRT[tree.index] = float(result['Variance'][result['Source'] == 'LRT'])

        self.V_G_SE[tree.index] = float(result['SE'][result['Source'] == 'V(G)'])
        self.V_e_SE[tree.index] = float(result['SE'][result['Source'] == 'V(e)'])
        self.Vp_SE[tree.index] = float(result['SE'][result['Source'] == 'Vp'])
        self.V_G_over_Vp_SE[tree.index] = float(result['SE'][result['Source'] == 'V(G)/Vp'])

    def write_to_file(self, ts_object, name, logfile):
        table = pd.DataFrame()
        table['start'] = ts_object.breakpoints(as_array=True)[
                         0:self.num_associations]  # otherwise the next start is included, i think this tree is removed due to incompleteness when taking tree subset
        table['end'] = ts_object.breakpoints(as_array=True)[1:]
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
        table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"

        table.to_csv(name + "_trees_REML_results.csv", index=False, header=True)
        logfile.info("- Wrote results from tree association tests to '" + name + "_trees_REML_results.csv'")

        stats = pd.DataFrame({'min_p_value': [np.nanmin(self.p_values)],
                              'max_p_value': [np.nanmax(self.p_values)]
                              })
        stats.to_csv(name + "_trees_REML_stats.csv", index=False, header=True)
        logfile.info("- Wrote stats from tree association tests to '" + name + "_trees_REML_stats.csv'")


class TTreeAssociation_Mantel(TAssociationTesting_trees):

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
            tree_obj = tt.TTree(tree, inds.num_haplotypes)
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
