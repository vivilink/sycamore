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
import subprocess
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt


class TGWAS:
    
    def __init__(self, phenotypes):
        self.name = phenotypes.name
        self.phenotypes = phenotypes
        self.num_associations = -1
        self.p_values = np.empty(0)
        # self.p_values_init = False
        
    # def _check_compatibility(self, ts_object, phenotypes):
    #     #check if ts_object is compatible with phenotype
    #     if phenotypes.num_variants != self.num_associations:
    #         raise ValueError("Phenotype object must contain same number of sampled variants as ts_object to create a GWAS object. Phenotype: " + str(phenotypes.num_variants) + " and ts_object: " + str(self.num_associations))
    #     if phenotypes.N != ts_object.num_samples:
    #         raise ValueError("Phenotype object must contain same number of samples as ts_object to create a GWAS object. Phenotype: " + str(phenotypes.N) + " and ts_object: " + str(ts_object.num_samples))

        
    def p_value_dist(self, subplot, num_bins):
        subplot.hist(self.p_values, bins=num_bins)
        subplot.set(xlabel='p-values', ylabel='density', title=self.phenotypes.name)
        
    def chiSquared(self, num_bins):
        h,bin_edges = np.histogram(self.p_values, num_bins)
        tmp = scipy.stats.chisquare(f_obs=h, f_exp=((len(self.p_values))/(num_bins)))
        return(tmp)

class TpGWAS(TGWAS):
    """
    SNP based association testing
    """
    
    def __init__(self, phenotypes, num_typed_variants):
        
        super().__init__(phenotypes)
        
        self.num_typed_variants = num_typed_variants
        self.num_associations = self.num_typed_variants
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        # self.q_values = np.empty(self.num_associations)
        
    def OLS(self, variants, inds, logfile):
        # i = 0
        for v, variant in enumerate(variants.variants):
            if variants.info.iloc[v]['typed'] == True:
                if inds.ploidy == 2:
                    genotypes = inds.get_diploid_genotypes(variant.genotypes)      
                else: 
                    genotypes = variant.genotypes
                
                #add intercept
                genotypes_test = sm.tools.add_constant(genotypes)
                
                # print("self.phenotypes.y", self.phenotypes.y)
                PVALUE = sm.OLS(self.phenotypes.y, genotypes_test).fit().pvalues[1] 
                # print("PVALUE", PVALUE)
                self.p_values[v] = PVALUE
                # i += 1
        logfile.info("- Ran OLS for all " + str(variants.number_typed) + " variants of " + self.name)
        

        
    def writeToFile(self, variants, name, logfile):        
        #results for each variant
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
        table.to_csv(name + "_variants_results.csv", index = False, header = True)        
        
        #summary statistics
        stats = pd.DataFrame()
        stats['min_p_value'] = min(self.p_values)
        stats['max_p_value'] = max(self.p_values) 
        logfile.info("- Writing stats from OLS to '" + name + "_variants_stats.csv'")           
        stats.to_csv(name + "_variants_stats.csv", index = False, header = True)        
        
            
    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, logfile, size=1, n_snps_lowess = 0, *args):
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
        p_values[(np.where(p_values == 0.0))] = np.nextafter(0,1)
        q_values = -np.log10(p_values)

        subplot.scatter(variant_positions[index_min:index_max], q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='variant position', ylabel='q-value', title=self.phenotypes.name)
        for v, var in enumerate(self.phenotypes.causal_variants):
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            colscale = self.phenotypes.causal_power[v] 
            subplot.axvline(x=var.site.position, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=var.site.position, color="black", lw=0.5)
        
        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess/len(q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=q_values[index_min:index_max], exog=variant_positions[index_min:index_max], frac=fraction, return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")
        
        subplot.axhline(y=8, color="red", lw=0.5)
        
    def manhattan_plot(self, variant_positions, logfile, plots_dir, *args):
        fig, ax = plt.subplots(1,figsize=(10,10))
        self.manhattan_plot_subset(variant_positions, ax, 0, len(self.p_values), logfile, *args)
        fig.tight_layout()
        fig.set_size_inches(30, 30)
        fig.savefig(plots_dir + 'OLS_GWAS.png', bbox_inches='tight')    
    
class TtGWAS(TGWAS):
    """
    base class for all tree-based association tests
    """
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(phenotypes)

        self.num_associations = ts_object.num_trees
        # self._check_compatibility(ts_object, phenotypes)
               
        
    def manhattan_plot(self, variant_positions, subplot, logfile, *args):
        self.manhattan_plot_subset(variant_positions, subplot, 0, self.num_associations, p_values = self.p_values, logfile = logfile)    
    
    def manhattan_plot_special_pvalues(self, variant_positions, p_values, subplot, logfile, title_supplement = "", *args):
        logfile.info("Plotting " + str(self.num_associations) + " associations")
        self.manhattan_plot_subset(variant_positions, subplot, 0, self.num_associations, p_values = p_values, logfile = logfile, title_supplement = title_supplement)    
        
    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, p_values, logfile, title_supplement = "", size=1, n_snps_lowess = 0, *args):
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
        
        p_values[(np.where(p_values == 0.0))] = np.nextafter(0,1)
        q_values = -np.log10(p_values)

        subplot.scatter(variant_positions[index_min:index_max], q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='tree index', ylabel='q-value', title=self.phenotypes.name + title_supplement)
        for t in self.phenotypes.causal_tree_indeces:
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            # colscale = self.phenotypes.causal_power[v] 
            # subplot.axvline(x=t, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=t, color="black", lw=0.5)
        
        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess/len(q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=q_values[index_min:index_max], exog=variant_positions[index_min:index_max], frac=fraction, return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")
            
        subplot.axhline(y=8, color="red", lw=0.5)


class HE_tGWAS(TtGWAS):
    """
    tree-based asssociation testing using GCTA Haseman-Elston algorithm
    """
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(ts_object, phenotypes)
        
        #p-value containers
        self.p_values_HECP_OLS = np.empty(self.num_associations)
        self.p_values_HECP_Jackknife = np.empty(self.num_associations)
        self.p_values_HESD_OLS = np.empty(self.num_associations)
        self.p_values_HESD_Jackknife = np.empty(self.num_associations)
        
        #other statistics
        self.V_G_over_Vp_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_HESD = np.empty(self.num_associations)
        
        self.V_G_over_Vp_SE_OLS_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_OLS_HESD = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_Jackknife_HECP = np.empty(self.num_associations)
        self.V_G_over_Vp_SE_Jackknife_HESD = np.empty(self.num_associations)

    def run_association(self, ts_object, inds, out, logfile, covariance_scaled):        
        self.phenotypes.write_to_file_gcta(out, logfile)        

        #log progress
        start = time.time()
        
        for tree in ts_object.trees():
            self.run_association_one_tree(tree, inds, out, logfile, covariance_scaled)  
            #log progress
            if tree.index % 100 == 0:
                end = time.time()
                logfile.info("- Ran HE for " + str(tree.index) + " trees in " + str(round(end-start)) + " s")

            
        logfile.info("- done running associations")
        
    def run_association_one_tree(self, tree, inds, out, logfile, covariance_scaled):      
        
        # logfile.info("- running association test on tree with interval: " + str(tree.interval.left) + "," + str(tree.interval.right))

        #calculate covariance and write to file
        tree_obj = tt.TTree(tree, inds.num_haplotypes)
        if covariance_scaled == True:
            covariance = tree_obj.covariance_scaled(inds, logfile)
        else:
            covariance = tree_obj.covariance(inds)     
            
        if covariance is not None:
            
            with open(out + "_GRM_covariance.txt", 'w') as f:
                np.savetxt(f, covariance)
            f.close()
            
    
            # create gcta input files, run gcta and parse output
            exit_code = subprocess.call([os.path.dirname(sys.argv[0]) + "/run_gcta_HE.sh", out])
            # exit_code = subprocess.call([os.getcwd() + "/run_gcta_HE.sh", out])
    
            # read results
            HE_CP = pd.read_table(out + "_HE-CP_result.txt")
            HE_SD = pd.read_table(out + "_HE-SD_result.txt")
            
            #p-values        
            self.p_values_HECP_OLS[tree.index] = HE_CP["P_OLS"][1]
            if(HE_CP["P_OLS"][1] < 0):
                raise ValueError("tree index", tree.index, "produced negative p-value for CP OLS")
                
            self.p_values_HECP_Jackknife[tree.index] = HE_CP["P_Jackknife"][1]        
            if(HE_CP["P_Jackknife"][1] < 0):
                raise ValueError("tree index", tree.index, "produced negative p-value for CP Jackknife")
    
            self.p_values_HESD_OLS[tree.index] = HE_SD["P_OLS"][1]
            if(HE_SD["P_OLS"][1] < 0):
                raise ValueError("tree index", tree.index, "produced negative p-value for SD OLS")
    
            self.p_values_HESD_Jackknife[tree.index] = HE_SD["P_Jackknife"][1]
            if(HE_SD["P_Jackknife"][1] < 0):
                raise ValueError("tree index", tree.index, "produced negative p-value for SD Jackknife")
                
            #other statistics
            self.V_G_over_Vp_HECP[tree.index] = HE_CP["Estimate"][1]
            self.V_G_over_Vp_HESD[tree.index] = HE_SD["Estimate"][1]
    
            self.V_G_over_Vp_SE_OLS_HECP[tree.index] = HE_CP["SE_OLS"][1]
            self.V_G_over_Vp_SE_OLS_HESD[tree.index] = HE_SD["SE_OLS"][1]
            self.V_G_over_Vp_SE_Jackknife_HECP[tree.index] = HE_CP["SE_Jackknife"][1]
            self.V_G_over_Vp_SE_Jackknife_HESD[tree.index] = HE_SD["SE_Jackknife"][1]
        

    def write_to_file(self, ts_object, name, logfile):
        table = pd.DataFrame()
        table['start'] = ts_object.breakpoints(as_array=True)[0:self.num_associations] #otherwise the next start is included, i think this tree is removed due to incompleteness when taking tree subset
        table['end'] = table['start']
        
        #p-values
        table['p_values_HECP_OLS'] = self.p_values_HECP_OLS
        table['p_values_HECP_Jackknife'] = self.p_values_HECP_Jackknife
        table['p_values_HESD_OLS'] = self.p_values_HESD_OLS
        table['p_values_HESD_Jackknife'] = self.p_values_HESD_Jackknife     

        #other stats
        table['V_G_over_Vp_HECP'] = self.V_G_over_Vp_HECP
        table['V_G_over_Vp_HESD'] = self.V_G_over_Vp_HESD
        table['V_G_over_Vp_SE_OLS_HECP'] = self.V_G_over_Vp_SE_OLS_HECP
        table['V_G_over_Vp_SE_OLS_HESD'] = self.V_G_over_Vp_SE_OLS_HESD
        table['V_G_over_Vp_SE_Jackknife_HECP'] = self.V_G_over_Vp_SE_Jackknife_HECP
        table['V_G_over_Vp_SE_Jackknife_HESD'] = self.V_G_over_Vp_SE_Jackknife_HESD

        #causal or not
        table['causal'] = np.repeat("FALSE", self.num_associations)
        table.loc[self.phenotypes.causal_tree_indeces, 'causal'] = "TRUE"
        
        table.to_csv(name + "_trees_HE_results.csv", index = False, header = True)        
        logfile.info("- Wrote results from tree association tests to '" + name + "_trees_HE_results.csv'")
        
        stats = pd.DataFrame()
        stats['min_p_value_HECP_OLS'] = min(self.p_values_HECP_OLS)
        stats['min_p_value_HECP_Jackknife'] = min(self.p_values_HECP_Jackknife)
        stats['min_p_value_HESD_OLS'] = min(self.p_values_HESD_OLS)
        stats['min_p_value_HESD_Jackknife'] = min(self.p_values_HESD_Jackknife)
        
        stats['max_p_value_HECP_OLS'] = max(self.p_values_HECP_OLS)
        stats['max_p_value_HECP_Jackknife'] = max(self.p_values_HECP_Jackknife)
        stats['max_p_value_HESD_OLS'] = max(self.p_values_HESD_OLS)
        stats['max_p_value_HESD_Jackknife'] = max(self.p_values_HESD_Jackknife)
 
        stats.to_csv(name + "_trees_HE_stats.csv", index = False, header = True)        
        logfile.info("- Wrote stats from HE to '" + name + "_trees_HE_stats.csv'")   


class REML_tGWAS(TtGWAS):
    """
    tree-based association testing using CGTA REML algorithm
    """
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(ts_object, phenotypes)
        
        # results containers
        self.p_values = np.empty(self.num_associations)
        self.V_G = np.empty(self.num_associations)
        self.V_e = np.empty(self.num_associations)
        self.Vp = np.empty(self.num_associations)
        self.V_G_over_Vp = np.empty(self.num_associations)
        self.logL = np.empty(self.num_associations)
        self.logL0 = np.empty(self.num_associations)
        self.LRT = np.empty(self.num_associations)
        
        self.V_G_SE = np.empty(self.num_associations)
        self.V_e_SE = np.empty(self.num_associations)
        self.Vp_SE = np.empty(self.num_associations)
        self.V_G_over_Vp_SE = np.empty(self.num_associations)


    def run_association_one_tree(self, tree, inds, out, logfile, covariance_scaled):  
        # logfile.info("starting association testing for tree with corrdinates: " + str(tree.interval.left) + ","  + str(tree.interval.right))

        #calculate covariance and write to file
        tree_obj = tt.TTree(tree, inds.num_haplotypes)
        if covariance_scaled == True:
            covariance = tree_obj.covariance_scaled(inds, logfile)
        else:
            covariance = tree_obj.covariance(inds)
            
        if covariance is not None:
        
            with open(out + '_GRM_covariance.txt', 'w') as f:
                np.savetxt(f, covariance)
            f.close()
                    
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
            if(result_pvalue < 0):
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


    def run_association(self, ts_object, inds, out, logfile, covariance_scaled):        
        self.phenotypes.write_to_file_gcta(out, logfile)        
        
        start = time.time()        
        for tree in ts_object.trees():
            self.run_association_one_tree(tree, inds, out, logfile, covariance_scaled)  
            # if tree.index == 0:
            #     raise ValueError("stop!!!")
            if tree.index % 100 == 0:
                end = time.time()
                logfile.info("- Ran REML for " + str(tree.index) + " trees in " + str(round(end-start)) + " s")


    def write_to_file(self, ts_object, name, logfile):
        table = pd.DataFrame()
        table['start'] = ts_object.breakpoints(as_array=True)[0:self.num_associations] #otherwise the next start is included, i think this tree is removed due to incompleteness when taking tree subset
        table['end'] = table['start']
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
        
        table.to_csv(name + "_trees_REML_results.csv", index = False, header = True)        
        logfile.info("- Wrote results from tree association tests to '" + name + "_trees_REML_results.csv'")
        
        stats = pd.DataFrame()
        stats['min_p_value'] = min(self.p_values)
        stats['max_p_value'] = max(self.p_values) 
        stats.to_csv(name + "_trees_stats.csv", index = False, header = True)        
        logfile.info("- Wrote stats from tree association tests to '" + name + "_trees_REML_stats.csv'")   


class Mantel_tGWAS(TtGWAS):
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(ts_object, phenotypes)
        
        # p-value containers
        self.p_values = np.empty(self.num_associations)

    def run_Mantel(self, ts_object, phenotypes, inds):
        #test for associations
        diffs = phenotypes.diffs()
        start = time.time()
        for tree in ts_object.trees():
            if tree.index % 100 == 0:
                end = time.time()
                print("Ran Mantel for", tree.index, "trees in ", round(end-start), "s")
            if tree.total_branch_length == 0: 
                print("tree's total branch length is zero")
                continue
            tree_obj = tt.TTree(tree, inds.num_haplotypes)
            tmrca = tree_obj.TMRCA(inds.num_haplotypes)
            # print("tmrca",tmrca)
            self.p_values[tree.index] = ut.mantel(tmrca, diffs)
            if(self.p_values[tree.index] < 0):
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
    






#-----------------   
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


       