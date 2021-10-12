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
from limix_lmm.lmm_core import LMMCore
import utils as ut
import time
import pandas as pd
import subprocess



class TGWAS:
    
    def __init__(self, phenotypes):
        # self.ts_object = ts_object
        self.name = phenotypes.name
        self.phenotypes = phenotypes
        self.num_associations = -1
        self.p_values = np.empty(0)
        self.q_values = np.empty(0)
        self.num_variants = phenotypes.num_variants
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
    
    def __init__(self, phenotypes):
        
        super().__init__(phenotypes)
        
        self.num_associations = self.num_variants
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        self.q_values = np.empty(self.num_associations)
        
    def OLS(self, variants):
        for v, variant in enumerate(variants.variants):
            self.p_values[v] = sm.OLS(self.phenotypes.y, sm.tools.add_constant(variant.genotypes)).fit().pvalues[1]
        self.q_values = -np.log10(self.p_values)
        print("ran OLS for all variants of " + self.name)
            
    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, size=1, n_snps_lowess = 0, *args):
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
        
        subplot.scatter(variant_positions[index_min:index_max], self.q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='variant position', ylabel='q-value', title=self.phenotypes.name)
        for v, var in enumerate(self.phenotypes.causal_variants):
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            colscale = self.phenotypes.causal_power[v] 
            subplot.axvline(x=var.site.position, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=var.site.position, color="black", lw=0.5)
        
        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess/len(self.q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=self.q_values[index_min:index_max], exog=variant_positions[index_min:index_max], frac=fraction, return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")
        
        subplot.axhline(y=8, color="red", lw=0.5)

        
    def manhattan_plot(self, variant_positions, subplot, *args):
        self.manhattan_plot_subset(variant_positions, subplot, 0, len(self.p_values), *args)
    
    
class TtGWAS(TGWAS):
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(phenotypes)

        self.num_associations = ts_object.num_trees
        # self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        self.q_values = np.empty(self.num_associations)
        
        self.p_values_HECP_OLS = np.empty(self.num_associations)
        self.p_values_HECP_Jackknife = np.empty(self.num_associations)
        self.p_values_HESD_OLS = np.empty(self.num_associations)
        self.p_values_HESD_Jackknife = np.empty(self.num_associations)
        
        self.lrt = np.empty(self.num_associations)

    def runMantel(self, ts_object, phenotypes, N):
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
            tree_obj = tt.TTree(tree, N)
            tmrca = tree_obj.TMRCA(N)
            # print("tmrca",tmrca)
            self.p_values[tree.index] = ut.mantel(tmrca, diffs)
            if(self.p_values[tree.index] < 0):
                print(tmrca)
                raise ValueError("p-value is negative")
            
    def runLimix(self, ts_object, N, y, F, random):   
        G = np.zeros(N).reshape(N,1) 
        # G = np.random.binomial(1, 0.5, N).reshape(N,1)
        # Inter = np.zeros(N).reshape(N,1)
        Inter = None

        for tree in ts_object.trees():
            # if tree.index == ts_object.num_trees-1:
            #     continue
            # if tree.index % 1000 == 0:
            print("tree index: ",tree.index)
            tree_obj = tt.TTree(tree, N)
            lmm = LMMCore(y, F, tree_obj.solving_function)    
            lmm.process(G, Inter) #this needs step param to produce only one p-value per tree. for this i need the number of sites per tree, or just use 1?
            self.p_values[tree.index] = lmm.getPv()
            print("p-value", self.p_values[tree.index])
            # raise ValueError("printing covariance")
            # beta = lmm.getBetaSNP()
            # beta_ste = lmm.getBetaSNPste()
            # self.lrt[tree.index] = lmm.getLRT() #likelihood ratio
            
    def runCGTA_HE(self, ts_object, N):
        
        # mimic example phenotype file from gcta (first column=family, second=ind id, third=pheno value)
        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.arange(1,N+1)
        tmp_pheno['2'] = tmp_pheno['1']
        tmp_pheno['3'] = self.phenotypes.y        
        tmp_pheno.to_csv("phenotypes.phen", sep=' ', index=False, header=False)
        
        # test each tree with gcta
        start = time.time()
        for tree in ts_object.trees():
            if tree.index % 100 == 0:
                end = time.time()
                print("Ran HE for", tree.index, "trees in ", round(end-start), "s")
            tree_obj = tt.TTree(tree, N)
            covariance = tree_obj.covariance(N)
            # print("covariance is symmetric",(covariance == covariance.T).all())
            # print(covariance[1,10])
            # print(covariance[10,1])

            with open('GRM_covariance.txt', 'w') as f:
                np.savetxt(f, covariance)
            f.close()
            
            # run gcta and parse output
            # TODO: locations of gcta and run_gcta_HE scripts only exist for me
            exit_code = subprocess.call('/data/ARGWAS/argwas/run_gcta_HE.sh')
            
            # read results
            HE_CP = pd.read_table("HE-CP_result.txt")
            HE_SD = pd.read_table("HE-SD_result.txt")

            self.p_values_HECP_OLS[tree.index] = HE_CP["P_OLS"][1]
            self.p_values_HECP_Jackknife[tree.index] = HE_CP["P_Jackknife"][1]
            self.p_values_HESD_OLS[tree.index] = HE_SD["P_OLS"][1]
            self.p_values_HESD_Jackknife[tree.index] = HE_SD["P_Jackknife"][1]

            
            # raise ValueError("stop after writing file")
            # subprocess.run(["/data/ARGWAS/gcta_1.93.2beta/gcta64 --HEreg --grm GRM_covariance --pheno phenotypes.phen --out GRM_covariance_tests"])
            
            
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

    def manhattan_plot(self, variant_positions, subplot, *args):
        self.manhattan_plot_subset(variant_positions, subplot, 0, self.num_associations, p_values = self.p_values)    

    def manhattan_plot_special_pvalues(self, variant_positions, p_values, subplot, *args):
        print("num associations",self.num_associations)
        self.manhattan_plot_subset(variant_positions, subplot, 0, self.num_associations, p_values = p_values)    
        
    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, p_values, size=1, n_snps_lowess = 0, *args):
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
        q_values = -np.log10(p_values)
        
        subplot.scatter(variant_positions[index_min:index_max], q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='tree index', ylabel='q-value', title=self.phenotypes.name)
        for t in self.phenotypes.causal_tree_indeces:
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            # colscale = self.phenotypes.causal_power[v] 
            # subplot.axvline(x=t, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=t, color="black", lw=0.5)
        
        if n_snps_lowess > 0:
            fraction = min(n_snps_lowess/len(q_values[index_min:index_max]), 1)
            low = sm.nonparametric.lowess(endog=q_values[index_min:index_max], exog=variant_positions[index_min:index_max], frac=fraction, return_sorted=False)
            subplot.plot(variant_positions[index_min:index_max], low, color="red")
