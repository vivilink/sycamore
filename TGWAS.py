#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:37:12 2021

@author: linkv
"""
import numpy as np
import utils as ut
import statsmodels.api as sm
import matplotlib.pyplot as plt
import scipy

class TGWAS:
    
    def __init__(self, phenotypes):
        # self.ts_object = ts_object
        self.name = phenotypes.name
        self.phenotypes = phenotypes
        self.num_associations = -1
        self.p_values = np.empty(0)
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
            self.p_values[v] = sm.OLS(variant.genotypes, self.phenotypes.y).fit().pvalues[0]
        self.q_values = -np.log10(self.p_values)
        print("ran OLS for all variants of " + self.name)
            
    def manhattan_plot_subset(self, variant_positions, subplot, index_min, index_max, size=1, n_snps_lowess = 100,*args):
        if index_min < 0 or index_max < 0:
            raise ValueError("data subset indeces for manhattan plotting must be positive")
        if index_max > len(self.p_values) or index_min > len(self.p_values):
            raise ValueError("data subset index cannot be larger than number of p-values")            
        
        subplot.scatter(variant_positions[index_min:index_max], self.q_values[index_min:index_max], s=size, *args)
        subplot.set(xlabel='position', ylabel='q-value', title=self.phenotypes.name)
        for v, var in enumerate(self.phenotypes.causal_variants):
            # print("power " + str(self.phenotypes.causal_power[v]) + " pos " + str(var.site.position))
            colscale = self.phenotypes.causal_power[v] 
            subplot.axvline(x=var.site.position, color=str(0.3), alpha = 0.5, lw=colscale*100)
            subplot.axvline(x=var.site.position, color="black", lw=0.5)
        
        fraction = min(n_snps_lowess/len(self.q_values[index_min:index_max]), 1)
        low = sm.nonparametric.lowess(endog=self.q_values[index_min:index_max], exog=variant_positions[index_min:index_max], frac=fraction, return_sorted=False)
        subplot.plot(variant_positions[index_min:index_max], low, color="red")
        
    def manhattan_plot(self, variant_positions, subplot, *args):
        self.manhattan_plot_subset(variant_positions, subplot, 0, len(self.p_values), *args)
    
    
class TtGWAS(TGWAS):
    
    def __init__(self, ts_object, phenotypes):
        
        super().__init__(ts_object, phenotypes)

        self.num_associations = ts_object.num_trees
        self._check_compatibility(ts_object, phenotypes)
        self.p_values = np.empty(self.num_associations)
        self.q_values = np.empty(self.num_associations)
        
    def mantel(self):
        diffs = ut.diff(self.phenotypes.y, self.phenotypes.N)
        for t,tree in enumerate(self.ts_object.trees()):
            if tree.total_branch_length == 0: 
                continue
            tmrca = ut.TMRCA(tree, self.phenotypes.N)
            self.p_values[t] = ut.mantel(ut.make(tmrca), diffs)
