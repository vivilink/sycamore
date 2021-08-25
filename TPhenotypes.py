#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np


class Phenotypes:
    def __init__(self, name, ts_object):
        self.name = name
        self.N = ts_object.num_samples
        self.samp_ids = ts_object.samples()
        self.samp_variants = list(ts_object.variants(samples=self.samp_ids))
        self.num_variants = len(list(self.samp_variants))
        self.y = np.zeros(self.N)
        self.betas = [0] * self.num_variants
        self.causal_variants = []
        self.causal_betas = []
        self.causal_power = []
        self.causal_trees = []
        self.causal_tree_indeces = []
        self.filled = False
        
        print("created '" + self.name + "' phenotype object")
        
        
    def simulateEnvNoise(self, sd_environmental_noise):
        """       
        simulate random noise around zero
        
        Parameters
        ----------
        sd_environmental_noise : float
            sd of normal distribution for environmental noise.

        Returns
        -------
        None.
        """
        self.y = np.random.normal(loc=0, scale=sd_environmental_noise, size=self.N)
        self.filled = True

    def simulateFixed(self, causal_variants, betas):
        """
        Simulate phenotypes based on predefined causal variant positions and effects

        Parameters
        ----------
        causal_mutations : TYPE
            DESCRIPTION.
        betas : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        
        if(len(causal_variants) != len(betas)):
            raise ValueError("must provide equal number of causal variants and betas to simulate fixed phenotype")
        
        for v, var in enumerate(causal_variants):
            self.betas[v] = betas[v]            
            self.y[var.genotypes == 1] += betas[v]
            self.causal_variants.append(var)
            self.causal_betas.append(betas[v])
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            self.causal_power.append(betas[v]**2 * allele_freq * (1-allele_freq))
        
    def simulateUniform(self, prop_causal_mutations, sd_beta_causal_mutations, mean_beta_causal_mutation = 0):
        """
        Parameters
        ----------
        ts_object_variants : TreeSequence.variants
            variants iterator from tskit.
        prop_causal_mutations : float
            proportion of variants that should be causal.
        sd_beta_causal_mutations : TYPE
            sd of normal distribution for betas.

        Returns
        -------
        None.

        """
        #add phenotypic effect to mutations that are uniformly distributed
        for v, var in enumerate(self.samp_variants): 
            r = np.random.uniform(0,1)
            if(r < prop_causal_mutations):
                                
                #define beta
                beta = np.random.normal(loc=0, scale=sd_beta_causal_mutations, size=1)[0]
                self.betas[v] = beta
                
                #simulate phenotype
                self.y[var.genotypes == 1] += beta
                
                #save causal position
                self.causal_variants.append(var)
                self.causal_betas.append(beta)
                allele_freq = sum(var.genotypes) / len(var.genotypes)
                self.causal_power.append(beta**2 * allele_freq * (1-allele_freq))
        
        print("simulated phenotypes based on " + str(len(self.causal_variants)) + " causal variants out of a total of " + str(self.num_variants) + ".")
        self.filled = True
        
    
    def findCausalTrees(self, ts_object):
        p = 0 
        t = 0
        tree = ts_object.first()
        while p < len(self.causal_variants):    
            ## Debugging:
            ##------------
            # print("p: " + str(p))
            # print("tree index " + str(t))
            # print("causal_variants[p] + " + str(causal_variants[p]))
            # print("tree.interval.left " + str(tree.interval.left))
            # print("tree.interval.right " + str(tree.interval.right)) 
            # print("trees.at(var.site.position).get_index() " + str(trees.at(var.site.position).get_index()))
                
            if tree.interval.left <= self.causal_variants[p].site.position <= tree.interval.right:        
                #save causal tree
                self.causal_tree_indeces.append(tree.get_index())
                p += 1
                
            elif self.causal_variants[p].site.position < tree.interval.left:
                p += 1        
            
            elif self.causal_variants[p].site.position > tree.interval.right:
                tree.next()
                t += 1
        