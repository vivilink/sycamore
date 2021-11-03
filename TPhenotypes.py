#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np
import pandas as pd


class Phenotypes:
    def __init__(self, name, variants, num_inds, logfile):
        self.name = name
        self.num_variants = variants.number
        self.N = num_inds
        self.y = np.zeros(self.N)
        self.betas = [0] * self.num_variants
        self.causal_variants = []
        self.causal_betas = []
        self.causal_power = []
        self.causal_trees = []
        self.causal_variant_indeces = []
        self.causal_tree_indeces = []
        self.filled = False
        
        logfile.info("created '" + self.name + "' phenotype object")
        
    def returnRandomState(self, random):
        print(random.random.get_state()[1][0])
        print(random.random.uniform(0,1,1))
       
        
    def simulateEnvNoise(self, sd_environmental_noise, random):
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
        self.y = random.random.normal(loc=0, scale=sd_environmental_noise, size=self.N)
        self.filled = True

    def simulateFixed(self, variants, causal_variant_indeces, betas, logfile):
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
        causal_variants = [variants.info['variant'][i] for i in causal_variant_indeces]
        if(len(causal_variants) != len(betas)):
            raise ValueError("must provide equal number of causal variants and betas to simulate fixed phenotype")
            
              
        for v, var in enumerate(causal_variants):
            self.betas[v] = betas[v]            
            self.y[var.genotypes == 1] += betas[v]
            self.causal_variants.append(var)
            self.causal_betas.append(betas[v])
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            self.causal_power.append(betas[v]**2 * allele_freq * (1-allele_freq))            
            self.causal_variant_indeces.append(causal_variant_indeces[v])
            
            logfile.info("Simulated causal variant at index " + str(causal_variant_indeces[v]) + " with beta " + str(betas[v]) + " and allele freq " + str(allele_freq) + " resulting in a power of " + str(betas[v]**2 * allele_freq * (1-allele_freq)))
        
    def simulateUniform(self, variants, prop_causal_mutations, sd_beta_causal_mutations, random, mean_beta_causal_mutation = 0):
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
        for v, var in enumerate(variants.info['variant']): 
            r = random.random.uniform(0,1,1)

            if(r < prop_causal_mutations):
                                
                #define beta
                beta = random.random.normal(loc=0, scale=sd_beta_causal_mutations, size=1)[0]
                self.betas[v] = beta
                
                #simulate phenotype
                self.y[var.genotypes == 1] += beta
                self.y[var.genotypes == 2] += 2 * beta

                #save causal position
                self.causal_variants.append(var)
                self.causal_betas.append(beta)
                allele_freq = variants.info['allele_freq'][v]
                self.causal_power.append(beta**2 * allele_freq * (1-allele_freq))
                self.causal_variant_indeces.append(v)
        
        print("simulated phenotypes based on " + str(len(self.causal_variants)) + " causal variants out of a total of " + str(self.num_variants) + ".")
        self.filled = True
        
    
    def findCausalTrees(self, ts_object):
        for v in self.causal_variants:
            causal_tree = ts_object.at(v.site.position)
            self.causal_tree_indeces.append(causal_tree.get_index())
        
        # p = 0 
        # t = 0
        # tree = ts_object.first()
        # while p < len(self.causal_variants):    
        #     # # Debugging:
        #     # #------------
        #     # print("p: " + str(p))
        #     # print("tree index " + str(t))
        #     # print("causal_variants[p] + " + str(self.causal_variants[p]))
        #     # print("tree.interval.left " + str(tree.interval.left))
        #     # print("tree.interval.right " + str(tree.interval.right)) 
        #     # print("trees.at(var.site.position).get_index() " + str(ts_object.at(self.causal_variants[p].site.position).get_index()))
                
        #     if tree.interval.left <= self.causal_variants[p].site.position <= tree.interval.right:        
        #         #save causal tree
        #         self.causal_tree_indeces.append(tree.get_index())
        #         p += 1
                
        #     elif self.causal_variants[p].site.position < tree.interval.left:
        #         p += 1        
            
        #     elif self.causal_variants[p].site.position > tree.interval.right:
        #         tree.next()
        #         t += 1

        
    def diffs(self):
        cols = np.tile(self.y, (self.N, 1))
        rows = cols.T
        buffer = cols - rows
        return np.abs(buffer)
    
    def write_to_file_gcta(self, out):
        """
        write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value)

        Returns
        -------
        None.

        """
    
        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.arange(1,self.N+1)
        tmp_pheno['2'] = tmp_pheno['1']
        tmp_pheno['3'] = self.y        
        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)
