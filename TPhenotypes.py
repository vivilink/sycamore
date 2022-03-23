#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np
import pandas as pd

#TODO: being typed or not should be an option for all causal variants

class Phenotypes:
    def __init__(self, name, variants, inds, logfile):
        self.name = name
        self.N = inds.num_inds
        self.y = np.zeros(self.N)
        self.betas = [0] * variants.number
        self.causal_variants = []
        self.causal_betas = []
        self.causal_power = []
        self.causal_trees = []
        self.causal_variant_indeces = []
        self.causal_tree_indeces = []
        self.filled = False
        
        logfile.info("- Created '" + self.name + "' phenotype object")
        
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
        
        causal_variants = [variants.variants[i] for i in causal_variant_indeces]
        causal_pos = [variants.variants[i].site.position for i in causal_variant_indeces]
                
        if(len(causal_variants) != len(betas)):
            raise ValueError("must provide equal number of causal variants and betas to simulate fixed phenotype")
                    
        for v, var in enumerate(causal_variants):
            self.betas[causal_variant_indeces[v]] = betas[v]            
            self.y[var.genotypes == 1] += betas[v]
            self.causal_variants.append(var)
            self.causal_betas.append(betas[v])
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            self.causal_power.append(betas[v]**2 * allele_freq * (1-allele_freq))            
            self.causal_variant_indeces.append(causal_variant_indeces[v])
            
            logfile.info("- Simulated causal variant at position " + str(causal_pos[v]) + " at index " + str(causal_variant_indeces[v]) + " with beta " + str(round(betas[v], 3)) + " and allele freq " + str(allele_freq) + " resulting in a power of " + str(round(betas[v]**2 * allele_freq * (1-allele_freq), 3)))
        
    def simulateUniform(self, variants, prop_causal_mutations, sd_beta_causal_mutations, random, logfile, mean_beta_causal_mutation = 0):
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
        for v, var in enumerate(variants.variants): 
            
            # if variants.info["typed"] == True:
            
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
        
        logfile.info("- Simulated phenotypes based on " + str(len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True
        
    def simulateCausalRegion(self, variants, left_bound, right_bound, sd_beta_causal_mutations, random, logfile):
        #add phenotypic effect to mutations that are uniformly distributed
        for v, var in enumerate(variants.variants): 
            # is the variant in the tree
            if left_bound <= variants.info.loc[v]['position'] <= right_bound:
                
                #define beta
                beta = random.random.normal(loc=0, scale = sd_beta_causal_mutations, size=1)[0]
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

        logfile.info("- Simulated phenotypes based on " + str(len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True

    
    def findCausalTrees(self, ts_object):
        for v in self.causal_variants:
            causal_tree = ts_object.at(v.site.position)
            self.causal_tree_indeces.append(causal_tree.get_index())

        
    def diffs(self):
        cols = np.tile(self.y, (self.N, 1))
        rows = cols.T
        buffer = cols - rows
        return np.abs(buffer)
    
    def write_to_file_gcta(self, out, logfile):
        """
        write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value)

        Returns
        -------
        None.

        """
        logfile.info("- Writing phenotype data in gcta format to '" + out + "_phenotypes.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.arange(1,self.N+1)
        tmp_pheno['2'] = tmp_pheno['1']
        tmp_pheno['3'] = self.y        
        
        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)
        

    def write_to_file(self, variants, out, logfile):
        """
        write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value)

        Returns
        -------
        None.

        """
    
        #results for each variant
        table = pd.DataFrame()
        table['start'] = variants.info['position']
        table['end'] = variants.info['position']
        table['allele_freq'] = variants.info['allele_freq']
        table['typed'] = variants.info['typed']
        table['causal'] = np.repeat("FALSE", variants.number)
        table.loc[self.causal_variant_indeces, 'causal'] = "TRUE"
        table['betas'] = self.betas 
        table['power'] = 0
        table.loc[self.causal_variant_indeces, 'power'] = self.causal_power
        table['genotypic_var'] = np.array(self.betas) * np.array(self.betas) * np.array(table['allele_freq']) * (np.repeat(1, variants.number) - np.array(table['allele_freq']))
        table['phenotypic_var'] = np.var(self.y)
      
        logfile.info("- Writing phenotype data '" + out + "_pheno_causal_vars.csv'")
        table.to_csv(out + "_pheno_causal_vars.csv", index = False, header = True)       
