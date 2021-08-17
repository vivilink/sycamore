#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np


class Phenotypes:
    def __init__(self, N, num_variants):
        self.N = N
        self.num_variants = num_variants
        self.y = [0] * N
        self.betas = [0] * num_variants
        self.causal_positions = []
        
        self.filled = False
        
    def simulateUniform(self, ts_object_variants, prop_causal_mutations, sd_environmental_noise, sd_beta_causal_mutations):
        y = np.random.normal(loc=0, scale=sd_environmental_noise, size=self.N)

        #add phenotypic effect to mutations that are uniformly distributed
        causal_positions = []
        variant_positions = []
        for v, var in enumerate(ts_object_variants):  

            variant_positions.append(var.site.position)
        
            #causal mutation
            r = np.random.uniform(0,1)
            if(r < prop_causal_mutations):
                
                #define beta
                beta = np.random.normal(loc=0, scale=sd_beta_causal_mutations, size=1)
                self.betas[v] = beta
                
                #simulate phenotype
                y[var.genotypes == 1] += beta
                
                #save causal position
                causal_positions.append(var.site.position)
            
        self.filled = True