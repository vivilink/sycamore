#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:44:45 2021

@author: linkv
"""
import numpy as np


class TVariants:
    def __init__(self):
        self.number = -1
        self.positions = np.empty(0)
        
class TVariantsSamples(TVariants):
    
    def __init__(self, ts_object, samp_ids, min_allele_freq, max_allele_freq):
        if min_allele_freq < 0 or max_allele_freq < 0 or min_allele_freq > 1 or max_allele_freq > 1 or min_allele_freq > max_allele_freq:
            raise ValueError("allele frequency filters are nonsensical")
        
        #initialize
        # self.number_tot = len(list(ts_object.variants(samples=samp_ids)))
        self.number = -1
        self.positions = []
        self.variants = []
        self.allele_frequencies = []
        
        #fill by filtering
        for v, var in enumerate(list(ts_object.variants(samples=samp_ids))):  
            tmp = sum(var.genotypes) / len(var.genotypes)
            af = min(tmp, 1-tmp)
            
            if af >= min_allele_freq and af <= max_allele_freq:
                self.variants.append(var)
                self.allele_frequencies.append(af)
                self.positions.append(var.site.position)

        self.number = len(self.variants)    
        self.allele_frequencies = np.array(self.allele_frequencies)
        
    
    def printGenotypes(self, index):        
        file = "genotypes_variant" + str(index) + ".txt"
        self.variants[index].genotypes.tofile(file=file)
        
        
