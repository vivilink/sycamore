#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:37:12 2021

@author: linkv
"""
import numpy as np

class GWAS:
    
    def __init__(self, genome):
        self.genome = genome
        self.p_values = []
        self.p_values_init = False
        
    def q_values(self):
        return -np.log10(self.p_values)
        

class pGWAS(GWAS):
    
    def length(self):
        return self.genome.num_variants
    
    
    
    
class GWAS(GWAS):
    
    def length(self):
        return self.genome.num_trees