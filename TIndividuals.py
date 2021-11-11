#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:43:16 2021

@author: linkv
"""

import pandas as pd
import numpy as np


class Individuals:
    def __init__(self, ploidy, N):
        self.ploidy = ploidy
        self.num_inds = int(N / ploidy)
        self.ind_assignment = pd.DataFrame()
        self.ind_assignment['haplotypes'] = range(0,N)
        self.ind_assignment['individual'] = np.repeat(-1, N)
        assignment =-1
        for i in range(N):    
            if i % 2 == 0:
                assignment += 1
            self.ind_assignment['individual'][i] = assignment
            
    def get_individual(self, haplotype):
        if haplotype > max(self.ind_assignment['haplotypes']) or haplotype < min(self.ind_assignment['haplotypes']):
            raise ValueError("Haplotype out of bounds")
        return(self.ind_assignment['individual'][haplotype])
    
    def writeShapeit2(self, out, logfile):
        
        logfile.info("- Writing individuals in Shapeit2 format to file '" + out + "_inds.sample'")
        
        haps = pd.DataFrame()  
        haps['ID_1'] = range(self.num_inds)
        haps['ID_2'] = "NA"
        haps['missing'] = np.repeat(0, self.num_inds)
        # haps.iloc[0] = [0, 'F','Piyu']
                
        haps.to_csv(out + "_inds.sample", sep=' ', header=True, index=False)
