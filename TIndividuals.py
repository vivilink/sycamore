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
        self.num_haplotypes = N
        self.num_inds = int(N / ploidy)
        self.ind_assignment = pd.DataFrame()
        self.ind_assignment['haplotypes'] = range(0,self.num_haplotypes)
        self.ind_assignment['individual'] = np.repeat(-1, self.num_inds)
        assignment =-1
        for i in range(N):    
            if i % 2 == 0:
                assignment += 1
            self.ind_assignment['individual'][i] = assignment
            
        print(self.ind_assignment)
            
    def get_individual(self, haplotype):
        if haplotype > max(self.ind_assignment['haplotypes']) or haplotype < min(self.ind_assignment['haplotypes']):
            raise ValueError("Haplotype out of bounds")
        return(self.ind_assignment['individual'][haplotype])
    
    def print_haplotype_assignment(self):
        print(self.ind_assignment)
    
    def writeShapeit2(self, out, logfile):
        
        logfile.info("- Writing individuals in Shapeit2 format to file '" + out + "_inds.sample'")
        
        haps = pd.DataFrame()  
        haps['ID_1'] = range(self.num_inds)
        if self.ploidy == 1:
            haps['ID_2'] = "NA"
        else:
            haps['ID_2'] = haps['ID_1']
        haps['missing'] = np.repeat(0, self.num_inds)
        
        #add top row of zeros, otherwise there will be one ind missing (https://myersgroup.github.io/relate/input_data.html)
        top_row = pd.DataFrame({'ID_1':[0],'ID_2':[0],'missing':[0]})
        haps = pd.concat([top_row, haps]).reset_index(drop = True)
                
        #write to file
        haps.to_csv(out + "_inds.sample", sep=' ', header=True, index=False)
