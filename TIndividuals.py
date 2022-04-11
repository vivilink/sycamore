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
        
        if ploidy is None:
            raise ValueError("Must provide ploidy with --ploidy")
        
        self.ploidy = ploidy
        self.num_haplotypes = N
        self.num_inds = int(self.num_haplotypes / ploidy)
        self.ind_assignment = pd.DataFrame()
        self.ind_assignment['haplotype'] = range(0,self.num_haplotypes)
        self.ind_assignment['individual'] = np.repeat(-1, self.num_haplotypes)
        assignment =-1
        for i in range(self.num_haplotypes):    
            if i % 2 == 0:
                assignment += 1
            self.ind_assignment['individual'][i] = assignment
                        
    def get_individual(self, haplotype):
        if haplotype > max(self.ind_assignment['haplotype']) or haplotype < min(self.ind_assignment['haplotype']):
            raise ValueError("Haplotype out of bounds")
        return(self.ind_assignment['individual'][haplotype])['individual']
    
    def get_haplotypes(self, individual):
        if individual > self.num_inds or individual < 0:
            raise ValueError("Individual out of bounds")        
    
        tmp = self.ind_assignment['haplotype'].values[self.ind_assignment['individual'] == individual]
        return(tmp)
                
    def get_diploid_genotypes(self, haploid_genotypes):
        table = pd.DataFrame()
        table['individual'] = self.ind_assignment['individual']
        table['haploid_genotypes'] = haploid_genotypes
        
        table = table.groupby('individual').agg(
            diploid_genotypes = pd.NamedAgg(column='haploid_genotypes', aggfunc=sum)
        )
        return(table['diploid_genotypes'])
        
    
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
