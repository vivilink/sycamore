#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:44:45 2021

@author: linkv
"""
import numpy as np
import pandas as pd


class TVariants:
    def __init__(self, ts_object, samp_ids):
        self.number = len(list(ts_object.variants(samples=samp_ids)))
        self.positions = np.empty(self.number)
        self.alleleFreq = np.empty(self.number)
        
        for v, var in enumerate(list(ts_object.variants(samples=samp_ids))):  
            tmp = sum(var.genotypes) / len(var.genotypes)
            af = min(tmp, 1-tmp)
            
            self.positions[v] = var.site.position
            self.alleleFreq[v] = af

        
    def writeVariantInfo(self, ts_object, samp_ids, name):
        table = pd.DataFrame()
        table['variant_index'] = np.arange(self.number)
        table['position'] = self.positions
        table['allele_freq'] = self.alleleFreq
        
        table.to_csv(name + "_simulated_sample_variants.csv", header=True, index=False)
        
class TVariantsFiltered(TVariants):
    
    def __init__(self, ts_object, samp_ids, min_allele_freq, max_allele_freq, prop_typed_variants, random):
        if min_allele_freq < 0 or max_allele_freq < 0 or min_allele_freq > 1 or max_allele_freq > 1 or min_allele_freq > max_allele_freq:
            raise ValueError("allele frequency filters are nonsensical")
        
        #initialize
        # self.number_tot = len(list(ts_object.variants(samples=samp_ids)))
        self.number = -1
        positions = np.empty(0)
        variants = np.empty(0)
        allele_frequencies = np.empty(0)
        typed = np.empty(0)
        
        
        #fill by filtering, can't directly fill into info table because number of variants is unknown
        for v, var in enumerate(list(ts_object.variants(samples=samp_ids))):  
            tmp = sum(var.genotypes) / len(var.genotypes)
            af = min(tmp, 1-tmp)
            
            if af >= min_allele_freq and af <= max_allele_freq:
                variants = np.append(variants, var)
                allele_frequencies = np.append(allele_frequencies, af)
                positions = np.append(positions, var.site.position)
                
                if prop_typed_variants == 1:
                    typed = np.append(typed, True)
                else:
                    r = random.random.uniform(0,1,1)
                    if r < prop_typed_variants:
                        typed = np.append(typed, True)
                    else: 
                        typed = np.append(typed, False)
                    
        
        self.number = len(variants)    
        
        # fill info into info table
        self.info = pd.DataFrame({'index': np.arange(0, self.number), 'position': positions, 'variant': variants, 'allele_freq': allele_frequencies, 'typed': typed})
                
    
    def print_genotypes(self, index):        
        file = "genotypes_variant" + str(index) + ".txt"
        self.info['variant'][index].genotypes.tofile(file=file)
        
    def fill_diploidGenotypes(self, individuals):
        for v in self.info['variant']:
            v.site.metadata = []
            
            # for h in individuals.ind_assignment['haplotypes']:
            
    def writeVariantInfo(self, name):    
        self.info.drop(columns='variant').to_csv(name + "_filtered_sample_variants.csv", header=True, index=False)
        
    def findVariant(self, typed, freq, interval, logfile):
        """
        Find a variant that fits criteria to simulate fixed genotype depending on only one variant

        Parameters
        ----------
        typed : bool
            Should the variant returned be typed or not.
        freq : float
            Requested allele frequency of the variant. If there is no variant with this exact frequency within the requested interval, 0.001 will be deducted from the freq repeatedly until a variant is found.
        interval : list of floats
            Requested genomic interval within which the variant should be.

        Raises
        ------
        ValueError
            If the search for a variant within the requested genomic interval ends up with the allele freq being negative, the search is stopped an an error is thrown.

        Returns
        -------
        Index of the variant found, can be used to simulate fixed phenotype.

        """
        # check if interval is valid
        if self.info[(self.info['position'] >= interval[0]) & (self.info['position'] <= interval[1])].shape[0] == 0:
            raise ValueError("The interval " + str(interval) + " contains no variants")
            
        # check if there are variants with requested typed status
        if self.info[self.info['typed'] == typed].shape[0] == 0:
            raise ValueError("There are no variants of typed status '" + str(typed) +"'")
 
        info = self.info[(self.info['typed'] == typed) & (self.info['allele_freq'] == freq) & (self.info['position'] >= interval[0]) & (self.info['position'] <= interval[1])]
        while info.shape[0] < 1:
            freq = round(freq - 0.001,3)
            info = self.info[(self.info['typed'] == typed) & (self.info['allele_freq'] == freq) & (self.info['position'] >= interval[0]) & (self.info['position'] <= interval[1])]
            if freq < 0:
                raise ValueError("allele frequency became negative while searching for one in interval " + str(interval))
                
        logfile.info("- Found variant with freq " + str(freq) + " within the following interval: " + str(interval))
        return info.iloc[0]['index']
        
