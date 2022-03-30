#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

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
        
        
    def simulate(self, args, r, logger, variants_orig, inds, trees, plots_dir):
        
        if args.pty_sim_method is None:
            raise ValueError("Must provide a phenotype simulation method with --pty_sim_method")

        if args.pty_sim_method == 'uniform':
            logger.info("- Simulating phenotypes based on uniformly chosen variants with prop_causal_mutations: " + str(args.pty_prop_causal_mutations) + " and sd_beta_causal_mutations: " + str(args.pty_sd_beta_causal_mutations)) 
            self.simulateUniform(variants_orig, inds, prop_causal_mutations=args.pty_prop_causal_mutations, sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r, logfile=logger)
     
        elif args.pty_sim_method == 'fixed':
            if args.pty_fixed_betas == None:
                raise ValueError("Must provide beta values provided for 'fixed' phenotype using '--pty_fixed_betas'")

            logger.info("- Simulating phenotypes based on the following indeces: " + str(args.pty_fixed_variant_indeces) + " and the following betas: " + str(args.pty_fixed_betas)) 
            self.simulateFixed(variants_orig, inds, args.pty_fixed_variant_indeces, args.pty_fixed_betas, logger)

        elif args.pty_sim_method == 'singleTyped':
            if args.pty_fixed_betas == None:
                raise ValueError("Must provide beta values for 'singleTyped' phenotype using '--pty_fixed_betas'")
            if args.single_variant_af == None:
                raise ValueError("Must provide allele freq values for 'singleTyped' phenotype using '--single_variant_af'")
                
            fig, ax = plt.subplots(1,figsize=(30,30))            
            var_index, pos = variants_orig.findVariant(typed=True, freq = args.single_variant_af, interval = args.single_variant_interval, out = args.out, subplot = ax, random = r, logfile = logger)
            fig.tight_layout()
            fig.set_size_inches(30, 30)
            fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logger.info("- Simulating a phenotypes based on the following typed variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
            self.simulateFixed(variants_orig, inds, [var_index], args.pty_fixed_betas, logger)
            
        elif args.pty_sim_method == 'singleUntyped':
            if args.pty_fixed_betas == None:
                raise ValueError("Must provide beta values for phenotype 'singleUntyped' using '--pty_fixed_betas'")
            if args.single_variant_af == None:
                raise ValueError("Must provide allele freq values for 'singleTyped' phenotype using '--single_variant_af'")
                
            fig, ax = plt.subplots(1,figsize=(30,30))            
            var_index, pos = variants_orig.findVariant(typed=False, freq = args.single_variant_af, interval = args.single_variant_interval, out = args.out, subplot = ax, random = r, logfile = logger)   
            fig.tight_layout()
            fig.set_size_inches(30, 30)
            fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')
     
            logger.info("- Simulating a phenotypes based on the following untyped variant index: " + str(var_index) + " at position " +  str(variants_orig.info['position'][var_index]) + " with allele freq " + str(variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(args.pty_fixed_betas)) 
            #to know which variants are untyped you need variants from simulated tree, not estimated tree
            if args.variants_file is None:
                raise ValueError("Must provide file with untyped variants to simulate phenotype with 'singleUntyped' model")
            self.simulateFixed(variants_orig, inds, [var_index], args.pty_fixed_betas, logger)
            
        elif args.pty_sim_method == 'oneTree':
            causal_tree = trees.at(args.causal_tree_pos)
            logger.info("- Simulating phenotypes based on all variants of the tree covering postion " + str(args.causal_tree_pos))
            if args.pty_sd_beta_causal_mutations is None:
                raise ValueError("pty_sd_beta_causal_mutations must be set to simulate phenotype with method 'oneTree'")
            self.simulateCausalRegion(variants_orig, inds, left_bound = causal_tree.interval.left, right_bound = causal_tree.interval.right, sd_beta_causal_mutations = args.pty_sd_beta_causal_mutations, random = r, logfile = logger)

        elif args.pty_sim_method == 'allelicHetero':
            if args.allelic_hetero_file == None:
                raise ValueError("No instruction file provided for allelic heterogeneity simulation")
            if args.single_variant_af != None or args.pty_fixed_betas != None:
                raise ValueError("Provided allele frequency or beta value as well as instruction file for allelic heterogeneity. Can accept only one type of instructions.")
            
            ah_info = pd.read_csv(args.allelic_hetero_file, delimiter = "\t")
            variant_indeces = []
            fixed_betas = []
            sum_betas = 0
            
            logger.info("- Searching for loci with requested allele frequencies.")
            logger.add()
            
            #start plot 
            # TODO: this plotting should not be done here but in a function instead
            fig, ax = plt.subplots(ah_info.shape[0],figsize=(30,30))        
            for index, row in ah_info.iterrows():
                #get allele freq
                f = -1
                if row['freq'] > 0.5:
                    f = 1 - row['freq']
                    logger.warning("- Allele frequencies above 0.5 are not allowed. Transformed " + str(row['freq'])  + " to " + str(f) + ".")
                else:
                    f = row['freq']
                
                #get beta
                if not np.isnan(row['beta']) and not np.isnan(row['power']):
                    raise ValueError("Cannot fix power and beta value. One value must be set to 'NA'. Beta is " + str(row['beta']) + " and power is " + str(row['power']))
                if np.isnan(row['beta']):
                    beta = np.sqrt(row['power'] / (f * (1 - f)))
                    #some betas should be negative
                    r_num = r.random.uniform(0,1,1)
                    if r_num < 0.5:
                        beta = -beta
                else:
                    beta = row['beta']


                var_index, pos = variants_orig.findVariant(typed=False, freq = f, interval = [row["interval_start"], row["interval_end"]], out = args.out, subplot = ax[index], random = r, logfile = logger)   
                variant_indeces.append(var_index)
                fixed_betas.append(beta)
                sum_betas += beta

                fig.tight_layout()
                fig.set_size_inches(30, 30)
                fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logger.sub()

            logger.info("- Simulating allelic heterogeneous phenotype with total beta " + str(sum_betas))
            
            logger.info("- Simulating phenotypes:")
            logger.add()
            self.simulateFixed(variants_orig, inds, variant_indeces, fixed_betas, logger)
            logger.sub()
        
        
        #write phenotypes to file
        logger.info("- Simulating random noise with sd " + str(args.pty_sd_envNoise))
        
        # count = (self.y == 0).sum()
        # print("number of zeros in y", count)
        
        self.simulateEnvNoise(args.pty_sd_envNoise, r)
        # self.standardize(logger)
        self.write_to_file(variants_orig, inds, args.out, logger)
        
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
        self.y += random.random.normal(loc=0, scale=sd_environmental_noise, size=self.N)
        self.filled = True

    def simulateFixed(self, variants, inds, causal_variant_indeces, betas, logfile):
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
            #define beta
            self.betas[causal_variant_indeces[v]] = betas[v]      
            
            #simulate phenotype
            if inds.ploidy == 1:
                self.y[var.genotypes == 1] += betas[v]
            else:
                genotypes = inds.get_diploid_genotypes(var.genotypes)
                self.y[genotypes == 1] += betas[v]
                self.y[genotypes == 2] += 2 * betas[v]
            
            #save causal position
            self.causal_variants.append(var)
            self.causal_betas.append(betas[v])
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            self.causal_power.append(betas[v]**2 * allele_freq * (1-allele_freq))            
            self.causal_variant_indeces.append(causal_variant_indeces[v])
            
            logfile.info("- Simulated causal variant at position " + str(causal_pos[v]) + " at index " + str(causal_variant_indeces[v]) + " with beta " + str(round(betas[v], 3)) + " and allele freq " + str(allele_freq) + " resulting in a power of " + str(round(betas[v]**2 * allele_freq * (1-allele_freq), 3)))
        
    def simulateUniform(self, variants, inds, prop_causal_mutations, sd_beta_causal_mutations, random, logfile, mean_beta_causal_mutation = 0):
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
                    if inds.ploidy == 1:
                        self.y[var.genotypes == 1] += beta
                        self.y[var.genotypes == 2] += 2 * beta
                    else:
                        genotypes = inds.get_diploid_genotypes(var.genotypes)
                        self.y[genotypes == 1] += beta
                        self.y[genotypes == 2] += 2 * beta
    
                    #save causal position
                    self.causal_variants.append(var)
                    self.causal_betas.append(beta)
                    allele_freq = variants.info['allele_freq'][v]
                    self.causal_power.append(beta**2 * allele_freq * (1-allele_freq))
                    self.causal_variant_indeces.append(v)
        
        logfile.info("- Simulated phenotypes based on " + str(len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True
        
    def simulateCausalRegion(self, variants, inds, left_bound, right_bound, sd_beta_causal_mutations, random, logfile):
        #add phenotypic effect to mutations that are uniformly distributed
        for v, var in enumerate(variants.variants): 
            # is the variant in the tree
            if left_bound <= variants.info.loc[v]['position'] <= right_bound:
                
                #define beta
                beta = random.random.normal(loc=0, scale = sd_beta_causal_mutations, size=1)[0]
                self.betas[v] = beta
                
                #simulate phenotype
                if inds.ploidy == 1:
                    self.y[var.genotypes == 1] += beta
                    self.y[var.genotypes == 2] += 2 * beta
                else:
                    genotypes = inds.get_diploid_genotypes(var.genotypes)
                    self.y[genotypes == 1] += beta
                    self.y[genotypes == 2] += 2 * beta

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
    
    def standardize(self, logger):
        logger.info("- Standardizing phenotypes")
        self.y = (self.y - np.mean(self.y)) / np.std(self.y)
    
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
        
    def write_to_file(self, variants, inds, out, logfile):
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
        table['genotypic_var'] = np.repeat(float(inds.ploidy), variants.number) * np.array(self.betas) * np.array(self.betas) * np.array(table['allele_freq']) * (np.repeat(1.0, variants.number) - np.array(table['allele_freq']))
        table['phenotypic_var'] = np.var(self.y)
      
        logfile.info("- Writing phenotype data '" + out + "_pheno_causal_vars.csv'")
        table.to_csv(out + "_pheno_causal_vars.csv", index = False, header = True)       
