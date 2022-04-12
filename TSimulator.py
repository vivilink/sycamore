#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:50:32 2021

@author: linkv
"""

import stdpopsim
import tskit
import msprime


class TSimulator:
    def __init__(self):
        self.trees = None
        
    
    def read_simulation(self, file):
        self.trees = tskit.load(file)

        
class TSimulatorStdPopsim(TSimulator):
    
    def __init__(self):
        super().__init__()        
        
    def run_simulation(self, arguments, randomGenerator, logfile):
        """
        simulate N haplotypes of chr1 of individuals from Africa, keep only 5 Mb. The species, contig, model, samples, engine and trees are all objects of stdpopsim

        Parameters
        ----------
        arguments : TParams
        randomGenerator : randomGenerator
        logfile : logging

        Returns
        -------
        None.

        """       
        #TODO: where do we specify ploidy?
        logfile.info("- Simulating trees using stdPopsim (msprime). Simulating chr1 of HomSap with model OutOfAfrica_3G09.")
        
        logfile.add()
        
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(arguments.N, 0, 0) # Returns a list of msprime.Sample objects, with the number of samples from each population determined by the positional arguments.
        engine = stdpopsim.get_engine("msprime") #returns an engine with a "simulate" method
        trees_full = engine.simulate(model, contig, samples, seed = randomGenerator.seed, discrete_genome=True) #this runs "msprime.sim_ancestry", default ploidy = 2. Extra arguments passed to simulate are passed to msprime.sim_ancestry
        trees_full.dump(arguments.out + "_full.trees")
        logfile.info("- Wrote trees of full chromosome to " + arguments.out + "_full.trees")

        #TODO: do not hardcode interval!
        interval = [49e6,50e6]
        self.trees = trees_full.keep_intervals([interval], simplify=True) 
        self.trees.dump(arguments.out + ".trees")
        logfile.info("- Wrote trees with interval " + str(interval) + " to " + arguments.out + ".trees")        
        logfile.sub()
        
        return(self.trees)
    
        
        
class TSimulatorMSPrime(TSimulator):
    def __init__(self):
        super().__init__()  
        
    def run_simulation(self, arguments, randomGenerator, logfile):
        logfile.info("- Simulating trees using msprime.")
        logfile.add()

        trees_msprime = msprime.sim_ancestry(samples=arguments.N, ploidy=arguments.ploidy, sequence_length=arguments.sequence_length, recombination_rate=10e-8)
        self.trees = msprime.sim_mutations(trees_msprime, rate=arguments.mu)
        self.trees.dump(arguments.out + ".trees")
        
        logfile.info("Writing trees to " + arguments.out + ".trees")
        logfile.sub()

        return(self.trees)