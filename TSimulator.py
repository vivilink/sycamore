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
        
    def run_simulation(self, out, randomGenerator, logfile):
        
        # simulate 500 haplotypes of chr1 of individuals from Europe, keep only 5 Mb. The species, contig, model, samples, engine and trees are all objects of stdpopsim
        # where do we specify ploidy?
        logfile.info("Simulating trees using stdPopsim (msprime). Simulating chr1 of HomSap with model OutOfAfrica_3G09.")
        
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(500, 0, 0) # Returns a list of msprime.Sample objects, with the number of samples from each population determined by the positional arguments.
        engine = stdpopsim.get_engine("msprime") #returns an engine with a "simulate" method
        trees_full = engine.simulate(model, contig, samples, seed = randomGenerator.seed) #this runs "msprime.sim_ancestry", default ploidy = 2. Extra arguments passed to simulate are passed to msprime.sim_ancestry
        trees_full.dump(out + "_full.trees")
        logfile.info("Wrote trees with interval " + str(interval) + " to " + out + "_full.trees")

        
        interval = [49e6,50e6]
        self.trees = trees_full.keep_intervals([interval], simplify=True) 
        self.trees.dump(out + ".trees")
        logfile.info("Wrote trees with interval " + str(interval) + " to " + out + ".trees")

        
        return(self.trees)
        
        
class TSimulatorMSPrime(TSimulator):
    def __init__(self):
        super().__init__()  
        
    def run_simulation(self, out, randomGenerator, logfile):
        trees_msprime = msprime.sim_ancestry(samples=250, ploidy=2, sequence_length=100000, recombination_rate=10e-8)
        self.trees = msprime.sim_mutations(trees_msprime, rate=0.01)
        self.trees.dump(out + ".trees")
        
        logfile.info("Writing trees to " + out + ".trees")

        return(self.trees)