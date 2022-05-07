#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 16:50:32 2021

@author: linkv
"""

import stdpopsim
import tskit
import msprime
import TTree as tt
import math
import numpy as np


class TSimulator:
    def __init__(self):
        self.trees = None

    def read_simulation(self, file):
        self.trees = tskit.load(file)

    @staticmethod
    def simulate_more_mutations(arguments, logfile):
        """
        Add more mutations to existing tskit object

        Parameters
        ----------
        arguments : TParams
        logfile : IndentedLoggerAdapter
        
        Returns
        -------
        tskit.TreeSequence

        """
        logfile.add()
        trees = tskit.load(arguments.tree_file)

        rate = arguments.mu
        if arguments.AH_tree_pos is not None:
            focal_tree = trees.at(arguments.AH_tree_pos)
            rate = msprime.RateMap(
                position=[0, focal_tree.interval.left, focal_tree.interval.right, trees.sequence_length],
                rate=[0, arguments.mu, 0]
            )
            logfile.info("- Adding more mutations to tree covering " + str(
                arguments.AH_tree_pos) + " in tree sequence read from " + arguments.tree_file + " with rate " + str(
                arguments.mu))
        else:
            logfile.info("- Adding more mutations to trees read from " + arguments.tree_file + " with rate " + str(
                arguments.mu) + " across the whole tree sequence")

        # set discrete genome to false to assume infinite sites model to avoid triallelic sites
        trees = msprime.sim_mutations(trees, rate=rate, random_seed=arguments.seed, discrete_genome=False)
        trees.dump(arguments.out + "_moreMutations.trees")
        logfile.info("- Wrote new tree file to " + arguments.out + "_moreMutations.trees")
        logfile.sub()

        tt.TTrees.writeStats(ts_object=trees, out=arguments.out + "_moreMutations", logfile=logfile)


class TSimulatorStdPopsim(TSimulator):

    def __init__(self):
        super().__init__()

    def run_simulation(self, arguments, randomGenerator, logfile):
        """
        simulate N haplotypes of chr1 of individuals from Africa, keep only 5 Mb. The species, contig, model,
        samples, engine and trees are all objects of stdpopsim

        Parameters
        ----------
        arguments : TParams
        randomGenerator : randomGenerator
        logfile : logging

        Returns
        -------
        None.

        """
        # TODO: where do we specify ploidy?
        logfile.info(
            "- Simulating trees using stdPopsim (msprime). Simulating chr1 of HomSap with model OutOfAfrica_3G09.")

        logfile.add()

        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr1")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(arguments.N, 0, 0)
        engine = stdpopsim.get_engine("msprime")
        trees_full = engine.simulate(model, contig, samples, seed=randomGenerator.seed, discrete_genome=True)
        trees_full.dump(arguments.out + "_full.trees")
        logfile.info("- Wrote trees of full chromosome to " + arguments.out + "_full.trees")

        # extract central interval of chr1
        interval = arguments.trees_interval
        self.trees = trees_full.keep_intervals([interval], simplify=True)
        self.trees.dump(arguments.out + ".trees")
        logfile.info("- Wrote trees with interval " + str(interval) + " to " + arguments.out + ".trees")

        logfile.sub()
        return self.trees


def step_mig_mat(m, nrow, ncol):
    pmat = np.arange(0, nrow * ncol).reshape(nrow, ncol)
    mmat = np.zeros(shape=[nrow * ncol, nrow * ncol])

    def contain(ix, max_ix):
        if ix < 0:
            return 0
        if ix > (max_ix - 1):
            return max_ix - 1
        else:
            return ix

    for ii in range(nrow * ncol):
        center_ix = np.where(pmat == ii)
        top_ix = pmat[contain(center_ix[0] - 1, nrow), contain(center_ix[1], ncol)]
        bottom_ix = pmat[contain(center_ix[0] + 1, nrow), contain(center_ix[1], ncol)]
        left_ix = pmat[contain(center_ix[0], nrow), contain(center_ix[1] - 1, ncol)]
        right_ix = pmat[contain(center_ix[0], nrow), contain(center_ix[1] + 1, ncol)]

        mmat[ii, top_ix] = mmat[ii, bottom_ix] = mmat[ii, left_ix] = mmat[
            ii, right_ix
        ] = m
        mmat[top_ix, ii] = mmat[bottom_ix, ii] = mmat[left_ix, ii] = mmat[
            right_ix, ii
        ] = m
        mmat[ii, ii] = 0

    return mmat


class TSimulatorMSPrime(TSimulator):
    def __init__(self):
        super().__init__()

    def run_simulation(self, arguments, randomGenerator, logfile):
        logfile.info("- Simulating trees using msprime.")
        logfile.add()

        trees_msprime = msprime.sim_ancestry(samples=arguments.N, ploidy=arguments.ploidy,
                                             sequence_length=arguments.sequence_length,
                                             recombination_rate=arguments.recomb_rate)
        self.trees = msprime.sim_mutations(trees_msprime, rate=arguments.mu, random_seed=arguments.seed)
        self.trees.dump(arguments.out + ".trees")

        logfile.info("Writing trees to " + arguments.out + ".trees")
        logfile.sub()

        return self.trees

