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
        logfile.info("- Simulating trees using msprime:")
        logfile.add()

        # get recombination rates
        try:
            recomb_obj = float(arguments.recomb_rate)
        except ValueError:
            recomb_obj = msprime.RateMap.read_hapmap(arguments.recomb_rate)
            # recomb_obj = recomb_obj.slice(left=190000000, right=200000000, trim=True)
            if arguments.recomb_map_start_random is False:
                recomb_map_start = recomb_obj.left[1]  # first position in map file
                recomb_map_end = recomb_obj.sequence_length  # last position in map file
            else:
                if recomb_obj.sequence_length - arguments.sequence_length < 0:
                    raise ValueError("Cannot extract slice with length " + str(arguments.sequence_length)
                                     + " from recombination map. Map too short")
                recomb_map_start = randomGenerator.random.randint(recomb_obj.left[1],
                                                                  recomb_obj.sequence_length - arguments.sequence_length)
                recomb_map_end = recomb_map_start + arguments.sequence_length

            recomb_obj = recomb_obj.slice(left=recomb_map_start, right=recomb_map_end, trim=True)
            logfile.info("- Simulating with recombination map read from " + arguments.recomb_rate + " with length "
                         + str(recomb_obj.sequence_length) + ", starting pos " + str(recomb_map_start)
                         + " and ending pos " + str(recomb_map_end))

        # get mutation rates
        if arguments.mu is not None:
            mut_obj = float(arguments.mu)
            logfile.info("- Simulating fixed mutation rate of " + str(arguments.mu))
        else:
            num_windows = int(np.floor(arguments.sequence_length / arguments.mut_window_size))
            positions = [0 + w * arguments.mut_window_size for w in range(num_windows + 1)]
            rates = [randomGenerator.random.beta(arguments.mut_beta_shape1, arguments.mut_beta_shape2) for w in range(num_windows)]
            # rates = arguments.mu
            # print(rates)
            logfile.info("- Simulating mutation rates in windows with length " + str(arguments.mut_window_size)
                         + ", mutation rate mean " + str(np.mean(rates)) + " and variance " + str(np.var(rates)))
            mut_obj = msprime.RateMap(
                position=positions,
                rate=rates
            )

        # simulate
        trees_msprime = msprime.sim_ancestry(samples=arguments.N, ploidy=1, # do not simulate ploidy here!
                                             sequence_length=arguments.sequence_length,
                                             recombination_rate=recomb_obj,
                                             population_size=arguments.population_size)
        self.trees = msprime.sim_mutations(trees_msprime, rate=mut_obj, random_seed=arguments.seed)
        self.trees.dump(arguments.out + ".trees")

        logfile.info("- Writing trees to " + arguments.out + ".trees")
        logfile.sub()
        return self.trees
