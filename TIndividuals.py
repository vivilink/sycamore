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

        self._ploidy = ploidy
        self._num_haplotypes = N
        self._num_inds = int(self._num_haplotypes / self._ploidy)
        self._ind_assignment = pd.DataFrame()
        self._ind_assignment['haplotype'] = range(0, self._num_haplotypes)
        self._ind_assignment['individual'] = np.repeat(-1, self._num_haplotypes)
        # TODO: population assignment should come from tree file!
        # self._ind_assignment['population'] = np.append(np.repeat("sample", N_sample_pop), np.repeat("ref", N_ref_pop))

        assignment = -1
        for i in range(self._num_haplotypes):
            if i % 2 == 0:
                assignment += 1
            self._ind_assignment['individual'][i] = assignment
        self._names = ["id_" + str(i) for i in np.arange(0, self._num_inds)]

    @property
    def ploidy(self):
        return self._ploidy

    @ploidy.setter
    def ploidy(self, ploidy):
        self._ploidy = ploidy

    @property
    def num_haplotypes(self):
        return self._num_haplotypes

    @num_haplotypes.setter
    def num_haplotypes(self, num_haplotypes):
        self._num_haplotypes = num_haplotypes

    @property
    def num_inds(self):
        return self._num_inds

    @num_inds.setter
    def num_inds(self, num_inds):
        self._num_inds = num_inds

    @property
    def names(self):
        return self._names

    @names.setter
    def names(self, names):
        self._names = names

    def get_individual(self, haplotype):
        if haplotype > max(self._ind_assignment['haplotype']) or haplotype < min(self._ind_assignment['haplotype']):
            raise ValueError("Haplotype out of bounds")
        return (self._ind_assignment['individual'][haplotype])['individual']

    def get_haplotypes(self, individual):
        if individual > self._num_inds or individual < 0:
            raise ValueError("Individual out of bounds")

        tmp = self._ind_assignment['haplotype'].values[self._ind_assignment['individual'] == individual]
        return (tmp)

    def get_diploid_genotypes(self, haploid_genotypes):
        table = pd.DataFrame()
        table['individual'] = self._ind_assignment['individual']
        table['haploid_genotypes'] = haploid_genotypes

        table = table.groupby('individual').agg(
            diploid_genotypes=pd.NamedAgg(column='haploid_genotypes', aggfunc=sum)
        )
        return table['diploid_genotypes']

    def write_shapeit2(self, out, logfile):
        logfile.info("- Writing individuals in Shapeit2 format to file '" + out + "_inds.sample'")

        haps = pd.DataFrame()
        haps['ID_1'] = range(self._num_inds)
        if self._ploidy == 1:
            haps['ID_2'] = "NA"
        else:
            haps['ID_2'] = haps['ID_1']
        haps['missing'] = np.repeat(0, self._num_inds)

        # add top row of zeros, otherwise there will be one ind missing (
        # https://myersgroup.github.io/relate/input_data.html)
        top_row = pd.DataFrame({'ID_1': [0], 'ID_2': [0], 'missing': [0]})
        haps = pd.concat([top_row, haps]).reset_index(drop=True)

        # write to file
        haps.to_csv(out + "_inds.sample", sep=' ', header=True, index=False)
