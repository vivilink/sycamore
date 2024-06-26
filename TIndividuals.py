#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 15:43:16 2021

@author: linkv
"""

import pandas as pd
import numpy as np
from python_log_indenter import IndentedLoggerAdapter

class Individuals:
    def __init__(self, ploidy: int, num_haplotypes: int, logfile: IndentedLoggerAdapter, relate_sample_names_file: str, phenotype_sample_file: str):
        """
        Initialize individuals. NAmes come either from relate file, from phen file or are generated automatically
        @param ploidy:
        @param num_haplotypes:
        @param logfile:
        @param relate_sample_names_file: only needed for real data, can be None if sample names should just be generated
        """

        logfile.info("- Initializing individuals with ploidy " + str(ploidy) + ". Found " + str(num_haplotypes) + " haplotypes.")

        if ploidy is None:
            raise ValueError("Must provide ploidy with --ploidy")

        self._ploidy = ploidy
        self._num_haplotypes = num_haplotypes
        self._num_inds = int(self._num_haplotypes / self._ploidy)
        self._ind_assignment = pd.DataFrame()
        self._ind_assignment['haplotype'] = range(0, self._num_haplotypes)
        self._ind_assignment['individual'] = np.repeat(-1, self._num_haplotypes)
        self._ind_has_phenotype = None
        # TODO: population assignment should come from tree file!
        # self._ind_assignment['population'] = np.append(np.repeat("sample", N_sample_pop), np.repeat("ref", N_ref_pop))

        assignment = -1
        for i in range(self._num_haplotypes):
            if i % 2 == 0:
                assignment += 1
            self._ind_assignment['individual'][i] = assignment
        if relate_sample_names_file:
            names_in_file = pd.read_csv(relate_sample_names_file, sep=' ').iloc[1:, :]
            if len(names_in_file['ID_1']) != self._num_inds:
                raise ValueError("There are " + str(self._num_inds) + " haplotypes in tree but " +
                                 str(len(names_in_file['ID_1'])) + " in " + relate_sample_names_file)
            self._names = np.array(names_in_file['ID_1'])
        elif phenotype_sample_file:
            pheno_df = pd.read_csv(phenotype_sample_file, sep=' ', header=None)
            num_phenotypes = pheno_df.shape[1] - 2
            header = ["pop", "ID"]
            for i in range(1, num_phenotypes + 1):
                header.append("phenotype" + str(i))
            pheno_df.columns = header
            self._names = np.array(pheno_df['ID'])
        else:
            self._names = np.array(["id_" + str(i) for i in np.arange(0, self._num_inds)])

        if len(self._names) != self._num_inds:
            logfile.warning("Length of ind names is " + str(len(self._names)) + " but number of individuals found in tree is " + str(self._num_inds) + ". This should only happen in task 'removeUnsampledInds'")




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

    @property
    def ind_has_phenotype(self):
        return self._ind_has_phenotype

    @ind_has_phenotype.setter
    def ind_has_phenotype(self, ind_has_phenotype):
        self._ind_has_phenotype = ind_has_phenotype

    def get_individual(self, haplotype):
        if haplotype > max(self._ind_assignment['haplotype']) or haplotype < min(self._ind_assignment['haplotype']):
            raise ValueError("Haplotype out of bounds")
        return self._ind_assignment['individual'][haplotype]

    def get_haplotypes(self, individual):
        if individual > self._num_inds or individual < 0:
            raise ValueError("Individual out of bounds")

        tmp = self._ind_assignment['haplotype'].values[self._ind_assignment['individual'] == individual]
        return tmp

    def get_diploid_genotypes(self, haploid_genotypes):
        table = pd.DataFrame()
        table['individual'] = self._ind_assignment['individual']
        table['haploid_genotypes'] = haploid_genotypes

        table = table.groupby('individual').agg(
            diploid_genotypes=pd.NamedAgg(column='haploid_genotypes', aggfunc="sum")
        )
        return table['diploid_genotypes']

    def get_indeces_inds_no_phenotype(self):
        # if self._ind_has_phenotype is None:
        #     raise ValueError("Do not know which individuals have missing phenotypes. Initiate _ind_has_phenotype!")
        return list(np.arange(0, self._num_inds, 1, dtype=int)[self._ind_has_phenotype == False])

    def write_shapeit2_relate(self, out, logfile):
        logfile.info("- Writing individuals in Shapeit2 format to file '" + out + "_relate.sample'")

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
        haps.to_csv(out + "_relate.sample", sep=' ', header=True, index=False)

    def write_sample_argNeedle(self, out, logfile):
        """
        https://www.cog-genomics.org/plink/2.0/formats#sample
        :param out:
        :param logfile:
        :return:
        """
        logfile.info("- Writing individuals in Oxford sample format to file '" + out + "_argNeedle.sample'")

        haps = pd.DataFrame()
        haps['ID_1'] = range(self._num_inds + 1)
        if self._ploidy == 1:
            haps['ID_2'] = "NA"
        else:
            haps['ID_2'] = haps['ID_1']
        haps['missing'] = np.repeat(0, self._num_inds + 1)

        # write to file
        haps.to_csv(out + "_argNeedle.sample", sep=' ', header=True, index=False)
