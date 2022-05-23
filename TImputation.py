#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 3 18:33:12 2022

@author: linkv
"""

import numpy as np
import os
import pandas as pd
import TTree as tt


def remove_monomorphic(trees):
    tables = trees.tables
    tables.sites.clear()
    tables.mutations.clear()
    n = trees.num_samples
    for tree in trees.trees():
        for site in tree.sites():
            visited = False
            for mutation in site.mutations:
                k = tree.num_samples(mutation.node)
                if k > 0 and k < n:
                    if not visited:
                        visited = True
                        site_id = tables.sites.add_row(
                            site.position, site.ancestral_state, metadata=site.metadata
                        )
                        tables.mutations.add_row(
                            site_id,
                            mutation.node,
                            mutation.derived_state,
                            parent=-1,
                            metadata=None,
                        )
    tables.compute_mutation_parents()
    return tables.tree_sequence()


class TImpute:

    def __init__(self):
        pass

    # TODO: getX and X2gen is creating the diploid genotype matrix. This should be produced by TTrees
    def getX(self, trees, idx):
        N = trees.num_samples
        print("idx", idx)
        M = idx
        X = np.zeros((N, M)).astype("int")
        i = 0
        num = 0
        for v in trees.variants():
            if i in idx:
                X[:, num] = v.genotypes
                num += 1
            i += 1
        return X

    def X2gen(self, X):
        n, m = X.shape
        maternals = np.array(range(0, n, 2))
        paternals = np.array(range(1, n, 2))
        X_dip = X[maternals, :] + X[paternals, :]

        tmp1 = (X_dip == 2).astype(int)
        tmp2 = (X_dip == 1).astype(int)
        tmp3 = (X_dip == 0).astype(int)

        n_gen = int(n * 3 / 2)
        buffer = np.zeros([n_gen, m]).astype(int)
        buffer[np.arange(0, n_gen, 3), :] = tmp1
        buffer[np.arange(1, n_gen, 3), :] = tmp2
        buffer[np.arange(2, n_gen, 3), :] = tmp3
        return buffer

    def gen2X(self, gen):
        n, m = gen.shape
        tmp1 = gen[np.arange(0, n, 3), :]
        tmp2 = gen[np.arange(1, n, 3), :]
        buffer = tmp1 * 2 + tmp2
        return buffer

    def run_impute_return_X(self, trees_sample, variants_ref, variants_sample, inds, inds_ref, genetic_map_file,
                            do_imputation, imputed_gen_file, out, logfile):
        """
        Impute
        @param genetic_map_file: str
        @param logfile:
        @param inds: TInds for sample
        @param inds_ref: TInds for reference panel
        @param variants_sample: TVariants for sample
        @param variants_ref: TVariants for reference panel
        @param trees_sample: tskit.treeSequence for sample (tree used for association testing)
        @param variants_ref: TVariants for reference panel
        @param out: str
        """
        if do_imputation:
            # write files in gen format
            name_imputation_output = out + "_imputed.gen"
            sample_gen_file = out + "_samples"
            reference_gen_file = out + "_reference"

            variants_sample.write_gen(sample_gen_file, inds, logfile)
            variants_ref.write_gen(reference_gen_file, inds_ref, logfile)

            # genetic map
            if genetic_map_file is None:
                genetic_map_file_name = variants_sample.write_genetic_map(out=out, logfile=logfile)
                logfile.info(
                    "- No genetic map file provided. Writing map with constant rate to " + genetic_map_file_name)
            else:
                genetic_map_file_name = genetic_map_file
                logfile.info("- Reading map with constant rate from " + genetic_map_file_name)

            logfile.info("- Starting imputation with impute2 for sample with " + str(variants_sample.num_typed)
                         + " typed variants using reference panel with " + str(variants_ref.num_typed) + " typed variants.")
            os.system(
                "./impute2 "
                + " -g_ref "
                + reference_gen_file + '.gen'
                + " -m "
                + genetic_map_file_name
                + " -g "
                + sample_gen_file + '.gen'
                + " -int 0 "
                + str(trees_sample.sequence_length)  # chromosome length
                + " -allow_large_regions "
                + " -o " + name_imputation_output
            )
        else:
            if imputed_gen_file is None:
                raise ValueError("When --do_imputation is set to False, the imputed genotypes must be provided with "
                                 "--imputed_gen_file parameter")
            logfile.info("- Assuming imputation was already run, reading imputed genotypes from " + imputed_gen_file)
            name_imputation_output = imputed_gen_file

        # read imputation results
        gen_imputed = pd.read_table(name_imputation_output, sep=" ", header=None).iloc[:, 5:]
        positions = np.array(pd.read_table(name_imputation_output, sep=" ", header=None).iloc[:, 2].values)
        gen_imputed = np.transpose(gen_imputed.values)
        X_imputed = self.gen2X(gen_imputed)

        # keep only polymorphic variants
        keep = np.logical_and(X_imputed.mean(axis=0) > 0, X_imputed.mean(axis=0) < 1)
        X_imputed = X_imputed[:, keep]
        logfile.info("- Done running impute2, imputed sample data set has " + str(X_imputed.shape[1])
                     + " variants, i.e. " + str(X_imputed.shape[1] - variants_sample.num_typed) + " more than before")

        return X_imputed, positions[keep]
