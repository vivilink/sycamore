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

    def run_impute_return_X(self, trees_ref, trees_sample, variants_ref, variants_sample, inds, genetic_map_file, out, logfile):
        """
        Impute
        @param genetic_map_file: str
        @param logfile:
        @param inds: TInds
        @param variants_sample: TVariants for sample
        @param variants_ref: TVariants for reference panel
        @param trees_ref: tskit.treeSequence for reference panel
        @param trees_sample: tskit.treeSequence for sample
        @param variants_ref: TVariants for reference panel
        @param out: str
        """
        # trees_ref = remove_monomorphic(trees_ref.simplify(trees_ref.samples()))
        # trees_sample = remove_monomorphic(trees_sample.simplify(trees_sample.samples()))
        #
        # MAFs = np.array([v.genotypes.mean() for v in trees_sample.variants()])
        # loci = np.array([v.position for v in trees_sample.variants()])

        # MAFs = variants.info['allele_freq']
        # loci = variants.info['position']
        # M = loci.shape[0]

        # # observation
        # idx_5 = MAFs
        #
        # loci_ceil = np.ceil(loci)
        # print("loci_ceil", loci_ceil)
        # overlapped = np.insert(loci_ceil[:-1] == loci_ceil[1:], 1, False)
        # non_overlapped = np.logical_not(overlapped)
        # print("overlapped", overlapped)
        # print("non_overlapped", non_overlapped)

        # observable = variants.info['allele_freq'].values[variants.info['typed'] == True]
        # M_observable = observable.sum()
        #
        # Beta = 1
        # obs_ratio = 0.2
        # weights = np.multiply(np.power(MAFs, Beta), np.power(1 - MAFs, Beta))
        # weights = np.multiply(weights, observable)
        # weights /= weights.sum()
        #
        # M_obs = int(round(M_observable * obs_ratio))
        # obss = np.random.choice(
        #     np.where(observable)[0], M_obs, replace=False, p=weights[observable]
        # )
        # obss.sort()
        #
        # print("obss", obss)
        #
        # print("computing K_imputed", flush=True)
        # os.mkdir("impute")
        # os.chdir("impute")

        # write files in gen format
        name_imputation_output = out + "_imputed.gen"
        sample_gen_file = out + "_samples"
        reference_gen_file = out + "_reference"
        variants_sample.write_gen(sample_gen_file, inds, logfile)
        variants_ref.write_gen(reference_gen_file, inds, logfile)
        # X_obs = self.getX(trees_sample, variants_sample.num_typed)
        # gen_obs = self.X2gen(X_obs)
        # loci_obs = variants_sample.info['positions'].values[variants_sample.info['typed'] is True]
        #
        # haps_file = open(out + ".gen", "a")
        # i = 0
        # for idx, obs in enumerate(trees_sample.num_mutations):
        #     string = "chr snp" + str(obs + 1) + " " + str(loci_obs[idx]) + " A" + " T "
        #     string = string + " ".join(map(str, gen_obs[:, idx])) + "\n"
        #     bytes = haps_file.write(string)
        #     i += 1
        # haps_file.close()
        logfile.info("- Starting imputation with impute2 for sample with " + str(variants_sample.num_typed)
                     + " typed variants using reference panel with " + str(variants_ref.num_typed) + " typed variants.")
        # os.system(
        #     "~/git/argwas/impute2 "
        #     + " -g_ref "
        #     + reference_gen_file + '.gen'
        #     + " -m "
        #     + genetic_map_file
        #     + " -g "
        #     + sample_gen_file + '.gen'
        #     + " -int 0 "
        #     + str(trees_sample.sequence_length)  # chromosome length
        #     + " -allow_large_regions "
        #     + " -o " + name_imputation_output
        # )

        # read imputation results
        gen_imputed = pd.read_table(name_imputation_output, sep=" ", header=None).iloc[:, 5:]
        gen_imputed = np.transpose(gen_imputed.values)
        X_imputed = self.gen2X(gen_imputed)

        # keep only polymorphic variants
        keep = np.logical_and(X_imputed.mean(axis=0) > 0, X_imputed.mean(axis=0) < 1)
        X_imputed = X_imputed[:, keep]
        logfile.info("- Done running impute2, imputed sample data set has " + str(X_imputed.shape[1]) + " variants, i.e. "
                     + str(X_imputed.shape[1] - variants_sample.num_typed) + " more than before")

        # Z_imputed = X_imputed
        # Z_imputed -= Z_imputed.mean(axis=0)
        # Z_imputed /= Z_imputed.std(axis=0)
        # K_imputed = np.dot(Z_imputed, np.transpose(Z_imputed))
        # os.chdir("..")
        return X_imputed
