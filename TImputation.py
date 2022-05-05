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

    def __init__(self, trees_ref, trees_sample, variants_ref, variants_sample, inds, genetic_map_file, out, logfile):
        self.run_impute(trees_ref=trees_ref, trees_sample=trees_sample, variants_ref=variants_ref,
                        variants_sample=variants_sample, inds=inds, genetic_map_file=genetic_map_file,
                        out=out, logfile=logfile)

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

    def run_impute(self, trees_ref, trees_sample, variants_ref, variants_sample, inds, genetic_map_file, out, logfile):
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
        trees_ref = remove_monomorphic(trees_ref.simplify(trees_ref.samples()))
        trees_sample = remove_monomorphic(trees_sample.simplify(trees_sample.samples()))

        MAFs = np.array([v.genotypes.mean() for v in trees_sample.variants()])
        loci = np.array([v.position for v in trees_sample.variants()])

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

        # write sample files
        sample_gen_file = out + "_samples"
        variants_sample.write_gen(sample_gen_file, inds, logfile)
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

        # write reference files
        reference_gen_file = out + "_reference"
        variants_ref.write_gen(out + "_reference", inds, logfile)
        inds.write_shapeit2(out, logfile)

        # X_ref = self.getX(trees_ref, np.arange(trees_ref.num_mutations))
        # gen_ref = self.X2gen(X_ref)
        # loci_ref = np.array([v.position for v in trees_ref.variants()])
        # loci_ref = np.ceil(loci_ref).astype(int)
        #
        # haps_file = open(out + "_ref.gen", "a")
        # i = 0
        # for idx, obs in enumerate(np.arange(trees_ref.num_mutations)):
        #     string = "1 refsnp" + str(obs + 1) + " " + str(loci_ref[idx]) + " A" + " T "
        #     string = string + " ".join(map(str, gen_ref[:, idx])) + "\n"
        #     bytes = haps_file.write(string)
        #     i += 1
        # haps_file.close()

        # map_file = open(out + ".map", "a")
        # map_file.write("pos COMBINED_rate Genetic_Map\n")
        # for idx, obs in enumerate(variants_ref.info['allele_freq']):
        #     string = str(loci_obs[idx]) + " " + str(1) + " "
        #     string = string + str(loci_obs[idx] / 1000000) + "\n"
        #     bytes = map_file.write(string)
        # map_file.close()

        os.system(
            "~/git/argwas/impute2 "
            + " -g_ref "
            + reference_gen_file + '.gen'
            + " -m "
            + genetic_map_file
            + " -g "
            + sample_gen_file + '.gen'
            + " -int 0 "
            + str(trees_sample.sequence_length)  # chromosome length
            + " -allow_large_regions "
            + " -o out.gen"
        )

        gen_imputed = pd.read_table("out.gen", sep=" ", header=None).iloc[:, 5:]
        gen_imputed = np.transpose(gen_imputed.values)
        X_imputed = self.gen2X(gen_imputed)

        keep = np.logical_and(X_imputed.mean(axis=0) > 0, X_imputed.mean(axis=0) < 1)
        X_imputed = X_imputed[:, keep]

        Z_imputed = X_imputed
        Z_imputed -= Z_imputed.mean(axis=0)
        Z_imputed /= Z_imputed.std(axis=0)
        K_imputed = np.dot(Z_imputed, np.transpose(Z_imputed))

        os.chdir("..")
