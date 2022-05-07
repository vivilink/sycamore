#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# TODO: being typed or not should be an option for all causal variants
import numpy as np


class Phenotypes:
    _y: np.ndarray
    _genetic_variance: float

    def __init__(self, variants, inds, logfile):
        self.num_inds = inds.num_inds
        self._y = np.zeros(self.num_inds)
        self._random_noise = np.zeros(self.num_inds)
        self._genetic_variance = -1.0
        self.betas = [0] * variants.number
        self.causal_variants = []
        self.causal_betas = []
        self.causal_power = []
        self.causal_trees = []
        self.causal_variant_indeces = []
        self.causal_tree_indeces = []
        self.filled = False

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y: np.ndarray):
        self._y = y

    @property
    def genetic_variance(self):
        return self._genetic_variance

    @genetic_variance.setter
    def genetic_variance(self, genetic_variance: float):
        self._genetic_variance = genetic_variance

    def simulate(self, args, r, logfile, variants_orig, inds, trees, plots_dir):
        """
        Simulate phenotypes
        :param args: TArgs
        :param r: TRandomGenerator
        :param logfile: IndentedLoggerAdapter
        :param variants_orig: TVariants
        :param inds: TInds
        :param trees: tskit.TreeSequence
        :param plots_dir: str
        :return:
        """
        # simulate trait architecture
        self.simulate_trait_architecture(args=args, r=r, logfile=logfile, variants_orig=variants_orig, inds=inds,
                                         trees=trees, plots_dir=plots_dir)

        # calculate genetic variance
        self._genetic_variance = np.var(self._y)
        logfile.info("- Simulated phenotypes with genetic variance component " + str(self._genetic_variance))

        # simulate noise
        if args.pty_sd_envNoise is None and args.pty_h_squared is not None:
            self._random_noise = self.simulate_env_noise_h(requested_hsquared=args.pty_h_squared, random=r)
            logfile.info("- Simulated random noise to result in h^2 " + str(args.pty_h_squared) +
                         ". The random variance component is thus " + str(np.var(self._random_noise)))
        elif args.pty_h_squared is None and args.pty_sd_envNoise is not None:
            self._random_noise = self.simulate_env_noise_sd(sd_random_noise=args.pty_sd_envNoise, random=r)
            logfile.info("- Simulated random noise with sd " + str(args.pty_sd_envNoise) +
                         ". The random variance component is thus " + str(np.var(self._random_noise)))
        else:
            raise ValueError("Must provide random noise distribution parameter. Either set noise sd with "
                             "'pty_sd_envNoise' or heritability with 'pty_h_squared'")

        self._y += self._random_noise

        # self.standardize(logfile)

        # write phenotypes to file
        self.write_to_file(variants_orig, inds, args.out, logfile)

        self.filled = True

    def simulate_trait_architecture(self, args, r, logfile, variants_orig, inds, trees, plots_dir):
        """
        Simulate phenotype's genetic architecture
        :param args: TArgs
        :param r: TRandomGenerator
        :param logfile: IndentedLoggerAdapter
        :param variants_orig: TVariants
        :param inds: TInds
        :param trees: tskit.TreeSequence
        :param plots_dir: str
        :return:
        """

        if args.pty_sim_method is None:
            raise ValueError("Must provide a phenotype simulation method with --pty_sim_method")

        if args.pty_sim_method == "null":
            logfile.info("- Simulating null phenotypes based only on random noise")
            self.simulate_null()

        if args.pty_sim_method == 'uniform':
            logfile.info(
                "- Simulating phenotypes based on uniformly chosen variants with prop_causal_mutations: " + str(
                    args.pty_prop_causal_mutations) + " and sd_beta_causal_mutations: " + str(
                    args.pty_sd_beta_causal_mutations))
            self.simulate_uniform(variants_orig, inds, prop_causal_mutations=args.pty_prop_causal_mutations,
                                  sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r, logfile=logfile)

        elif args.pty_sim_method == 'fixed':
            if args.pty_fixed_betas is None:
                raise ValueError("Must provide beta values provided for 'fixed' phenotype using '--pty_fixed_betas'")

            logfile.info("- Simulating phenotypes based on the following indeces: " + str(
                args.pty_fixed_variant_indeces) + " and the following betas: " + str(args.pty_fixed_betas))
            self.simulate_fixed(variants_orig, inds, args.pty_fixed_variant_indeces, args.pty_fixed_betas, logfile)

        elif args.pty_sim_method == 'singleTyped':
            if args.pty_fixed_betas is None:
                raise ValueError("Must provide beta values for 'singleTyped' phenotype using '--pty_fixed_betas'")
            if args.single_variant_af is None:
                raise ValueError(
                    "Must provide allele freq values for 'singleTyped' phenotype using '--single_variant_af'")

            fig, ax = plt.subplots(1, figsize=(30, 30))
            var_index, pos = variants_orig.find_variant(typed=True, freq=args.single_variant_af,
                                                        interval=args.single_variant_interval, subplot=ax,
                                                        random=r, logfile=logfile)
            fig.tight_layout()
            fig.set_size_inches(30, 30)
            fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logfile.info("- Simulating a phenotypes based on the following typed variant index: " + str(
                var_index) + " at position " + str(
                variants_orig.info['position'][var_index]) + " with allele freq " + str(
                variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(
                args.pty_fixed_betas))
            self.simulate_fixed(variants_orig, inds, [var_index], args.pty_fixed_betas, logfile)

        elif args.pty_sim_method == 'singleUntyped':
            if args.pty_fixed_betas is None:
                raise ValueError("Must provide beta values for phenotype 'singleUntyped' using '--pty_fixed_betas'")
            if args.single_variant_af is None:
                raise ValueError(
                    "Must provide allele freq values for 'singleTyped' phenotype using '--single_variant_af'")

            fig, ax = plt.subplots(1, figsize=(30, 30))
            var_index, pos = variants_orig.find_variant(typed=False, freq=args.single_variant_af,
                                                        interval=args.single_variant_interval, subplot=ax,
                                                        random=r, logfile=logfile)
            fig.tight_layout()
            fig.set_size_inches(30, 30)
            fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logfile.info("- Simulating a phenotypes based on the following untyped variant index: " + str(
                var_index) + " at position " + str(
                variants_orig.info['position'][var_index]) + " with allele freq " + str(
                variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(
                args.pty_fixed_betas))
            # to know which variants are untyped you need variants from simulated tree, not estimated tree
            if args.variants_file is None:
                raise ValueError(
                    "Must provide file with untyped variants to simulate phenotype with 'singleUntyped' model")
            self.simulate_fixed(variants_orig, inds, [var_index], args.pty_fixed_betas, logfile)

        elif args.pty_sim_method == 'oneTree':
            # TODO: This should be the tree from the original simulated trees, so that at different propTyped the
            #  experiments are comparable
            if args.causal_tree_pos is None:
                raise ValueError("Must provide causal tree position for phenotype 'oneTree' using '--causal_tree_pos'")
            if args.pty_sd_beta_causal_mutations is None:
                raise ValueError(
                    "Must provide effect size sd for phenotype 'oneTree' using --pty_sd_beta_causal_mutations")

            causal_tree = trees.at(args.causal_tree_pos)
            logfile.info("- Simulating phenotypes based on all variants of the tree covering postion " + str(
                args.causal_tree_pos))
            self.simulate_causal_region(variants_orig, inds, left_bound=causal_tree.interval.left,
                                        right_bound=causal_tree.interval.right,
                                        sd_beta_causal_mutations=args.pty_sd_beta_causal_mutations, random=r,
                                        logfile=logfile)

        elif args.pty_sim_method == 'allelicHetero':
            if args.allelic_hetero_file is None:
                raise ValueError("No instruction file provided for allelic heterogeneity simulation")
            if args.single_variant_af is not None or args.pty_fixed_betas is not None:
                raise ValueError(
                    "Provided allele frequency or beta value as well as instruction file for allelic heterogeneity. "
                    "Can accept only one type of instructions.")

            ah_info = pd.read_csv(args.allelic_hetero_file, delimiter="\t")
            variant_indeces = []
            fixed_betas = []
            sum_betas = 0

            logfile.info("- Searching for loci with requested allele frequencies.")
            logfile.add()

            # start plot
            # TODO: this plotting should not be done here but in a function instead
            fig, ax = plt.subplots(ah_info.shape[0], figsize=(30, 30))
            for index, row in ah_info.iterrows():
                # get allele freq
                f = -1
                if row['freq'] > 0.5:
                    f = 1 - row['freq']
                    logfile.warning("- Allele frequencies above 0.5 are not allowed. Transformed " + str(
                        row['freq']) + " to " + str(f) + ".")
                else:
                    f = row['freq']

                # get beta
                if not np.isnan(row['beta']) and not np.isnan(row['power']):
                    raise ValueError("Cannot fix power and beta value. One value must be set to 'NA'. Beta is " + str(
                        row['beta']) + " and power is " + str(row['power']))
                if np.isnan(row['beta']):
                    beta = np.sqrt(row['power'] / (f * (1 - f)))
                    # some betas should be negative
                    r_num = r.random.uniform(0, 1, 1)
                    if r_num < 0.5:
                        beta = -beta
                else:
                    beta = row['beta']

                var_index, pos = variants_orig.find_variant(typed=False, freq=f,
                                                            interval=[row["interval_start"], row["interval_end"]],
                                                            out=args.out, subplot=ax[index], random=r, logfile=logfile)
                variant_indeces.append(var_index)
                fixed_betas.append(beta)
                sum_betas += beta

                fig.tight_layout()
                fig.set_size_inches(30, 30)
                fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logfile.sub()

            logfile.info("- Simulating allelic heterogeneous phenotype with total beta " + str(sum_betas))

            logfile.info("- Simulating phenotypes:")
            logfile.add()
            self.simulate_fixed(variants=variants_orig, inds=inds, causal_variant_indeces=variant_indeces,
                                betas=fixed_betas, logfile=logfile)
            logfile.sub()

    @staticmethod
    def return_random_state(random):
        print(random.random.get_state()[1][0])
        print(random.random.uniform(0, 1, 1))

    def simulate_env_noise_sd(self, sd_random_noise, random):
        """
        Simulate random noise according to N(0, sd_random_noise)
        :param sd_random_noise: float
        :param random: TRandomGenerator
        :return noise: np.array
        """
        noise = random.random.normal(loc=0, scale=sd_random_noise, size=self.num_inds)
        return noise

    def simulate_env_noise_h(self, requested_hsquared, random):
        """
        Simulate random noise according to requested heritability and actual genetic variance
        :param requested_hsquared: float
        :param random: TRandomGenerator
        :return noise: np.array
        """
        V_G = self._genetic_variance
        V_E = V_G * (1 - requested_hsquared) / requested_hsquared
        sd_random_noise = np.sqrt(V_E)
        noise = random.random.normal(loc=0, scale=sd_random_noise, size=self.num_inds)
        return noise

    def simulate_null(self):
        """
        Do nothing, phenotype will only contain random effect
        :return:
        """
        return

    def simulate_fixed(self, variants, inds, causal_variant_indeces, betas, logfile):
        """
        Simulate phenotypes based on predefined causal variant positions and effects

        :param variants: TVariants
        :param inds: TInds
        :param causal_variant_indeces: list
        :param betas: list
        :param logfile:
        :return:
        """

        causal_variants = [variants.variants[i] for i in causal_variant_indeces]
        causal_pos = [variants.variants[i].site.position for i in causal_variant_indeces]

        if len(causal_variants) != len(betas):
            raise ValueError("must provide equal number of causal variants and betas to simulate fixed phenotype")

        for v, var in enumerate(causal_variants):
            # define beta
            self.betas[causal_variant_indeces[v]] = betas[v]

            # simulate phenotype
            if inds.ploidy == 1:
                self._y[var.genotypes == 1] += betas[v]
            else:
                genotypes = inds.get_diploid_genotypes(var.genotypes)
                self._y[genotypes == 1] += betas[v]
                self._y[genotypes == 2] += 2.0 * betas[v]

            # save causal position
            self.causal_variants.append(var)
            self.causal_betas.append(betas[v])
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            self.causal_power.append(betas[v] ** 2.0 * allele_freq * (1.0 - allele_freq))
            self.causal_variant_indeces.append(causal_variant_indeces[v])

            logfile.info("- Simulated causal variant at position " + str(causal_pos[v]) + " at index " + str(
                causal_variant_indeces[v]) + " with beta " + str(round(betas[v], 3)) + " and allele freq " + str(
                allele_freq) + " resulting in a power of " + str(
                round(betas[v] ** 2 * allele_freq * (1 - allele_freq), 3)))

    def simulate_uniform(self, variants, inds, prop_causal_mutations, sd_beta_causal_mutations, random, logfile,
                         mean_beta_causal_mutation=0.0):
        """
        Simulate phenotypes based on uniformly distributed causal variant positions and normally distributed effect sizes

        Parameters
        ----------
        variants : TVariants.py
            variants iterator from tskit.
        prop_causal_mutations : float
            proportion of variants that should be causal.
        sd_beta_causal_mutations : float
            sd of normal distribution for betas.
        mean_beta_causal_mutation : float

        Returns
        -------
        None.
        """

        # add phenotypic effect to mutations that are uniformly distributed
        for v, var in enumerate(variants.variants):

            # if variants.info["typed"] == True:

            r = random.random.uniform(0, 1, 1)

            if r < prop_causal_mutations:

                # define beta
                beta = random.random.normal(loc=mean_beta_causal_mutation, scale=sd_beta_causal_mutations, size=1)[0]
                self.betas[v] = beta

                # simulate phenotype
                if inds.ploidy == 1:
                    self._y[var.genotypes == 1] += beta
                    self._y[var.genotypes == 2] += 2.0 * beta
                else:
                    genotypes = inds.get_diploid_genotypes(var.genotypes)
                    self._y[genotypes == 1] += beta
                    self._y[genotypes == 2] += 2.0 * beta

                # save causal position
                self.causal_variants.append(var)
                self.causal_betas.append(beta)
                allele_freq = variants.info['allele_freq'][v]
                self.causal_power.append(beta ** 2 * allele_freq * (1 - allele_freq))
                self.causal_variant_indeces.append(v)

        logfile.info("- Simulated phenotypes based on " + str(
            len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True

    def simulate_causal_region(self, variants, inds, left_bound, right_bound, sd_beta_causal_mutations, random,
                               logfile):
        # add phenotypic effect to mutations that are uniformly distributed
        for v, var in enumerate(variants.variants):
            # is the variant in the tree
            if left_bound <= variants.info.loc[v]['position'] <= right_bound:

                # define beta
                beta = random.random.normal(loc=0, scale=sd_beta_causal_mutations, size=1)[0]
                self.betas[v] = beta

                # simulate phenotype
                if inds.ploidy == 1:
                    self._y[var.genotypes == 1] += beta
                    self._y[var.genotypes == 2] += 2.0 * beta
                else:
                    genotypes = inds.get_diploid_genotypes(var.genotypes)
                    self._y[genotypes == 1] += beta
                    self._y[genotypes == 2] += 2.0 * beta

                # save causal position
                self.causal_variants.append(var)
                self.causal_betas.append(beta)
                allele_freq = variants.info['allele_freq'][v]
                self.causal_power.append(beta ** 2 * allele_freq * (1 - allele_freq))
                self.causal_variant_indeces.append(v)

        logfile.info("- Simulated phenotypes based on " + str(
            len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True

    def find_causal_trees(self, ts_object):
        for v in self.causal_variants:
            causal_tree = ts_object.at(v.site.position)
            self.causal_tree_indeces.append(causal_tree.get_index())

    def diffs(self):
        cols = np.tile(self._y, (self.num_inds, 1))
        rows = cols.T
        buffer = cols - rows
        return np.abs(buffer)

    def standardize(self, logfile):
        logfile.info("- Standardizing phenotypes")
        self._y = (self._y - np.mean(self._y)) / np.std(self._y)

    def write_to_file_gcta_eGRM(self, inds, out, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format will match the binary output created with plinkFile R package.

        Returns
        -------
        None.

        """
        logfile.info("- Writing phenotype data in gcta format to '" + out + "_phenotypes.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.repeat(0, self.num_inds)
        tmp_pheno['2'] = inds.names
        tmp_pheno['3'] = self._y

        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)

    def write_to_file_gcta_scaled(self, out, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format
        will match the binary output created with egrm.

        Returns
        -------
        None.
        """
        logfile.info("- Writing phenotype data in gcta format to '" + out + "_phenotypes.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.arange(1, self.num_inds + 1)
        tmp_pheno['2'] = tmp_pheno['1']
        tmp_pheno['3'] = self._y

        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)

    def write_to_file(self, variants, inds, out, logfile):
        """
        Provide phenotypic variance partitioning information, and information related to each variant's phenotypic effect

        Returns
        -------
        None.

        """
        logfile.info("- Writing phenotype data '" + out + "_pheno_causal_vars.csv'")

        # results for each variant
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
        table['var_genotypic_from_betas'] = np.repeat(float(inds.ploidy), variants.number) * np.array(
            self.betas) * np.array(
            self.betas) * np.array(table['allele_freq']) * (np.repeat(1.0, variants.number) - np.array(
            table['allele_freq']))
        table['var_genotypic_empiric'] = np.repeat(self._genetic_variance, variants.number)
        table['var_random'] = np.repeat(np.var(self._random_noise), variants.number)
        table['var_phenotypic'] = np.var(self._y)

        table.to_csv(out + "_pheno_causal_vars.csv", index=False, header=True)
