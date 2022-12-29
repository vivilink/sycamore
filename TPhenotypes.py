#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 16 17:52:46 2021

@author: linkv
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import norm

# TODO: being typed or not should be an option for all causal variants


class Phenotypes:
    def __init__(self):
        self._y = np.ndarray
        self._num_inds = -1
        self._genetic_variance = -1.0
        self.causal_variants = []
        self.causal_betas = []
        self.causal_power = []
        self.causal_trees = []
        self.causal_variant_indeces = []
        self.causal_tree_indeces = []
        self.causal_window_indeces = []
        self.filled = False

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y: np.ndarray):
        self._y = y

    def diffs(self):
        cols = np.tile(self._y, (self._num_inds, 1))
        rows = cols.T
        buffer = cols - rows
        return np.abs(buffer)

    def standardize(self, out, inds, logfile):
        logfile.info("- Standardizing phenotypes")
        self._y = (self._y - np.mean(self._y)) / np.std(self._y)

    def write_to_file_gcta_eGRM(self, inds, out, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format
        will match the binary output created with plinkFile R package.

        Returns
        -------
        None.

        """
        logfile.info("- Writing phenotype data in gcta format to '" + out + "_phenotypes.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.repeat(0, inds.num_inds)
        tmp_pheno['2'] = inds.names
        tmp_pheno['3'] = self._y

        indeces_to_remove = inds.get_indeces_inds_no_phenotype()

        if np.count_nonzero(np.isnan(self._y)) > 0 and len(indeces_to_remove) == 0:
            raise ValueError("There are phenotypes that are nan")

        # remove missing data
        if len(indeces_to_remove) > 0:
            tmp_pheno.drop(axis=0, index=indeces_to_remove, inplace=True)

        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)

    def write_to_file_gcta_scaled(self, out, inds, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format
        will match the binary output created with egrm.

        Returns
        -------
        None.
        """
        logfile.info("- Writing phenotype data in gcta format to '" + out + "_phenotypes.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.arange(1, inds.num_inds + 1)
        tmp_pheno['2'] = tmp_pheno['1']
        tmp_pheno['3'] = self._y

        # TODO: remove inds with missing phenotype data!

        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)


class PhenotypesBMI(Phenotypes):
    _sample_IDs: np.ndarray
    _pheno_df: pd.DataFrame

    def __init__(self, filename, inds, out, logfile):
        super().__init__()
        self.initialize_from_file(filename=filename, inds=inds, out=out, logfile=logfile)

    def initialize_from_file(self, filename, out, inds, logfile):
        if filename is None:
            raise ValueError("Provide file with BMI phenotype information using 'filename'")
        logfile.info("- Reading BMI phenotype information from " + filename)
        pheno_df = pd.read_csv(filename, names=["ID", "sex", "age", "BMI"])
        missing_in_phenotypes, added_in_phenotypes = self.find_missing_individuals(inds_tree=inds.names, inds_phenotype=pheno_df['ID'])
        logfile.info("- There are " + str(len(missing_in_phenotypes)) + " individuals missing from the phenotypes file "
                        "and " + str(len(added_in_phenotypes)) + " individuals added. Will add missing ones with NA and "
                        "remove added ones.")

        for i in missing_in_phenotypes:
            pheno_df.loc[len(pheno_df.index)] = [i, np.nan, np.nan, np.nan]

        for i in added_in_phenotypes:
            indexInd = pheno_df[(pheno_df['ID'] == i)].index
            pheno_df.drop(indexInd, inplace=True)

        pheno_df = self.sortPhenotypes(names_correct_order=inds.names, pheno_df=pheno_df)
        self._pheno_df = pheno_df

        self._num_inds = len(pheno_df['ID'])
        self._sample_IDs = np.array(pheno_df['ID'])

        # self._y = np.array(pheno_df['BMI'])
        # inform inds object about which inds have missing phenotypes
        # print(pheno_df)
        # print(pheno_df[pheno_df.isna().any(axis=1)])

        self.standardize(out=out, inds=inds, logfile=logfile)

    def set_missing_phenotype_status(self, inds):
        tmp = np.repeat(True, inds.num_inds)
        tmp[self._pheno_df.isna().any(axis=1)] = False
        if "outlier" in self._pheno_df.columns:
            tmp[self._pheno_df["outlier"] == True] = False
        inds.ind_has_phenotype = tmp

    @staticmethod
    def sortPhenotypes(names_correct_order, pheno_df):
        """
        Sort read-in phenotypes file according to sample names in sample file used to run Relate
        @param names_correct_order: sample names in order of file used to run Relate
        @param pheno_df: data frame from read-in phenotype file
        @return: sorted pheno_df
        """
        df_mapping = pd.DataFrame({'names_correct_order': names_correct_order, })
        sort_mapping = df_mapping.reset_index().set_index('names_correct_order')
        pheno_df['ind_num'] = pheno_df['ID'].map(sort_mapping['index'])
        pheno_df = pheno_df.sort_values('ind_num')

        return pheno_df

    @staticmethod
    def find_missing_individuals(inds_tree, inds_phenotype):
        set_tree = set(inds_tree)
        set_phenotype = set(inds_phenotype)
        missing_in_phenotypes = list(sorted(set_tree - set_phenotype))
        added_in_phenotypes = list(sorted(set_phenotype - set_tree))

        return missing_in_phenotypes, added_in_phenotypes

    @property
    def sample_IDs(self):
        return self._sample_IDs

    def regression_on_age_sex_find_outliers(self, inds, out):
        """
        Normalize residuals according to McCaw et al. (2019) Operating characteristics of the rank-based inversenormal
        transformation for quantitative trait analysisin genome-wide association studies :return:
        """
        self._pheno_df['residual'] = np.nan
        self._pheno_df['outlier'] = False
        self._pheno_df['rank_inv_transform'] = np.nan

        for sex in [1.0, 2.0]:
            # perform linear regression and save residuals
            results = smf.ols('BMI ~ sex + age +' + 'I(age**2)',
                              data=self._pheno_df[self._pheno_df['sex'] == sex]).fit()
            residuals = results.resid
            self._pheno_df.loc[
                (self._pheno_df['sex'] == sex) & (self._pheno_df['BMI'].notnull()), 'residual'] = residuals

            # define outliers
            sd = np.std(residuals)
            mn = np.mean(residuals)
            self._pheno_df.loc[(self._pheno_df['sex'] == sex) & (self._pheno_df['residual'] > (mn + 6.0 * sd))
                               | ((self._pheno_df['sex'] == sex) & (self._pheno_df['residual'] < (mn - 6.0 * sd)))
                                , "outlier"] = True

            # do rank-based inversement transformation
            # TODO: should I exclude outliers from this?
            ecdf = sm.distributions.ECDF(residuals)
            quants = ecdf(residuals)
            self._pheno_df.loc[(self._pheno_df['sex'] == sex) & (self._pheno_df['BMI'].notnull()), 'rank_inv_transform'] = norm.ppf(quants)

        # write to file and set necessary parameters
        self.set_missing_phenotype_status(inds=inds)
        self._pheno_df.to_csv(out + "_standardized_pheno_df.csv")
        self._y = np.array(self._pheno_df['rank_inv_transform'])

    def standardize(self, out, inds, logfile):
        """
        @param logfile:
        @return:
        """
        logfile.info("- standardizing phenotypes with indirect rank inverse transformation")
        self.regression_on_age_sex_find_outliers(out=out, inds=inds)

    def write_to_file_gcta_eGRM(self, inds, out, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format
        will match the binary output created with plinkFile R package.

        Returns
        -------
        None.

        """
        print("in write to file gcta eGRM of BMI phenotypes")

        self.standardize(out=out, inds=inds, logfile=logfile)

        super().write_to_file_gcta_eGRM(inds, out, logfile)


class PhenotypesSimulated(Phenotypes):
    _genetic_variance: float

    def __init__(self, variants, num_inds):
        super().__init__()

        self._random_noise = np.zeros(num_inds)
        self.betas = [0] * variants.number
        self._y = np.empty(num_inds)

    @property
    def genetic_variance(self):
        return self._genetic_variance

    @genetic_variance.setter
    def genetic_variance(self, genetic_variance: float):
        self._genetic_variance = genetic_variance

    @staticmethod
    def set_missing_phenotype_status(inds):
        tmp = np.repeat(True, inds.num_inds)
        inds.ind_has_phenotype = tmp

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

        print("_y before adding noise 1", self._y)

        # calculate genetic variance
        self._genetic_variance = float(np.var(self._y))
        logfile.info("- Simulated phenotypes with genetic variance " + str(self._genetic_variance))

        if args.pty_h_squared is not None and args.pty_h_squared > 0 and self._genetic_variance == 0:
            raise ValueError("Genetic variance is zero and heritability is not zero. Make sure you are simulating "
                             "causal variants.")

        # simulate noise
        if args.pty_sd_envNoise is None and args.pty_h_squared is not None:
            self._random_noise = self.simulate_env_noise_h(requested_hsquared=args.pty_h_squared, inds=inds, random=r)
            logfile.info("- Simulated random noise to result in h^2 " + str(args.pty_h_squared) +
                         ". The variance of the random noise is thus " + str(np.var(self._random_noise)))
        elif args.pty_h_squared is None and args.pty_sd_envNoise is not None:
            self._random_noise = self.simulate_env_noise_sd(sd_random_noise=args.pty_sd_envNoise, inds=inds, random=r)
            logfile.info("- Simulated random noise with sd " + str(args.pty_sd_envNoise) +
                         ". The variance of the random noise is thus " + str(np.var(self._random_noise)))
        else:
            raise ValueError("Must provide random noise distribution parameter. Either set noise sd with "
                             "'pty_sd_envNoise' or heritability with 'pty_h_squared'")

        print("_y before adding noise 2", self._y)

        self._y += self._random_noise

        # self.standardize(logfile)

        # write phenotypes to file
        self.write_sim_params_to_file(variants_orig, inds, args.out, logfile)

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

        elif args.pty_sim_method == "null":
            logfile.info("- Simulating null phenotypes based only on random noise")
            self.simulate_null()

        elif args.pty_sim_method == 'uniform':
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
            if args.single_variant_interval[0] < 0:
                raise ValueError("single_variant_interval start cannot be negative")
            if args.single_variant_interval[1] > trees.sequence_length:
                raise ValueError("single_variant_interval end cannot be larger than tree sequence length, which is "
                                 + str(trees.sequence_length))

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
                    "Must provide allele freq values for 'singleUntyped' phenotype using '--single_variant_af'")
            if args.single_variant_interval[1] > trees.sequence_length:
                raise ValueError("single_variant_interval end cannot be larger than tree sequence length, which is "
                                 + str(trees.sequence_length))

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
                                        causal_mutations_effect_size_def=args.pty_sd_beta_causal_mutations,
                                        random=r,
                                        local_heritability=args.pty_h_squared,
                                        min_allele_freq_causal=args.min_allele_freq_causal,
                                        prop_causal_mutations=args.pty_prop_causal_mutations,
                                        max_allele_freq_causal=args.max_allele_freq_causal,
                                        allow_typed_causal_variants=args.allow_typed_causal_variants,
                                        logfile=logfile)

        elif args.pty_sim_method == "oneRegion":
            if args.causal_region_coordinates is None or len(args.causal_region_coordinates) != 2:
                raise ValueError("Must provide start and end coordinates of causal window using "
                                 "'causal_region_coordinates'")
            if args.pty_prop_causal_mutations is None or 1 < args.pty_prop_causal_mutations < 0:
                raise ValueError("Must provide proportion of causal variants between 0 and 1 using "
                                 "'pty_prop_causal_mutations'")
            if args.pty_sd_beta_causal_mutations is None:
                raise ValueError("Must provide sd for causal mutation effect size using "
                                 "'pty_sd_beta_causal_mutations'")
            logfile.info("- Simulating phenotypes based on variants within region with coordinates " + str(
                args.causal_region_coordinates))

            self.simulate_causal_region(variants=variants_orig,
                                        inds=inds,
                                        left_bound=args.causal_region_coordinates[0],
                                        right_bound=args.causal_region_coordinates[1],
                                        causal_mutations_effect_size_def=args.pty_sd_beta_causal_mutations,
                                        local_heritability=args.pty_h_squared,
                                        random=r,
                                        prop_causal_mutations=args.pty_prop_causal_mutations,
                                        min_allele_freq_causal=args.min_allele_freq_causal,
                                        max_allele_freq_causal=args.max_allele_freq_causal,
                                        allow_typed_causal_variants=args.allow_typed_causal_variants,
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

    def simulate_env_noise_sd(self, sd_random_noise, random, inds):
        """
        Simulate random noise according to N(0, sd_random_noise)
        :param sd_random_noise: float
        :param random: TRandomGenerator
        :return noise: np.array
        """
        noise = random.random.normal(loc=0, scale=sd_random_noise, size=inds.num_inds)

        print("nas in noise simulate_env_noise_sd", np.count_nonzero(np.isnan(noise)))

        return noise

    def simulate_env_noise_h(self, requested_hsquared, inds, random):
        """
        Simulate random noise according to requested heritability and actual genetic variance
        :param requested_hsquared: float
        :param random: TRandomGenerator
        :return noise: np.array
        """
        if requested_hsquared == 0.0:
            V_E = 1.0
            sd_random_noise = np.sqrt(V_E)
            noise = random.random.normal(loc=0, scale=sd_random_noise, size=inds.num_inds)

        else:
            V_G = self._genetic_variance
            if V_G == 0:
                raise ValueError("Genetic variance is zero")
            V_E = V_G * (1 - requested_hsquared) / requested_hsquared
            sd_random_noise = np.sqrt(V_E)
            noise = random.random.normal(loc=0, scale=sd_random_noise, size=inds.num_inds)

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
        print(self._y)

    def simulate_causal_region(self, variants, inds, left_bound, right_bound, causal_mutations_effect_size_def,
                               local_heritability, prop_causal_mutations, random, min_allele_freq_causal,
                               max_allele_freq_causal, logfile, allow_typed_causal_variants):
        """
        Simulate causal effect sizes for variants with
        @param allow_typed_causal_variants: If true, only also typed variants can be causal
        @param prop_causal_mutations: float [0,1]
        @param local_heritability: float [0,1]
        @param variants: TVariants
        @param inds: TInds
        @param left_bound: float. Left bound of causal region (included)
        @param right_bound: Right bound of causal region (included)
        @param causal_mutations_effect_size_def: str, Std. dev. for betas of causal mutations . If it can be converted to a
                                        float, betas will sampled from N(0, pty_sd_beta_causal_mutations). If set to
                                        'standardized', betas will be sampled from
                                        N(0, [2 * f * (1 - f)]^{-0.5} * h2g / p),
                                        where h2g is the heritability of the trait and p is the number of causal SNPs.
        @param random: TRandom
        @param min_allele_freq_causal: float. Only variants with freq higher than this can be causal
        @param max_allele_freq_causal: Only variants with freq lower than this can be causal
        @param logfile:
        @return:
        """

        # define causal variants
        info_window = variants.info.loc[(left_bound <= variants.info['position'])
                                        & (variants.info['position'] < right_bound)
                                        & (min_allele_freq_causal <= variants.info['allele_freq'])
                                        & (variants.info['allele_freq'] <= max_allele_freq_causal)]

        variants.info.loc[(left_bound <= variants.info['position'])
                          & (variants.info['position'] < right_bound), "causal_region"] = "TRUE"

        # remove typed variants
        if not allow_typed_causal_variants:
            info_window = info_window.loc[info_window['typed'] == False]

        if len(info_window.index) < 1:
            raise ValueError("Found no variants in causal region. Are you restricting to only typed or untyped "
                             "variants?")

        tmp = random.random.uniform(0, 1, len(info_window.index))
        causal = np.zeros(len(info_window.index))
        causal[tmp < prop_causal_mutations] = 1
        num_causal_vars = np.count_nonzero(causal == 1)

        if num_causal_vars < 1:
            logfile.info("- WARNING: Number of causal variants is 0!")

        # add phenotypic effect to mutations that are uniformly distributed
        for v_i, v in enumerate(info_window['var_index']):
            if causal[v_i] == 1:
                # get beta
                sd = None
                # get effect size sd
                try:
                    sd = float(causal_mutations_effect_size_def)
                except ValueError:
                    causal_mutations_effect_size_def = causal_mutations_effect_size_def
                    if causal_mutations_effect_size_def != "standardized":
                        raise ValueError("'sd_beta_causal_mutations' is set to unknown option. Must be a float or "
                                         "'standardized'")
                    f = variants.info['allele_freq'][v]
                    sd = (local_heritability / num_causal_vars) / np.sqrt(2 * f * (1 - f))

                beta = self.get_beta_normal(random, sd)
                self.betas[v] = beta

                # simulate phenotype
                if inds.ploidy == 1:
                    self._y[variants.variants[v].genotypes == 1] += beta
                else:
                    genotypes = inds.get_diploid_genotypes(variants.variants[v].genotypes)
                    self._y[genotypes == 1] += beta
                    self._y[genotypes == 2] += 2.0 * beta

                # save causal position
                self.causal_variants.append(variants.variants[v])
                self.causal_betas.append(beta)
                allele_freq = variants.info['allele_freq'][v]
                self.causal_power.append(beta ** 2 * allele_freq * (1 - allele_freq))
                self.causal_variant_indeces.append(v)

        logfile.info("- Simulated phenotypes based on " + str(
            len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True

    @staticmethod
    def get_beta_normal(random, sd_beta_causal_mutations):
        """
        Draw effect size from normal distribution

        @param sd_beta_causal_mutations: float
        @param random: TRandom
        @return: float
        """
        beta = random.random.normal(loc=0, scale=sd_beta_causal_mutations, size=1)[0]
        return beta

    def find_causal_trees(self, ts_object):
        for v in self.causal_variants:
            causal_tree = ts_object.at(v.site.position)
            self.causal_tree_indeces.append(causal_tree.get_index())
        # remove duplicates
        self.causal_tree_indeces = list(set(self.causal_tree_indeces))

    def find_causal_windows(self, window_ends, window_starts):
        for v in self.causal_variants:
            for w in range(len(window_ends)):
                if window_starts[w] <= v.site.position < window_ends[w]:
                    self.causal_window_indeces.append(w)

    def write_sim_params_to_file(self, variants, inds, out, logfile):
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
        table['causal_region'] = variants.info['causal_region']
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

        print("nas in write_sim_params_to_file", np.count_nonzero(np.isnan(self._y)))

        table.to_csv(out + "_pheno_causal_vars.csv", index=False, header=True)

