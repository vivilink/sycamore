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
from scipy.stats import norm
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.stats import norm
import TVariants as tvar
import tskit
import TRandomGenerator as rg
import TParameters as tparams
import TIndividuals as tind
from python_log_indenter import IndentedLoggerAdapter
import TTree as tt


def transform_to_binary(pheno_file: str, population_disease_prevalence: list[float], out: str,
                        logger: IndentedLoggerAdapter, pheno_column=2):
    """
    Transform phenotypes of a sample in liability scale to binary scale by using the population_disease_prevalence as
    cutoff. It expects a .phen file with three columns. 1: family ID, 2: sample ID, 3: phenotype on liability scale (
    continuous). It produces a file where the liability phenotypes are replaced by binary phenotype, defined based on
    population_disease_prevalence cutoff.

     :param pheno_column: column in pheno file with liability phenotypes
     :param pheno_file:
     :param population_disease_prevalence:
     :param out:
     :param logger:

    """
    if pheno_file is None:
        raise ValueError("Provide pheno_file file with 'pheno_file'")
    if population_disease_prevalence is None:
        raise ValueError("Provide expected population prevalence with 'population_disease_prevalence'")

    # standardize the phenotype
    pheno_file = pd.read_csv(pheno_file, delim_whitespace=True, header=None)
    pheno_standardized = pheno_file.iloc[:, pheno_column]
    pheno_standardized = (pheno_standardized - np.mean(pheno_standardized)) / np.std(pheno_standardized)

    # find cutoffs for standardized phenotypes
    qnorm_values = [norm.ppf(1 - i) for i in population_disease_prevalence]
    # assign cases
    for pp in qnorm_values:
        tmp = pheno_standardized
        tmp.loc[tmp < pp] = 0
        tmp.loc[tmp >= pp] = 1
        pheno_file[pheno_file.shape[1]] = tmp.astype(int)

    logger.info("- Writing population's disease status to '" + out + "_disease_status.phen'")
    pheno_file.to_csv(out + "_disease_status.phen", sep=' ', index=False, header=False)


def ascertain_sample(pheno_file: str, sample_size: int, sample_prevalence: float, out: str, random: rg,
                     logger: IndentedLoggerAdapter, pheno_column:int=3):
    """
    Creates ascertained sample from the population file of phenotypes (e.g. .phen file produced by task
    'transformToBinaryPhenotype'). It samples cases and controls according to satisfy the sample_prevalence parameter.

    :param pheno_column: column in pheno file with binary phenotypes
    :param pheno_file:
    :param sample_size:
    :param sample_prevalence:
    :param out:
    :param random:
    :param logger:
    :return:
    """
    if pheno_file is None:
        raise ValueError("Provide pheno_file file with 'pheno_file'")
    if sample_size is None:
        raise ValueError("Provide expected sample size with 'sample_size'")
    if sample_prevalence is None:
        raise ValueError("Provide expected sample prevalence with 'sample_prevalence'")

    # read file and find population size
    pheno_file = pd.read_csv(pheno_file, delim_whitespace=True, header=None)
    # pheno_file.columns = ['1', '2', '3']
    population_size = pheno_file.shape[0]
    if population_size < sample_size:
        raise ValueError(
            "Sample size cannot be larger than number of phenotypes in phen file (" + str(population_size) + ")")

    # determine number of cases and controls
    num_cases_required = round(np.ceil(sample_prevalence * sample_size))
    num_controls_required = sample_size - num_cases_required
    logger.info(
        "- Keeping " + str(num_cases_required) + " cases and " + str(num_controls_required) + " controls from a total of " +
        str(population_size) + " individuals")
    if num_controls_required <= 0:
        raise ValueError("Sample size needs to be larger than number of cases")

    # remove surplus cases
    cases_rows = pheno_file[pheno_file.iloc[:, pheno_column] == 1].index.tolist()
    num_cases_available = len(cases_rows)

    if num_cases_available < num_cases_required:
        raise ValueError("Cannot satisfy sample_prevalence of " + str(num_cases_required) + ". Only "
                         + str(num_cases_available) + " cases are available.")

    if num_cases_available > num_cases_required:
        random_rows = np.random.choice(cases_rows, num_cases_available - num_cases_required, replace=False)
        for row in random_rows:
            pheno_file.iloc[row, pheno_column] = -9

    # remove surplus controls
    controls_rows = pheno_file[pheno_file.iloc[:, pheno_column] == 0].index.tolist()
    num_controls_available = len(controls_rows)

    if num_controls_available < num_controls_required:
        raise ValueError("Cannot satisfy requested number " + str(num_controls_required) + " of controls. Only "
                         + str(num_controls_available) + " controls are available.")

    if num_controls_available > num_controls_required:
        random_rows = np.random.choice(controls_rows, num_controls_available - num_controls_required, replace=False)
        for row in random_rows:
            pheno_file.iloc[row, pheno_column] = -9

    # keep only lines without missing phenotype
    pheno_file = pheno_file[pheno_file.iloc[:, pheno_column] != -9]

    # write
    logger.info(
        "- Writing ascertained sample's disease status to '" + out + "_disease_status_ascertained.phen'")
    pheno_file.to_csv(out + "_disease_status_ascertained.phen", sep=' ', index=False, header=False)


# TODO: being typed or not should be an option for all causal variants
def make_phenotypes(args: tparams, trees: tskit.trees, sample_ids, inds: tind, plots_dir: str, random: rg,
                    logfile: IndentedLoggerAdapter):
    if args.simulate_phenotypes or args.task == "simulatePhenotypes":
        if args.tree_file_simulated is None:
            raise ValueError("To simulate phenotypes based on untyped variants the simulated trees need to be "
                             "provided with 'tree_file_simulated'.")

        logfile.info("- Reading simulated tree used for simulating phenotypes from " + args.tree_file_simulated)
        trees_object = tt.TTrees(tree_file=args.tree_file,
                                 trees_interval=args.trees_interval,
                                 trees_interval_start=args.trees_interval_start,
                                 trees_interval_end=args.trees_interval_end,
                                 skip_first_tree=args.skip_first_tree,
                                 logfile=logfile)
        trees_orig = trees_object.trees

        if trees_orig.num_samples != trees.num_samples:
            raise ValueError("The trees provided with params 'tree_file_simulated' and 'tree_file' must have the same "
                             "number of samples")

        # variants_orig are used to simulate phenotypes. They need to be consistent with original tree and the typed
        # status that might have been defined earlier with a variants file. The causal mutation should not be affected
        # by a freq filter

        variants_orig = tvar.TVariantsFiltered(tskit_object=trees_orig,
                                               samp_ids=sample_ids,
                                               min_allele_freq=0,
                                               max_allele_freq=1,
                                               prop_typed_variants=1,
                                               pos_float=args.pos_float,
                                               random=random,
                                               logfile=logfile,
                                               filtered_variants_file=args.variants_file)

        logfile.info("- Phenotypes:")
        logfile.add()

        pheno = PhenotypesSimulated(variants=variants_orig, num_inds=inds.num_inds)

        pheno.simulate(args=args, r=random, logfile=logfile, variants_orig=variants_orig, trees=trees, inds=inds,
                       plots_dir=plots_dir, samp_ids=sample_ids, trees_orig=trees_orig)
        logfile.sub()

    else:
        if args.pheno_file_BMI:
            pheno = PhenotypesBMI()
            pheno.initialize_from_file(filename=args.pheno_file_BMI, inds=inds, out=args.out, logfile=logfile)
            # TODO: maybe restrict tree to inds for which we have phenotypes here
        else:
            pheno = Phenotypes()
            pheno.initialize_from_file(filename=args.pheno_file, inds=inds, out=args.out, logfile=logfile,
                                       num_phenotypes=args.num_phenotypes,
                                       phenotype_number_of_interest=args.phenotype_number_of_interest)

    return pheno


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


def find_missing_individuals(inds_tree, inds_phenotype):
    set_tree = set(inds_tree)
    set_phenotype = set(inds_phenotype)
    missing_in_phenotypes = list(sorted(set_tree - set_phenotype))
    added_in_phenotypes = list(sorted(set_phenotype - set_tree))

    return missing_in_phenotypes, added_in_phenotypes


class Phenotypes:
    _pheno_df: pd.DataFrame

    def __init__(self):
        self._sample_IDs = None
        self._y = np.ndarray
        self._disease_status = np.ndarray
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
        self._pheno_df = None  # this is needed when the phentypes are not simulated but initialized from file (there might be missing individuals)

    @property
    def y(self):
        return self._y

    @y.setter
    def y(self, y: np.ndarray):
        self._y = y

    @property
    def num_inds(self):
        return self._num_inds

    @num_inds.setter
    def num_inds(self, num_inds: int):
        self._num_inds = num_inds

    def set_missing_phenotype_status(self, inds):
        if self._pheno_df is None:
            tmp = np.repeat(True, inds.num_inds)
            inds.ind_has_phenotype = tmp
        else:
            tmp = np.repeat(True, inds.num_inds)
            tmp[self._pheno_df.isna().any(axis=1)] = False
            if "outlier" in self._pheno_df.columns:
                tmp[self._pheno_df["outlier"] == True] = False
            inds.ind_has_phenotype = tmp

    # @staticmethod
    # def set_missing_phenotype_status(inds):
    #     tmp = np.repeat(True, inds.num_inds)
    #     inds.ind_has_phenotype = tmp

    def diffs(self):
        cols = np.tile(self._y, (self._num_inds, 1))
        rows = cols.T
        buffer = cols - rows
        return np.abs(buffer)

    def standardize(self, out, inds, logfile):
        logfile.info("- Standardizing phenotypes")
        self._y = (self._y - np.mean(self._y)) / np.std(self._y)

    def write_to_file_fam(self, inds, out):
        """
        Write fam file (plink file format) with missing phenotypes (-9)
        :param inds:
        :param out:
        :param logfile:
        :return:
        """
        # logfile.info("- Writing phenotype data in fam format to '" + out + ".fam")
        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.repeat(0, inds.num_inds)
        tmp_pheno['2'] = inds.names
        tmp_pheno['3'] = np.repeat(0, inds.num_inds)
        tmp_pheno['4'] = np.repeat(0, inds.num_inds)
        tmp_pheno['5'] = np.repeat(-9, inds.num_inds)
        tmp_pheno['6'] = np.repeat(-9, inds.num_inds)

        self.set_missing_phenotype_status(inds=inds)
        indeces_to_remove = inds.get_indeces_inds_no_phenotype()

        if np.count_nonzero(np.isnan(self._y)) > 0 and len(indeces_to_remove) == 0:
            raise ValueError("There are phenotypes that are nan")

        # remove missing data
        if len(indeces_to_remove) > 0:
            tmp_pheno.drop(axis=0, index=indeces_to_remove, inplace=True)

        tmp_pheno.to_csv(out + ".fam", sep=' ', index=False, header=False)

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

        self.set_missing_phenotype_status(inds=inds)
        indeces_to_remove = inds.get_indeces_inds_no_phenotype()

        if np.count_nonzero(np.isnan(self._y)) > 0 and len(indeces_to_remove) == 0:
            raise ValueError("There are phenotypes that are nan")

        # remove missing data
        if len(indeces_to_remove) > 0:
            tmp_pheno.drop(axis=0, index=indeces_to_remove, inplace=True)

        tmp_pheno.to_csv(out + "_phenotypes.phen", sep=' ', index=False, header=False)

    def write_to_file_gcta_eGRM_disease_status(self, inds, out, logfile):
        """
        Write phenotypes to file in gtca format (first column=family, second=ind id, third=pheno value). This format
        will match the binary output created with plinkFile R package.

        Returns
        -------
        None.

        """
        logfile.info("- Writing disease status data in gcta format to '" + out + "_disease_status.phen'")

        tmp_pheno = pd.DataFrame()
        tmp_pheno['1'] = np.repeat(0, inds.num_inds)
        tmp_pheno['2'] = inds.names
        tmp_pheno['3'] = self._disease_status

        self.set_missing_phenotype_status(inds=inds)
        indeces_to_remove = inds.get_indeces_inds_no_phenotype()

        if np.count_nonzero(np.isnan(self._disease_status)) > 0 and len(indeces_to_remove) == 0:
            raise ValueError("There are phenotypes that are nan")

        # remove missing data
        if len(indeces_to_remove) > 0:
            tmp_pheno.drop(axis=0, index=indeces_to_remove, inplace=True)

        tmp_pheno.to_csv(out + "_disease_status.phen", sep=' ', index=False, header=False)

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

    def initialize_from_file(self, filename, out, inds, logfile, num_phenotypes=1, phenotype_number_of_interest=1):
        if filename is None:
            raise ValueError("Provide file with phenotype information in gcta .phen format using 'pheno_file' or "
                             "simulate phenotypes using 'simulate_phenotypes'")
        logfile.info("- Reading phenotype information from " + filename)
        header = ["0", "ID"]
        for i in range(1, num_phenotypes + 1):
            header.append("phenotype" + str(i))
        pheno_df = pd.read_csv(filename, names=header, sep=' ')

        missing_in_phenotypes, added_in_phenotypes = find_missing_individuals(inds_tree=inds.names,
                                                                              inds_phenotype=pheno_df['ID'])
        if len(missing_in_phenotypes) == len(inds.names):
            raise ValueError("There are no individual names in common between the ARG and phenotype file")
        logfile.info("- There are " + str(len(missing_in_phenotypes)) + " individuals missing from the phenotypes file "
                                                                        "and " + str(
            len(added_in_phenotypes)) + " individuals added. Will add missing ones with NA and "
                                        "remove added ones.")

        for i in missing_in_phenotypes:
            pheno_df.loc[len(pheno_df.index)] = [0, i, np.nan]

        for i in added_in_phenotypes:
            indexInd = pheno_df[(pheno_df['ID'] == i)].index
            pheno_df.drop(indexInd, inplace=True)

        self.set_missing_phenotype_status(inds=inds)

        pheno_df = sortPhenotypes(names_correct_order=inds.names, pheno_df=pheno_df)
        self._pheno_df = pheno_df

        self._num_inds = len(pheno_df['ID'])
        self._sample_IDs = np.array(pheno_df['ID'])
        self._y = np.array(pheno_df['phenotype' + str(phenotype_number_of_interest)])

    def scale(self):
        self._y = (self._y - np.mean(self._y)) / np.std(self._y)

    def add_disease_status(self, prevalence, logfile):
        """
        Define binary phenotype based on liability score. Individuals with a liability score above threshold given by
        prevalence will have disease
        :param prevalence:
        :return:
        """
        logfile.info("- Scaling phenotypes to add disease status")
        self.scale()
        liability_cutoff = norm.ppf(1 - prevalence)
        logfile.info("- The liability cutoff is set to " + str(liability_cutoff))
        self._disease_status = np.zeros(self._num_inds)
        self._disease_status[self._y > liability_cutoff] = 1


class PhenotypesBMI(Phenotypes):
    _sample_IDs: np.ndarray
    _pheno_df: pd.DataFrame

    def __init__(self):
        super().__init__()

    def initialize_from_file(self, filename, out, inds, logfile, num_phenotypes=1, phenotype_number_of_interest=1):
        if filename is None:
            raise ValueError("Provide file with BMI phenotype information using 'filename'")
        logfile.info("- Reading BMI phenotype information from " + filename)
        pheno_df = pd.read_csv(filename, names=["ID", "sex", "age", "BMI"])
        missing_in_phenotypes, added_in_phenotypes = find_missing_individuals(inds_tree=inds.names,
                                                                              inds_phenotype=pheno_df['ID'])
        logfile.info("- There are " + str(len(missing_in_phenotypes)) + " individuals missing from the phenotypes file "
                                                                        "and " + str(
            len(added_in_phenotypes)) + " individuals added. Will add missing ones with NA and "
                                        "remove added ones.")

        for i in missing_in_phenotypes:
            pheno_df.loc[len(pheno_df.index)] = [i, np.nan, np.nan, np.nan]

        for i in added_in_phenotypes:
            indexInd = pheno_df[(pheno_df['ID'] == i)].index
            pheno_df.drop(indexInd, inplace=True)

        pheno_df = sortPhenotypes(names_correct_order=inds.names, pheno_df=pheno_df)
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
            self._pheno_df.loc[
                (self._pheno_df['sex'] == sex) & (self._pheno_df['BMI'].notnull()), 'rank_inv_transform'] = norm.ppf(
                quants)

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

        self.standardize(out=out, inds=inds, logfile=logfile)

        super().write_to_file_gcta_eGRM(inds, out, logfile)


class PhenotypesSimulated(Phenotypes):
    _genetic_variance: float

    def __init__(self, variants, num_inds):
        super().__init__()

        self._random_noise = np.zeros(num_inds)
        self.betas = [0] * variants.number
        self._y = np.zeros(
            num_inds)  # must be zero and not np.empty because we add noise afterwards! np-empty leads to nan
        self._num_inds = num_inds

    @property
    def genetic_variance(self):
        return self._genetic_variance

    @genetic_variance.setter
    def genetic_variance(self, genetic_variance: float):
        self._genetic_variance = genetic_variance.find_missing_individuals

    def simulate(self, args: tparams, r: rg, logfile, variants_orig: tvar, inds: tind, trees: tskit, trees_orig: tskit, plots_dir: str, samp_ids: list):
        """
        Simulate phenotypes. If this is run on simulated trees for which the variants have not yet been downsampled
        (or 'typed'), trees and trees_orig can both be the original simulated trees. variants_orig is then the variants
        file produced when simulating. The downsampling will then remove some of the causal variants, so
        the causal variants will automatically be a mixture of typed and untyped variants.

        If you want to have a handle on how many causal variants are typed or untyped, you rely on the random
        downsampling process to decide which variants are typed. You must then define trees_orig to be the original
        simulated trees, trees to be the trees estimated based on the downsampled typed variants, and variants_orig
        to be the variants file that is produced by task 'downsampleVariantsWriteShapeit'

        :param samp_ids: output of tskit.sample() function. Necessary for tskit to know which variants are in the variants iterator
        :param trees_orig: the simulated ARG
        :param args: TArgs
        :param r: TRandomGenerator
        :param logfile: IndentedLoggerAdapter
        :param variants_orig: the original variants file produced during ARG simulation, or the result of task 'downsampleVariantsWriteShapeit'. it provides 'typed' status information of variants.
        :param inds: TInds
        :param trees: the simulated or estimated ARG
        :param plots_dir: str
        :return:
        """
        # simulate trait architecture
        self.simulate_trait_architecture(args=args, r=r, logfile=logfile, variants_orig=variants_orig,
                                         trees_orig=trees_orig, inds=inds, trees=trees, plots_dir=plots_dir,
                                         samp_ids=samp_ids)

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
            self._random_noise = self.simulate_env_noise_sd(sd_random_noise=args.pty_sd_envNoise,
                                                            pty_mean_envNoise=args.pty_mean_envNoise,
                                                            inds=inds,
                                                            random=r)
            logfile.info("- Simulated random noise with mean " + str(args.pty_mean_envNoise) + " and sd " +
                         str(args.pty_sd_envNoise) + ". The variance of the random noise is thus " +
                         str(np.var(self._random_noise)))
        else:
            raise ValueError("Must provide random noise distribution parameter. Either set noise sd with "
                             "'pty_sd_envNoise' or heritability with 'pty_h_squared'")

        self._y += self._random_noise

        if args.add_1_to_half_of_inds:
            half = int(inds.num_inds / 2)
            logfile.info("Adding 1 to phenotypes of first " + str(half) + " individuals")
            self._y[0:half] += 1
        # self.standardize(logfile)

        # write phenotypes to file
        self.write_sim_params_to_file(variants_orig, inds, args.out, logfile)

        self.filled = True

    def simulate_trait_architecture(self, args: tparams, r: rg, logfile: IndentedLoggerAdapter, variants_orig: tvar,
                                    inds: tind, trees: tskit, trees_orig: tskit, plots_dir: str, samp_ids: list):
        """
        Simulate phenotype's genetic architecture
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
            self.simulate_fixed(trees=trees_orig, inds=inds, causal_variant_indeces=args.pty_fixed_variant_indeces,
                                betas=args.pty_fixed_betas, logfile=logfile)

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

            # fig, ax = plt.subplots(1, figsize=(30, 30))
            var_index, pos = variants_orig.find_variant(typed=True, freq=args.single_variant_af,
                                                        interval=args.single_variant_interval,
                                                        random=r, logfile=logfile)
            # fig.tight_layout()
            # fig.set_size_inches(30, 30)
            # fig.savefig(plots_dir + 'allele_freq_spectrum.png', bbox_inches='tight')

            logfile.info("- Simulating a phenotypes based on the following typed variant index: " + str(
                var_index) + " at position " + str(
                variants_orig.info['position'][var_index]) + " with allele freq " + str(
                variants_orig.info['allele_freq'][var_index]) + " and the following betas: " + str(
                args.pty_fixed_betas))
            self.simulate_fixed(trees=trees_orig, inds=inds, causal_variant_indeces=[var_index],
                                betas=args.pty_fixed_betas, logfile=logfile)
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
            self.simulate_fixed(trees=trees_orig, inds=inds, causal_variant_indeces=[var_index],
                                betas=args.pty_fixed_betas, logfile=logfile)

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
            self.simulate_causal_region(variants=variants_orig, inds=inds, left_bound=causal_tree.interval.left,
                                        right_bound=causal_tree.interval.right,
                                        trees=trees,
                                        causal_mutations_effect_size_def=args.pty_sd_beta_causal_mutations,
                                        random=r,
                                        local_heritability=args.pty_h_squared,
                                        min_allele_freq_causal=args.min_allele_freq_causal,
                                        prop_causal_mutations=args.pty_prop_causal_mutations,
                                        max_allele_freq_causal=args.max_allele_freq_causal,
                                        allow_typed_causal_variants=args.allow_typed_causal_variants,
                                        logfile=logfile,
                                        samp_ids=samp_ids)

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
                                        trees=trees,
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
                                        logfile=logfile,
                                        samp_ids=samp_ids)

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
            self.simulate_fixed(trees=trees_orig, inds=inds, causal_variant_indeces=variant_indeces,
                                betas=fixed_betas, logfile=logfile)
            logfile.sub()

    @staticmethod
    def return_random_state(random):
        print(random.random.get_state()[1][0])
        print(random.random.uniform(0, 1, 1))

    def simulate_env_noise_sd(self, pty_mean_envNoise, sd_random_noise, random, inds):
        """
        Simulate random noise according to N(0, sd_random_noise)
        :param sd_random_noise: float
        :param pty_mean_envNoise:
        :param random: TRandomGenerator
        :return noise: np.array
        """
        noise = random.random.normal(loc=pty_mean_envNoise, scale=sd_random_noise, size=inds.num_inds)

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

    def simulate_fixed(self, trees: tskit, inds, causal_variant_indeces, betas, logfile):
        """
        Simulate phenotypes based on predefined causal variant positions and effects.
        This function is not efficient if there are many causal variants because the tskit variants iterator is not used

        :param variants: TVariants
        :param inds: TInds
        :param causal_variant_indeces: list
        :param betas: list
        :param logfile:
        :return:
        """

        causal_variants = []
        causal_pos = []
        for var in trees.variants():
            if var.site.id in causal_variant_indeces:
                causal_variants.append(var)
                causal_pos.append(var.site.position)

        if len(causal_variants) != len(betas):
            raise ValueError("must provide equal number of causal variants and betas to simulate fixed phenotype")

        for v, var in enumerate(causal_variants):
            # define beta
            self.betas[causal_variant_indeces[v]] = betas[v]
            allele_freq = sum(var.genotypes) / len(var.genotypes)
            allele_freq = min(allele_freq, 1 - allele_freq)

            self.add_one_variant_effect(tskit_variant_object=var,
                                        variant_index=var.site.id,
                                        inds=inds,
                                        beta=betas[v],
                                        allele_freq=allele_freq)

            logfile.info("- Simulated causal variant at position " + str(causal_pos[v]) + " at index " + str(
                causal_variant_indeces[v]) + " with beta " + str(round(betas[v], 3)) + " and allele freq " + str(
                allele_freq) + " resulting in a power of " + str(
                round(betas[v] ** 2 * allele_freq * (1 - allele_freq), 3)))

    def simulate_uniform(self, variants: tvar, inds: tind, prop_causal_mutations: float,
                         sd_beta_causal_mutations: float,
                         random: rg, logfile, mean_beta_causal_mutation=0.0):
        """
        DEPRECATED! Simulate phenotypes based on uniformly distributed causal variant positions and normally distributed effect sizes

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

        logfile.warning("Deprecated!")

        for v, var in enumerate(variants.variants):

            # if variants.info["typed"] == True:

            r = random.random.uniform(0, 1, 1)

            if r < prop_causal_mutations:
                allele_freq = variants.info['allele_freq'][var.site.id]

                # define beta
                beta = np.random.normal(loc=mean_beta_causal_mutation, scale=sd_beta_causal_mutations, size=None)
                self.betas[v] = beta

                self.add_one_variant_effect(tskit_variant_object=var,
                                            variant_index=var.site.id,
                                            inds=inds,
                                            beta=beta,
                                            allele_freq=allele_freq)

        logfile.info("- Simulated phenotypes based on " + str(
            len(self.causal_variants)) + " causal variants out of a total of " + str(variants.number) + ".")
        self.filled = True

    def add_one_variant_effect(self, tskit_variant_object: tskit.Variant, variant_index: int,
                               inds: tind, beta: float, allele_freq: float):

        # simulate phenotype
        if inds.ploidy == 1:
            self._y[tskit_variant_object.genotypes == 1] += beta
        else:
            genotypes = inds.get_diploid_genotypes(tskit_variant_object.genotypes)
            self._y[genotypes == 1] += beta
            self._y[genotypes == 2] += 2.0 * beta

        # save causal position
        self.causal_variants.append(tskit_variant_object)
        self.causal_betas.append(beta)
        self.causal_power.append(beta ** 2 * allele_freq * (1 - allele_freq))
        self.causal_variant_indeces.append(variant_index)

    def simulate_beta(self, variants: tvar, causal_mutations_effect_size_def: str, random: rg, variant_index: int,
                      local_heritability: float, num_causal_vars: int, ):
        # get beta
        sd = None
        # get effect size sd
        try:
            sd = float(causal_mutations_effect_size_def)
        except ValueError:
            causal_mutations_effect_size_def = causal_mutations_effect_size_def
            if causal_mutations_effect_size_def != "freq_dependent":
                raise ValueError("'sd_beta_causal_mutations' is set to unknown option. Must be a float or "
                                 "'freq_dependent'")
            f = variants.info['allele_freq'][variant_index]
            sd = (local_heritability / num_causal_vars) / np.sqrt(2 * f * (1 - f))

        beta = self.get_beta_normal(random, sd)
        self.betas[variant_index] = beta

        return beta

    def simulate_causal_region(self, trees: tskit, variants: tvar, inds, left_bound: float, right_bound: float,
                               causal_mutations_effect_size_def: str, local_heritability, prop_causal_mutations: float,
                               random, min_allele_freq_causal: float,
                               max_allele_freq_causal, logfile, allow_typed_causal_variants: bool, samp_ids: list):
        """
        Simulate causal effect sizes for variants within a region
        :param samp_ids:
        :param trees:
        :param variants:
        :param inds:
        :param left_bound:
        :param right_bound:
        :param causal_mutations_effect_size_def: Std. dev. for betas of causal mutations . If it can be converted to a
                                        float, betas will sampled from N(0, pty_sd_beta_causal_mutations). If set to
                                        'freq_dependent', betas will be sampled from
                                        N(0, [2 * f * (1 - f)]^{-0.5} * h2g / p),
                                        where h2g is the heritability of the trait and p is the number of causal SNPs.
        :param local_heritability:
        :param prop_causal_mutations: [0,1]
        :param random:
        :param min_allele_freq_causal: [0,1] Only variants with freq higher than this can be causal
        :param max_allele_freq_causal: [0,1] Only variants with freq lower than this can be causal
        :param logfile:
        :param allow_typed_causal_variants: If true, only also typed variants can be causal
        :return:
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
            logfile.info("- WARNING: Number of causal variants is 0! Selecting one causal variant at random.")
            indeces = range(len(info_window.index))
            c = random.random.choice(indeces)
            causal[c] = 1
            num_causal_vars = 1

        # add phenotypic effect to mutations that are uniformly distributed
        v_i = 0
        for var in trees.variants(samples=samp_ids):
            if var.site.id in info_window['var_index']:
                if causal[v_i] == 1:
                    allele_freq = variants.info['allele_freq'][var.site.id]

                    beta = self.simulate_beta(causal_mutations_effect_size_def=causal_mutations_effect_size_def,
                                              random=random,
                                              local_heritability=local_heritability,
                                              variant_index=var.site.id,
                                              variants=variants,
                                              num_causal_vars=num_causal_vars)

                    self.add_one_variant_effect(tskit_variant_object=var,
                                                variant_index=var.site.id,
                                                inds=inds,
                                                beta=beta,
                                                allele_freq=allele_freq)
                v_i += 1

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
        table['power'] = 0.0
        table.loc[self.causal_variant_indeces, 'power'] = self.causal_power
        table['var_genotypic_from_betas'] = np.repeat(float(inds.ploidy), variants.number) * np.array(
            self.betas) * np.array(
            self.betas) * np.array(table['allele_freq']) * (np.repeat(1.0, variants.number) - np.array(
            table['allele_freq']))
        table['var_genotypic_empiric'] = np.repeat(self._genetic_variance, variants.number)
        table['var_random'] = np.repeat(np.var(self._random_noise), variants.number)
        table['var_phenotypic'] = np.var(self._y)
        table.to_csv(out + "_pheno_causal_vars.csv", index=False, header=True)
