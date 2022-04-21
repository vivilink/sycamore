#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 17:44:45 2021

@author: linkv
"""
import numpy as np
import pandas as pd


class TVariants:
    """
    A class containing all variants of a tree, all variants are typed
    """

    def __init__(self, ts_object, samp_ids, pos_int):
        self._variants = list(ts_object.variants(samples=samp_ids))
        self._number = len(self._variants)
        self._number_typed = self._number
        self._positions = np.empty(self._number)
        self._alleleFreq = np.empty(self._number)
        self._info_columns = ['var_index', 'position', 'allele_freq', 'num_alleles', 'typed', 'tree_index']
        self._info = pd.DataFrame(index=range(self._number), columns=self._info_columns)

        for v, var in enumerate(list(ts_object.variants(samples=samp_ids))):
            tmp = sum(var.genotypes) / len(var.genotypes)

            if len(np.bincount(var.genotypes)) > 2:
                print(var.genotypes)
                print(var.alleles)
                print("af is: ", tmp, "expected homo alt: ", tmp * tmp)
                print("Invalid genotypes encountered at variant index ", v,
                      ". Currently we cannot accept simulated individuals to be diploid. Genotype counts are ",
                      np.bincount(var.genotypes))
            af = min(tmp, 1 - tmp)

            # determine position
            pos = -1
            if pos_int:
                pos = round(var.site.position)
                if v > 0 and pos == self._info['position'][v - 1]:
                    # print("am in special case for v", v)
                    pos += 1
            else:
                pos = var.site.position

            self._positions[v] = pos
            self._alleleFreq[v] = af

        self._info['var_index'] = np.arange(self._number)
        self._info['position'] = self._positions
        self._info['allele_freq'] = self._alleleFreq
        self._info['typed'] = "True"

        # set indeces of trees the variants belong to, digitize starts at 1 but tree indeces at 0
        self._info['tree_index'] = np.digitize(self._info['position'], ts_object.breakpoints(as_array=True)) - 1

    @property
    def variants(self):
        return self._variants

    @variants.setter
    def variants(self, variants):
        self._variants = variants

    @property
    def number(self):
        return self._number

    @number.setter
    def number(self, number):
        self._number = number

    @property
    def number_typed(self):
        return self._number_typed

    @number_typed.setter
    def number_typed(self, number_typed):
        self._number_typed = number_typed

    @property
    def info(self):
        return self._info

    @info.setter
    def info(self, info):
        self._info = info

    def write_variant_info(self, ts_object, samp_ids, out, logfile):
        logfile.info("- Writing variant info to file '" + out + "_sample_variants.csv'")
        self._info.to_csv(out + "_sample_variants.csv", header=True, index=False)


class TVariantsFiltered(TVariants):
    """
    A class containing all variants of a tree, typed status depends on either what is given by the
    filtered_variants_file, or if none given it depends on the allele frequency filteres and prop_typed_variants filters
    """

    def __init__(self, ts_object, samp_ids, min_allele_freq, max_allele_freq, prop_typed_variants, pos_int, random,
                 logfile, filtered_variants_file=None):
        """
        Parameters
        ----------
        ts_object : TreeSequence
            This tree sequence must match the file provided with filtered_variants_file, if provided.
        samp_ids : numpy.ndarray (dtype=np.int32)
            Samples given by ts_object.samples()
        min_allele_freq : float
        max_allele_freq : float
        prop_typed_variants : float
        pos_int : bool
            Defines if positions of variants are to be saved as integers.
        random : TRandomGenerator
        logfile : logger
        filtered_variants_file : str, optional
            If provided, it determines the typed status of the variants. Must match ts_object. 
        """

        # TODO: When should the filtered_variants_file be provided? If it is provided, the allele frequencies will be
        #  the original allele frequencies and the typed status will be predefined
        # TODO: Also they may also already be filtered for min and max freq. However, they may not match the tree file?
        # TODO: I don't understand if I should use the variants as a list or not. I don't know how to save them if not
        #  as a list (self._variants = trees._variants() or trees._variants does not work. maybe as np.array? but then
        #  enumerate in TPhenotypes does not work
        # TODO: develop an iterator for variants that only goes over the typed ones

        # build variant object from tree file -> filter!
        if filtered_variants_file is None:

            super().__init__(ts_object=ts_object, samp_ids=samp_ids, pos_int=pos_int)

            self._number_typed = -1

            logfile.info("- Building variant information from scratch based on simulated trees")

            if min_allele_freq < 0 or max_allele_freq < 0 or min_allele_freq > 1 or max_allele_freq > 1 or min_allele_freq > max_allele_freq:
                raise ValueError("allele frequency filters are nonsensical")

            # determine typed status for each variant
            for v, var in enumerate(ts_object.variants(samples=samp_ids)):
                # is variant typed? Depends on freq filter and proportion typed
                af = self._info['allele_freq'].iloc[v]
                if af >= min_allele_freq and af <= max_allele_freq:
                    if prop_typed_variants == 1:
                        typed = True
                    else:
                        r = random.random.uniform(0, 1, 1)
                        if r < prop_typed_variants:
                            typed = True
                        else:
                            typed = False
                else:
                    typed = False

                # change in table
                self._info["typed"].iloc[v] = typed

        # variants are already filtered -> read from file!
        else:
            logfile.info("- Reading variant information from " + filtered_variants_file)
            self._info = pd.read_csv(filtered_variants_file)

            if prop_typed_variants is not None:
                logfile.info("WARNING: 'prop_typed_variants' has no effect when variants file is already provided")

            # rename if variant file is in old format
            if 'index' in self._info.columns:
                self._info.rename(columns={'index': 'var_index'}, inplace=True)

                # check if file is ok
            if self._info.columns != self._info_columns:
                raise ValueError("Columns of variants file do not match the current standard")
            if len(self._info['var_index']) != self._number != len(self._variants):
                raise ValueError("Variant file " + filtered_variants_file + " contains " + str(
                    len(self._info['var_index'])) + " variants, expected " + str(self._number))
            if 'tree_index' not in self._info.columns:
                # set indeces of trees the variants belong to, digitize starts at 1 but tree indeces at 0
                self._info['tree_index'] = np.digitize(self._info['position'], ts_object.breakpoints(as_array=True)) - 1

            # remove variants that now do not pass the frequency filter
            self._info['typed'].iloc[self._info['allele_freq'] < min_allele_freq] = False
            self._info['typed'].iloc[self._info['allele_freq'] > max_allele_freq] = False

        # set number typed for both cases (built from scratch or file)
        self._number_typed = self._info['typed'].value_counts()[True]

    def print_genotypes(self, index):
        file = "genotypes_variant" + str(index) + ".txt"
        self._info['variant'][index].genotypes.tofile(file=file)

    # def fill_diploidGenotypes(self, individuals):
    #     for v in self._info['variant']:
    #         v.site.metadata = []

    def write_variant_info(self, name, logfile):
        logfile.info("- Writing variant info to file '" + name + "_filtered_sample_variants.csv'")
        self._info.to_csv(name + "_filtered_sample_variants.csv", header=True, index=False)

    def write_shapeit2(self, name, inds, logfile):
        """
        Write files in SHAPEIT2 format, to be used as input by RELATE (https://myersgroup.github.io/relate/input_data.html)

        Parameters
        ----------
        name : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        haps = pd.DataFrame(index=range(self._number_typed), columns=range(5 + inds.num_haplotypes))
        info_typed = self._info.loc[self._info['typed'] == True]
        info_typed['index'] = range(self._number_typed)
        info_typed.set_index(info_typed['index'], drop=True, inplace=True)

        haps.iloc[:, 0] = np.repeat(1, self._number_typed)
        haps.iloc[:, 1] = '.'
        haps.iloc[0:self._number_typed, 2] = info_typed['position']
        haps.iloc[:, 3] = 'A'
        haps.iloc[:, 4] = 'T'

        logfile.info("- Building haplotypes for typed variants")

        index = 0
        for v, var in enumerate(self._variants):
            # print(v, self._info.iloc[v]['typed'])
            if self._info.iloc[v]['typed']:
                # if self._info.iloc[v]['position'] == None :
                #     print(self._info.iloc[v])
                haps.iloc[index,
                5:] = var.genotypes  # can't use v for index because it loops over all variants, not only typed ones
                # if index in [9,10, 11, 12]:
                #     print("v", v, "positions\n", self._info.iloc[v])
                index += 1

        logfile.info("- Writing haplotypes in Shapeit2 format to file '" + name + "_variants.haps'")
        haps.to_csv(name + "_variants.haps", sep=' ', header=False, index=False)

    def find_variant(self, typed, freq, interval, out, subplot, random, logfile):
        """
        Find a variant that fits criteria to simulate fixed genotype depending on only one variant

        Parameters
        ----------
        typed : bool
            Should the variant returned be typed or not.
        freq : float
            Requested allele frequency of the variant. If there is no variant with this exact frequency within the requested interval, 0.001 will be deducted from the freq repeatedly until a variant is found.
        interval : list of floats
            Requested genomic interval within which the variant should be.

        Raises
        ------
        ValueError
            If the search for a variant within the requested genomic interval ends up with the allele freq being negative, the search is stopped an an error is thrown.

        Returns
        -------
        Index of the variant found, can be used to simulate fixed phenotype.

        """

        # check if interval is valid
        num_lines = self._info[(self._info['position'] >= interval[0]) & (self._info['position'] <= interval[1])].shape[
            0]
        if num_lines == 0:
            raise ValueError("The interval " + str(interval) + " contains no variants")

        # make subset of variants in interval
        info_interval = self._info.loc[
            (self._info['position'] >= interval[0]) & (self._info['position'] <= interval[1])]
        # info_interval['index'] = range(self._info[(self._info['position'] >= interval[0]) & (self._info['position'] <= interval[1])].shape[0])
        # info_interval.set_index(info_interval['index'], drop=True, inplace=True)

        # check if there are variants with requested typed status
        num_lines = info_interval[info_interval['typed'] == typed].shape[0]
        if num_lines == 0:
            raise ValueError("There are no variants of typed status '" + str(typed) + "'")

        # make subset of variants in interval with correct typed status
        info_interval = info_interval.loc[info_interval['typed'] == typed]
        # info_interval['index'] = range(num_lines)
        # info_interval.set_index(info_interval['index'], drop=True, inplace=True)

        # plot allele freq spectrum
        n, bins, patches = subplot.hist(info_interval['allele_freq'], density=False)
        subplot.axvline(x=freq, color="red", lw=1)
        subplot.set(xlabel='Allele freq', ylabel='count',
                    title='Histogram of af in requested interval and typed status')

        # find first variant with requested allele frequency
        info = info_interval[np.isclose(np.array(info_interval['allele_freq'], dtype=float), freq, rtol=0.01)]

        if info.shape[0] < 1:
            logfile.info(
                "- Did not find locus with requested af " + str(freq) + ". Adapting af in increments of 0.001.")

            # set direction of search
            r = random.random.uniform(0, 1, 1)
            if r < 0.5:
                step = -0.001
            else:
                step = 0.001

            freq_orig = freq

            while info.shape[0] < 1 and 0 <= freq <= 0.5:
                # remove or add small value to freq until a locus is found
                freq = round(freq + step, 3)
                info = info_interval[np.isclose(np.array(info_interval['allele_freq'], dtype=float), freq, rtol=0.01)]

            # if loop was left because out of bounds, search in other direction
            if freq < 0 or freq > 0.5:
                logfile.warning(
                    "Allele frequency became negative or exceeded 0.5 while searching for locus with requested af " + str(
                        freq_orig) + " in interval " + str(interval) + ". Starting search in opposite direction.")
                step = -step

                # set search back to starting freq and go in other direction
                freq = freq_orig

                while info.shape[0] < 1 and 0 <= freq <= 0.5:
                    freq = round(freq + step, 3)
                    info = info_interval[
                        np.isclose(np.array(info_interval['allele_freq'], dtype=float), freq, rtol=0.01)]

            if freq < 0 or freq > 0.5:
                raise ValueError("Could not find locus with requested allele frequency")

        subplot.axvline(x=freq, color="black", lw=1)

        # logfile.info("--> Found variant with freq " + str(freq) + " within the following interval: " + str(interval))
        return info.iloc[0]['var_index'], info.iloc[0]['position']
