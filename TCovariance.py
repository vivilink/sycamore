#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 4 10:19:55 2022

@author: linkv
"""

import numpy as np
import subprocess
import sys
import os
import struct
import math


def get_covariance_object(covariance_name):
    covariance_obj = None
    if covariance_name == "eGRM":
        covariance_obj = TCovarianceeGRM()
    elif covariance_name == "GRM":
        covariance_obj = TCovarianceGRM()
    elif covariance_name == "scaled":
        covariance_obj = TCovarianceScaled()
    else:
        raise ValueError("Did not recognize " + str(covariance_name) + " as a covariance type")

    return covariance_obj


class TCovariance:
    def __init__(self):
        self.covariance_type = "base"
        self.covariance_matrix = None

    def write_for_gcta(self, covariance_matrix, mu, inds, out):
        """
        Write out eGRM in GCTA binary format.
        :param: covariance numpy.ndarray of expected relatedness
        :param: mu floating point number of expected mutations
        :param: numpy.ndarray/list of individual IDs
        :param: str of output
        :returns: None
        """

        # K = prefix_path.grm.bin; relatedness diagonal + lower diagonal
        # mu = prefix_path.grm.N.bin; number of shared mutations between individuals on diagonal + lower diagonal
        # samples = prefix_path.grm.id; 2 column text = family_id individual_id
        n, n = covariance_matrix.shape
        with open("{}.grm.bin".format(out), "wb") as grmfile:
            for idx in range(n):
                for jdx in range(idx + 1):
                    val = struct.pack("f", covariance_matrix[idx, jdx])
                    grmfile.write(val)

        with open("{}.grm.N.bin".format(out), "wb") as grmfile:
            for idx in range(n):
                for jdx in range(idx + 1):
                    val = struct.pack("f", mu)
                    grmfile.write(val)

        with open("{}.grm.id".format(out), "w") as grmfile:
            for idx in range(n):
                fid = 0
                iid = inds.names[idx]
                grmfile.write("\t".join([str(fid), str(iid)]) + os.linesep)

    def clear(self):
        self.covariance_matrix = None


class TCovarianceeGRM(TCovariance):
    def __init__(self):
        super().__init__()

        self.covariance_type = "eGRM"
        self.covariance_matrix = None
        self.mu = None

    def clear(self):
        self.covariance_matrix = None
        self.mu = None

    def write(self, out, inds):
        if self.covariance_matrix is None:
            return None
        if self.mu is None:
            raise ValueError("mu is not defined but covariance matrix is, this should not happen")
        self.write_for_gcta(covariance_matrix=self.covariance_matrix, mu=self.mu, inds=inds, out=out)

    def add_tree(self, tree_obj, proportion, inds):
        cov, mu = tree_obj.get_eGRM(tree_obj=tree_obj, inds=inds)
        if cov is None:
            return None
        if self.covariance_matrix is None:
            self.covariance_matrix = proportion * cov
        else:
            self.covariance_matrix += proportion * cov


class TCovarianceGRM(TCovariance):
    def __init__(self):
        super().__init__()

        self.covariance_type = "GRM"
        self.covariance_matrix = None

    def write(self, out, inds, mu, logfile):
        if self.covariance_matrix is None:
            return None
        if np.trace(self.covariance_matrix) < 0:
            raise ValueError("Trace of matrix cannot be negative")
        if inds.ploidy == 1 and not math.isclose(np.trace(self.covariance_matrix), inds.num_inds):
            # trace for haploids is expected to be equal to number of individuals (not true for diploids if they
            # are not in perfect HWE)
            logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
                inds.num_inds) + " but obtained " + str(np.trace(self.covariance_matrix)))
        self.write_for_gcta(covariance_matrix=self.covariance_matrix, mu=mu, inds=inds, out=out)

    def get_GRM(self, window_beginning, window_end, variants, inds):
        """
        Calculate GRM matrix based on variants

        Parameters
        ----------
        window_beginning : int
        window_end : int
        variants : TVariantsFiltered
        inds : TInds

        Returns
        -------
        np.array if there are variants that are typed and have allele freq > 0. Otherwise None.

        """
        # tree_variant_info = variants.info[(variants.info['tree_index'] == self.index) & (variants.info['typed'] == True) & (variants.info['allele_freq'] > 0.0)]
        tree_variants = np.array(variants.variants)[
            (variants.info['position'] >= window_beginning) &
            (variants.info['position'] < window_end) &
            (variants.info['typed'] == True) &
            (variants.info['allele_freq'] > 0.0)]
        num_vars = tree_variants.shape[0]
        if num_vars == 0:
            return None, None

        # loop over variants
        M_sum = np.zeros(shape=(inds.num_inds, inds.num_inds))
        for v_i in range(num_vars):  # range(num_vars)
            gt_haploid = tree_variants[v_i].genotypes
            if inds.ploidy == 1:
                gt = gt_haploid
            else:
                gt = inds.get_diploid_genotypes(gt_haploid)
            # need unmirrored allele freq
            af = np.sum(gt) / (2 * inds.num_inds)
            first = np.array([gt - inds.ploidy * af]).T
            second = np.array([gt - inds.ploidy * af])
            M = np.dot(first, second)
            M = M / (inds.ploidy * af * (1 - af))
            M_sum += M
        M_window = M_sum / float(num_vars)
        return M_window, num_vars


class TCovarianceScaled(TCovariance):
    def __init__(self):
        super().__init__()

        self.covariance_type = "scaled"
        self.covariance_matrix = None

    def write_for_gcta(self, covariance_matrix, mu, inds, out):
        if self.covariance_matrix is None:
            return None

        with open(out + '_GRM_covariance.txt', 'w') as f:
            np.savetxt(f, self.covariance_matrix)
        f.close()

        subprocess.call(["Rscript", os.path.dirname(sys.argv[0]) + "/create_gcta_GRM.R", out])

    def add_tree(self, tree_obj, proportion, inds):
        cov = tree_obj.get_covariance_scaled(inds=inds)
        self.covariance_matrix += proportion * cov

    # -----------------------------

    def calculate_and_write_covariance_matrix_to_gcta_file(self, ts_object, variants, tree_obj, inds, covariance_type,
                                                           out, skip_first_tree, logfile):
        """
        Writes covariance and other files necessary to run gcta. The program egrm does that automatically, the scaled

        Parameters
        ----------
        ts_object : tskit.treeSequence
        variants : TVariants
        tree_obj : TTree.py
        inds : TInds
        covariance_type : str
        out : str
        skip_first_tree: bool
        logfile : IndentedLoggerAdapter

        Raises
        ------
        ValueError
            If the covariance_type is not a recognized method.

        Returns
        -------
        Covariance: ndarray(inds.num_inds, inds.num_inds).

        """
        if covariance_type == "scaled":
            covariance = tree_obj.get_covariance_scaled(inds=inds)
            write_covariance_matrix_R(covariance=covariance, out=out)

        elif covariance_type == "eGRM":
            # trees = ts_object.keep_intervals(np.array([[tree_obj.start, tree_obj.end]]), simplify=True)
            covariance, mu = tree_obj.get_eGRM(tskit_obj=ts_object, tree_obj=tree_obj, inds=inds)
            if covariance is None:
                return None
            write_covariance_matrix_bin(covariance=covariance, mu=mu, inds=inds, out=out)

            # if np.trace(covariance) != inds.num_inds:
            # raise ValueError("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
            # inds.num_inds) + " but obtained " + str(np.trace(covariance)))
            # logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
            #     inds.num_inds) + " but obtained " + str(np.trace(covariance)))

        elif covariance_type == "GRM":
            covariance, mu = tree_obj.get_GRM(variants=variants, inds=inds)
            if covariance is None:
                return None
            if np.trace(covariance) < 0:
                raise ValueError("Trace of matrix cannot be negative")
            if inds.ploidy == 1 and not math.isclose(np.trace(covariance), inds.num_inds):
                # trace for haploids is expected to be equal to number of individuals (not true for diploids if they
                # are not in perfect HWE)
                logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
                    inds.num_inds) + " but obtained " + str(np.trace(covariance)))
            write_covariance_matrix_bin(covariance=covariance, mu=mu, inds=inds, out=out)

        else:
            raise ValueError("Did not recognize " + str(covariance_type) + " as a covariance type")

        return covariance
