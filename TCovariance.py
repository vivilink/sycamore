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
import pickle


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
        self._covariance_matrix_haploid = None
        self._covariance_matrix_diploid = None
        self._covariance_matrix = None  # this is either the diploid or haploid matrix, depending on the inds
        self._mu = None

    @property
    def covariance_matrix(self):
        # do not check for this because in GRM it may be none because there are no variants!
        # if self._covariance_matrix is None:
        #     raise ValueError("Covariance matrix has not been calculated yet!")
        return self._covariance_matrix

    @covariance_matrix.setter
    def covariance_matrix(self, covariance_matrix):
        self.covariance_matrix = covariance_matrix

    @property
    def mu(self):
        return self._mu

    @mu.setter
    def mu(self, mu):
        self._mu = mu

    @staticmethod
    def remove_inds_with_missing_phenotypes(covariance_matrix, inds):
        print("indeces_to_remove", inds.get_indeces_inds_no_phenotype())
        indeces_to_remove = inds.get_indeces_inds_no_phenotype()
        cleaned_covariance_matrix = (
            np.delete(np.delete(covariance_matrix, indeces_to_remove, 0), indeces_to_remove, 1))

        return cleaned_covariance_matrix

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
        # print("covariance", covariance_matrix)
        # np.savetxt(out + '_covariance_matrix.csv', covariance_matrix, delimiter=',')
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
        self._covariance_matrix_haploid = None
        self._covariance_matrix_diploid = None
        self._covariance_matrix = None
        self._mu = None


class TCovarianceeGRM(TCovariance):
    def __init__(self):
        super().__init__()
        self.covariance_type = "eGRM"

    def finalize(self, inds):
        if self._covariance_matrix_haploid is None:
            return False
        if self._mu is None:
            raise ValueError("mu is not defined but covariance matrix is, this should not happen")

        self.normalize()

        if inds.ploidy == 2:
            self.calculate_diploid()
            self._covariance_matrix = self._covariance_matrix_diploid
        else:
            self._covariance_matrix = self._covariance_matrix_haploid

        return True

    def write(self, out, inds, covariances_picklefile):
        """
        Write covariance matrix to file so it can be used to test for association. In the case of eGRM, association
        testing is done with GCTA.

        @param covariances_picklefile: pickle file with covariances calculated for all windows
        @param write_covariance_picklefiles:
        @param out: str
        @param inds: TInds
        """
        if self._covariance_matrix_haploid is None:
            return False
        if inds.ploidy == 2:
            self.write_for_gcta(covariance_matrix=self._covariance_matrix_diploid, mu=self._mu, inds=inds, out=out)
            if covariances_picklefile is not None:
                pickle.dump(self._covariance_matrix_diploid, covariances_picklefile)
        else:
            self.write_for_gcta(covariance_matrix=self._covariance_matrix_haploid, mu=self._mu, inds=inds, out=out)
            if covariances_picklefile is not None:
                pickle.dump(self._covariance_matrix_haploid, covariances_picklefile)

        return True

    def add_tree(self, tree_obj, proportion, inds):
        """
        Calculate the unnormalized and haploid eGRM (=cov) so that it can be added to window eGRM, and also its
        number of expected mutations.

        @param tree_obj: TTree
        @param proportion: proportion of tree that is within window (necessary to calculate l(e), see Fan et al.)
        @param inds: TInds
        @return: np.array, float
        """

        cov, mu = tree_obj.get_unnormalized_eGRM(tree_obj=tree_obj, inds=inds)
        if cov is None:
            return None
        if self._covariance_matrix_haploid is None:
            self._covariance_matrix_haploid = proportion * cov
            self._mu = proportion * mu
        else:
            self._covariance_matrix_haploid += proportion * cov
            self._mu += proportion * mu

    def calculate_diploid(self):
        """
        As in Fan et al. 2022
        @return:
        """
        N = self._covariance_matrix_haploid.shape[0]
        maternals = np.array(range(0, N, 2))
        paternals = np.array(range(1, N, 2))
        self._covariance_matrix_diploid = 0.5 * (self._covariance_matrix_haploid[maternals, :][:, maternals]
                                                 + self._covariance_matrix_haploid[maternals, :][:, paternals]
                                                 + self._covariance_matrix_haploid[paternals, :][:, maternals]
                                                 + self._covariance_matrix_haploid[paternals, :][:, paternals])

    def normalize(self):
        """
        As in Fan et al. 2022, divide haploid matrix by total number of expected mutations in window and center by
        column and row.
        """
        self._covariance_matrix_haploid /= self._mu
        self._covariance_matrix_haploid -= self._covariance_matrix_haploid.mean(axis=0)
        self._covariance_matrix_haploid -= self._covariance_matrix_haploid.mean(axis=1, keepdims=True)


class TCovarianceGRM(TCovariance):
    def __init__(self):
        super().__init__()
        self.covariance_type = "GRM"

    def write(self, out, inds, covariances_picklefile):
        # checks
        if self._covariance_matrix is None:
            return False
        if np.trace(self._covariance_matrix) < 0:
            raise ValueError("Trace of matrix cannot be negative")
       # if inds.ploidy == 1 and not math.isclose(np.trace(self._covariance_matrix), inds.num_inds):
            # trace for haploids is expected to be equal to number of individuals (not true for diploids if they
            # are not in perfect HWE)
            # logfile.info("Trace of matrix is not equal to the number of individuals. Was expecting " + str(
            #    inds.num_inds) + " but obtained " + str(np.trace(self._covariance_matrix)))

        # write matrix
        self.write_for_gcta(covariance_matrix=self._covariance_matrix, mu=self._mu, inds=inds, out=out)

        # write matrix to picklefile if debugging
        if covariances_picklefile is not None:
            pickle.dump(self._covariance_matrix, covariances_picklefile)

        return True

    def calculate_GRM(self, window_beginning, window_end, variants, inds):
        """
        Calculate GRM matrix based on variants in a window (can bee tree interval)

        Parameters
        ----------
        window_beginning : int
        window_end : int
        variants : TVariantsFiltered
        inds : TInds

        """
        window_variants = np.array(variants.variants)[
            (variants.info['position'] >= window_beginning) &
            (variants.info['position'] < window_end) &
            (variants.info['typed'] == True) &
            (variants.info['allele_freq'] > 0.0)]
        num_vars = window_variants.shape[0]
        if num_vars == 0:
            return None, None

        # loop over variants
        M_sum = np.zeros(shape=(inds.num_inds, inds.num_inds))
        for v_i in range(num_vars):  # range(num_vars)
            gt_haploid = window_variants[v_i].genotypes
            if inds.ploidy == 1:
                gt = gt_haploid
            else:
                gt = inds.get_diploid_genotypes(gt_haploid)
            # need unmirrored allele freq
            af = np.sum(gt) / (inds.num_haplotypes)
            first = np.array([gt - inds.ploidy * af]).T
            second = np.array([gt - inds.ploidy * af])
            M = np.dot(first, second)
            M = M / (inds.ploidy * af * (1 - af))
            M_sum += M
        M_window = M_sum / float(num_vars)
        # M_window = self.remove_inds_with_missing_phenotypes(covariance_matrix=M_window, inds=inds)
        self._covariance_matrix = M_window
        self._mu = num_vars

        return self._covariance_matrix, self._mu

    def finalize(self):
        pass


class TCovarianceScaled(TCovariance):
    def __init__(self):
        super().__init__()

        self.covariance_type = "scaled"
        self._covariance_matrix_haploid = None
        self._covariance_matrix_diploid = None

    def clear(self):
        self._covariance_matrix_haploid = None
        self._covariance_matrix_diploid = None

    def write(self, covariance_matrix, mu, inds, out):
        if self._covariance_matrix_haploid is None:
            return False

        with open(out + '_GRM_covariance.txt', 'w') as f:
            np.savetxt(f, self._covariance_matrix_haploid)
        f.close()

        subprocess.call(["Rscript", os.path.dirname(sys.argv[0]) + "/create_gcta_GRM.R", out])

    def add_tree(self, tree_obj, proportion, inds):
        cov = tree_obj.get_covariance_scaled(inds=inds)
        self._covariance_matrix_haploid += proportion * cov

    def finalize(self):
        pass
