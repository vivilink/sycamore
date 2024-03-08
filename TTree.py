#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 17:43:29 2021

@author: linkv
"""
import numpy as np
import pandas as pd
import tskit
import TCovariance
from egrm import varGRM_C
from egrm import egrm_one_tree_no_normalization_C
from arg_needle_lib import deserialize_arg
from arg_needle_lib import arg_to_tskit
import arg_needle_lib

# Silence other loggers
import logging
logging.getLogger(arg_needle_lib.__name__).setLevel(logging.WARNING)

logging.basicConfig(level=logging.DEBUG)

logger_blocklist = [
    "arg_needle_lib",
]
for module in logger_blocklist:
    logging.getLogger(module).setLevel(logging.WARNING)
# ---------------------
# read tree file
# ---------------------


class TTrees:
    def __init__(self, tree_file, skip_first_tree, trees_interval, trees_interval_start, trees_interval_end, logfile):
        self.trees, self.actual_trees_interval = self.read_trees(tree_file,
                                                                 trees_interval,
                                                                 trees_interval_start,
                                                                 trees_interval_end,
                                                                 logfile)

        self.number = self.trees.num_trees  # this contains potential empty border trees
        self.actual_trees_interval = self.remove_border_trees_if_empty(skip_first_tree)

    def read_trees(self, tree_file, trees_interval, trees_interval_start, trees_interval_end, logfile):
        """
        Read in tskit file and extract only a certain region, if defined by user

        @param tree_file: tskit file
        @param trees_interval:
        @param trees_interval_start:
        @param trees_interval_end:
        @param logfile:
        @return:
        """
        logfile.info("- Reading tree from " + tree_file + ". Assuming it is in tskit format unless ending with '.argn'")

        if tree_file.endswith(".argn"):
            logfile.add()
            logfile.info("- Converting arg-needle format to tskit")
            arg = deserialize_arg(tree_file)
            trees_full = arg_to_tskit(arg, mutations=True, sample_permutation=None)
        else:
            trees_full = tskit.load(tree_file)

        # define interval
        if trees_interval is None:
            trees_interval = [0, trees_full.sequence_length]
        if trees_interval_start:
            trees_interval[0] = trees_interval_start
        if trees_interval_end:
            trees_interval[1] = trees_interval_end

        if trees_interval != [0, trees_full.sequence_length]:
            logfile.info(
                "- Extracting trees overlapping the following interval: " + str(trees_interval))
            trees_extract = trees_full.keep_intervals([trees_interval], simplify=True)
        else:
            trees_extract = trees_full.simplify()

        trees_extract = TTrees.remove_monomorphic(trees_extract)

        return trees_extract, trees_interval

    def remove_border_trees_if_empty(self, skip_first_tree):
        """
        @param skip_first_tree:
        remove untestable trees (can happen for first and last tree if tskit object was extracted from a larger one)
        """
        actual_trees_interval = self.actual_trees_interval
        first = self.trees.first()
        first_tree_obj = TTree(first)
        if not first_tree_obj.is_testable(skip_first_tree=skip_first_tree):
            # tree sequence starts at second tree
            actual_trees_interval[0] = self.trees.breakpoints(as_array=True)[1]

        last_tree_obj = TTree(self.trees.last())
        if not last_tree_obj.is_testable(skip_first_tree=skip_first_tree):
            # tree sequence ends at second to last tree
            actual_trees_interval[1] = self.trees.breakpoints(as_array=True)[-2]

        return actual_trees_interval

    @staticmethod
    def writeStats(ts_object, out, logfile):
        """
        Parameters
        ----------
        ts_object : tskit.treeSequence
            A tskit tree sequence of interest.
        out : str
            Output prefix.
        logfile : IndentedLoggerAdapter

        Returns
        -------
        None.

        """
        info = pd.DataFrame(index=range(ts_object.num_trees),
                            columns=['index', 'start', 'end', 'length', 'num_mutations', "total_branch_length",
                                     "num_roots"])
        for tree in ts_object.trees():
            info.iloc[tree.index] = [tree.index, tree.interval.left, tree.interval.right,
                                     tree.interval.right - tree.interval.left, len(list(tree.mutations())),
                                     tree.total_branch_length, len(tree.roots)]

        logfile.info("- Writing tree info to file '" + out + "_trees_statistics.csv'")
        info.to_csv(out + "_trees_statistics.csv", header=True, index=False)

    @staticmethod
    def extract_single_tree(ts_object, out, logfile, position):
        focal_tree = ts_object.at(position)
        trees = ts_object.keep_intervals([np.array(focal_tree.interval)], simplify=True)
        trees.dump(out + "_focal.trees")
        logfile.info("- Wrote trees with " + str(focal_tree.interval) + " to " + out + "_focal.trees")

    @staticmethod
    def remove_monomorphic(trees):
        """
        From egrm/manuscript/simulate, needed to remove variants that are not polymorphic in one of the populations
        if original tree simulation contained multiple populations
        @param trees: tskit.treeSequence
        @return: tskit.treeSequence
        """
        tables = trees.tables
        tables.sites.clear()
        tables.mutations.clear()
        n = trees.num_samples
        for tree in trees.trees():
            for site in tree.sites():
                visited = False
                for mutation in site.mutations:
                    k = tree.num_samples(mutation.node)
                    if 0 < k < n:
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


class TTree:
    def __init__(self, tree_iterator):
        self.tree: tskit.TreeIterator = tree_iterator
        self.start: float = tree_iterator.interval.left
        self.end: float = tree_iterator.interval.right
        self.length: float = tree_iterator.length
        self.index: int = tree_iterator.index
        self.height: float = -1.0
        # the first marginal tree has no root and no height in ARG simulated with stdpopsim
        if len(self.tree.roots) == 1:
            self.height = self.tree.time(self.tree.root)

        self.covariance: np.array = None
        self.covariance_scaled: TCovariance.TCovarianceScaled = TCovariance.TCovarianceScaled()
        # self.covariance_diploid: np.array = None
        # self.covariance_scaled_diploid: np.array = None
        self.eGRM: TCovariance.TCovarianceeGRM = TCovariance.TCovarianceeGRM()
        self.EK_relate_mu = None
        # self.eGRM_diploid: np.array = None

    def test_PSD(self, A, tol=1e-8):
        E = np.linalg.eigvalsh(A)
        if not np.all(E > -tol):
            raise ValueError("Covariance matrix has negative eigenvalues")
        return np.all(E > -tol)

    def TMRCA(self, num_haplotypes):
        """
        see tskit to understand these classes: https://tskit.dev/tskit/docs/stable/python-api.html#the-tree-class
        
        - the goal of this function is to get pairwise branch lengths between all nodes of the tree
        - we go through all internal nodes and update the cells of all their descendants at once
        - the tree.time of a node gets larger backwards in time
        - the idea is to fill the TMRCA matrix with the total height of the tree and then to 
        deduct for all descendants of one node the branches from farther back in time
        for example, the terminal nodes will have all branch lengths that connect their ancestors deducted 
        -> they are only separated by the terminal branches
        Since the height of the tree can only be known by finding the oldest node, which we don't know a priori 
        we deduct the branch lengths from zero first, using the loop over nodes to find the height, and then add the height
    
        Parameters
        ----------
        num_haplotypes : int
            number of haplotypes.
    
        Returns
        -------
        tmrca : pairwise distance between two haplotypes. square matrix of dimension num_haplotypes and type float.
    
        """
        tmrca = np.zeros([num_haplotypes, num_haplotypes])
        # self.height = 0
        for c in self.tree.nodes():
            descendants = list(self.tree.samples(c))
            n = len(descendants)
            if (n == 0 or n == num_haplotypes or self.tree.time(
                    c) == 0):  # The branch length for a node that has no parent (e.g., a root) is defined as zero.
                continue
            t = self.tree.time(self.tree.parent(c)) - self.tree.time(c)
            tmrca[np.ix_(descendants, descendants)] -= t
            # self.height = max(self.height, self.tree.time(self.tree.parent(c))) #time returns the time of a node
        tmrca += self.height
        np.fill_diagonal(tmrca, 0)

        # tmrca = (tmrca + tmrca.T) / 2

        return tmrca

    def get_covariance(self, inds):
        """
        Calculate variance-covariance between haplotypes. Total height of tree = distance_ij + covariance_ij.
        Variance of one haplotype = total height of tree.

        Returns
        -------
        Variance-covariance matrix.

        """

        if self.height == -1:
            raise ValueError("Cannot calculate covariance from tree with multiple roots")

        if self.covariance is None:
            TMRCA = self.TMRCA(inds.num_haplotypes)
            self.covariance = -TMRCA + self.height
        return self.covariance

    def get_covariance_scaled(self, inds) -> np.array:
        """
        Caclulate scaled variance-covariance between haplotypes. This allows gcta REML to run without numeric issues
        such as singular Information matrix.

        Returns
        -------
        Scaled variance-covariance matrix

        """
        # calculate haploid covariance and covariance scaled
        if self.covariance is None:
            TMRCA = self.TMRCA(inds.num_haplotypes)
            self.covariance = -TMRCA + self.height

        if self.covariance_scaled._covariance_matrix_haploid is None:
            self.covariance_scaled._covariance_matrix_haploid = self.covariance * float(inds.num_haplotypes) / np.trace(
                self.covariance)

        if inds.ploidy == 1:
            return self.covariance_scaled._covariance_matrix_haploid

        # calculate diploid covariance scaled
        else:
            if self.covariance_scaled._covariance_matrix_diploid is None:
                # logfile.add()

                # add together unscaled covariance of haplotypes of one individual
                self.covariance_scaled._covariance_matrix_diploid = np.zeros([inds.num_inds, inds.num_inds])

                # off-diagonals upper triangle (this only works if ind assignment is equal to neighboring pairs!)
                for i in range(inds.num_inds):
                    # if i % 100 == 0:
                    #     logfile.info("- Filling diploid covariance matrix for individual " + str(i) + " of " + str(inds.num_inds))
                    i1 = i * 2
                    i2 = i1 + 1
                    for j in range(i + 1, inds.num_inds):
                        j1 = j * 2
                        j2 = j1 + 1
                        self.covariance_scaled._covariance_matrix_diploid[i, j] = self.covariance[i1, j1] + \
                                                                                  self.covariance[i1, j2] + \
                                                                                  self.covariance[i2, j1] + \
                                                                                  self.covariance[i2, j2]

                # lower triangle
                self.covariance_diploid = self.covariance_scaled._covariance_matrix_diploid \
                                          + self.covariance_scaled._covariance_matrix_diploid.T

                # diagonals
                for ii in range(inds.num_inds):
                    ii1 = inds.get_haplotypes(ii)[0]
                    ii2 = inds.get_haplotypes(ii)[1]
                    self.covariance_diploid[ii, ii] = 2.0 * self.height + 2.0 * self.covariance[ii1, ii2]

                self.covariance_scaled._covariance_matrix_diploid = self.covariance_scaled._covariance_matrix_diploid \
                                                                    * float(inds.num_inds) \
                                                                    / np.trace(
                    self.covariance_scaled._covariance_matrix_diploid)

                # logfile.sub()

            return self.covariance_diploid

    def get_eGRM(self, tree_obj, inds):
        """       
        Parameters
        ----------
        tree_obj : TTree
        tskit_obj : tskit.treeSequence
            the tskit tree sequence that contains the single tree to be tested.
        inds : TInds
            Which haplotypes are assigned to the same individual.

        Returns
        -------
        local eGRM as calculated by egrm (Fan et al. 2022).
        """

        if self.eGRM._covariance_matrix_haploid is None:
            # extract tree and write to file
            # TTrees.extract_single_tree(tree_obj=tree_obj, out=out, logfile=logfile, position=self.start)

            EK_relate, EK_relate_mu = varGRM_C(tree=tree_obj.tree, num_samples=inds.num_haplotypes)
            self.eGRM = EK_relate
            self.EK_relate_mu = EK_relate_mu

        if inds.ploidy == 2:
            if self.eGRM._covariance_matrix_diploid is None:
                self.eGRM.calculate_diploid()
            return self.eGRM._covariance_matrix_diploid, self.EK_relate_mu
        else:
            return self.eGRM._covariance_matrix_haploid, self.EK_relate_mu

    @staticmethod
    def get_unnormalized_eGRM(tree_obj, inds):
        # cov, mu = egrm_one_tree_no_normalization(tree=tree_obj.tree, num_samples=inds.num_haplotypes)
        cov, mu = egrm_one_tree_no_normalization_C(tree=tree_obj.tree, num_samples=inds.num_haplotypes)
        return cov, mu

    def solving_function(self, array, inds):
        covariance = self.covariance(inds.num_haplotypes)
        # covariance = (X+X.T)/2
        # np.fill_diagonal(covariance, 0)
        # print(np.diagonal(covariance))
        # print(np.shape(covariance))
        # print(covariance)
        self.test_PSD(covariance)
        inv = np.linalg.inv(covariance)
        tmp = np.dot(inv, array)
        # print("shape of my dot product",np.shape(tmp))
        return (tmp)

    def is_testable(self, skip_first_tree):
        # TODO: this condition is because if you extract a region from tskit sequence, the first tree goes
        #  from zero to the first tree. This causes problems with eGRM. Needs to be investigated what the
        #  problem is and a better condition needs to be found!
        if self.height != -1 and not (skip_first_tree and self.index == 0):
            return True
        else:
            return False
