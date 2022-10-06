import unittest
import pytest
import sys

sys.path.append("/home/linkv/git/argwas")
import TAssociationTesting as at
import stdpopsim
import numpy as np


class TestAssociationTesting(unittest.TestCase):

    def test_get_window_ends(self):
        # simulate full trees
        species = stdpopsim.get_species("HomSap")
        contig = species.get_contig("chr22")
        model = species.get_demographic_model("OutOfAfrica_3G09")
        samples = model.get_samples(20, 0, 0)
        engine = stdpopsim.get_engine("msprime")
        trees_full = engine.simulate(model, contig, samples, discrete_genome=True)

        # extract an interval
        trees_interval = [100000, 200000]
        interval = np.array(trees_interval)
        trees = trees_full.keep_intervals([interval], simplify=True)

        # get windows when skipping first tree, sequence length is multiple of window_size
        window_size = 20000
        window_ends = at.get_window_ends(ts_object=trees, window_size=window_size, trees_interval=trees_interval)
        assert window_ends == [120000, 140000, 160000, 180000, 200000]

        # get windows when skipping first tree, sequence length is not multiple of window_size
        window_size = 30000
        window_ends = at.get_window_ends(ts_object=trees, window_size=window_size, trees_interval=trees_interval)
        assert window_ends == [130000, 160000, 190000, 200000]


if __name__ == '__main__':
    unittest.main()
