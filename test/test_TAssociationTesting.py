import unittest
import pytest
import sys

sys.path.append("/home/linkv/git/argwas")
import TAssociationTesting as at
import stdpopsim
import numpy as np


@pytest.fixture(scope="package")
def simulate():
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

    return trees, trees_interval


class TestAssociationTesting:
    @pytest.mark.parametrize(
        "window_size, window_ends_true",
        [(20000, [120000, 140000, 160000, 180000, 200000]),
         (30000, [130000, 160000, 190000, 200000]),
         ],
    )
    def test_get_window_ends(self, simulate, window_size, window_ends_true):
        trees, trees_interval = simulate

        # get windows when skipping first tree, sequence length is multiple of window_size
        window_ends = at.get_window_ends(window_size=window_size, trees_interval=trees_interval)
        # get windows when skipping first tree, sequence length is not multiple of window_size
        assert window_ends == window_ends_true

    @pytest.mark.parametrize(
        "window_start, window_end, tree_start, tree_end, proportion_true",
        [(2000, 3000, 2000, 2020, 1.0),  # overlaps start
         (2000, 3000, 1980, 2000, 0.0),  # before start
         (2000, 3000, 1980, 2001, 1 / 21),  # overlaps start
         (2000, 3000, 1980, 3020, 1.0),  # overlaps whole window
         (2000, 3000, 3000, 3020, 0.0),  # starts after window
         (2000, 3000, 2500, 3500, 0.5),  # overlaps part of window
         (2000, 3000, 2500, 3000, 1.0),  # goes up to end of window
         ],
    )
    def test_get_proportion_of_tree_within_window(self, window_start, window_end, tree_start, tree_end,
                                                  proportion_true):
        proportion = at.get_proportion_of_tree_within_window(window_start=window_start,
                                                             window_end=window_end,
                                                             tree_start=tree_start,
                                                             tree_end=tree_end)

        assert proportion == proportion_true
