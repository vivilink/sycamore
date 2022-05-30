from numpy.random import RandomState


class TRandomGenerator:
    def __init__(self, seed):
        self.seed = seed
        self.random = RandomState(seed)
