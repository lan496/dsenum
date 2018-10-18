import numpy as np

from permutation import Permutation


class Labeling(object):

    def __init__(self, hnf, num_type):
        self.hnf = hnf
        self.num_type = num_type
        self.num_site = np.prod(self.hnf.diagonal())

        self.permutation = Permutation(self.hnf)
        self.list_prm = self.permutation.get_permutation_list()

    def generate_possible_labelings(self):
        pass

    def reduce_by_translation_symmetry(self):
        pass

    def reduce_by_rotation_symmetry(self):
        pass
