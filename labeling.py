import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import numpy as np

from permutation import Permutation
from smith_normal_form import smith_normal_form


class Labeling(object):

    def __init__(self, hnf, num_type, list_rotations=None):
        self.hnf = hnf
        self.num_type = num_type
        self.num_site = np.prod(self.hnf.diagonal())

        self.permutation = Permutation(self.hnf, list_rotations)
        # assuming that the 0-th element of permutaions is identity operation
        self.prm_t = self.permutation.get_translation_permutations()
        self.prm_all = self.permutation.get_symmetry_operation_permutaions()

    def act_permutation(self, labeling, prm):
        ret = tuple([labeling[prm[i]] for i in range(self.num_site)])
        return ret

    def generate_possible_labelings(self):
        list_labelings = []

        for labeling in itertools.product(range(self.num_type), repeat=self.num_site):
            if len(set(labeling)) == self.num_type:
                list_labelings.append(labeling)

        return list_labelings

    def remove_duplicates(self, labelings):
        unique_labelings = []
        flags = ["unvisited" for _ in labelings]
        for i, lbl in enumerate(labelings):
            if flags[i] != "unvisited":
                continue

            # remove superperiodic
            if any([(self.act_permutation(lbl, prm) == lbl) for prm in self.prm_t[1:]]):
                flags[i] = "superperiodic"
                continue

            # remove symmetry operation and exchanging duplicates
            for type_prm in itertools.permutations(range(self.num_type)):
                lbl_ex = [type_prm[e] for e in lbl]
                for prm in self.prm_all:
                    if (lbl_ex == lbl) and (prm == range(self.num_site)):
                        continue
                    operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                    idx = labelings.index(operated_lbl_ex)
                    flags[idx] = "duplicate"

            flags[i] = "distinct"
            unique_labelings.append(lbl)

        return unique_labelings

    def get_inequivalent_labelings(self):
        labelings = self.generate_possible_labelings()
        labelings = self.remove_duplicates(labelings)
        assert self.check_uniqueness(labelings)
        return labelings

    def check_uniqueness(self, labelings):
        for i, lbl in enumerate(labelings):
            for type_prm in itertools.permutations(range(self.num_type)):
                lbl_ex = [type_prm[e] for e in lbl]
                for prm in self.prm_all:
                    operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                    try:
                        idx = labelings.index(operated_lbl_ex)
                    except:
                        continue

                    if idx != i:
                        return False
        return True


class DerivativeStructure(object):

    def __init__(self, hnf, num_type, lattice_vectors, labeling):
        self.hnf = hnf
        self.num_type = num_type
        self.lattice_vectors = lattice_vectors
        self.labelings = labeling

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R

    def draw(self, ax=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for indices in itertools.product(*[range(e + 1) for e in self.hnf.diagonal().tolist()]):
            src = np.dot(self.lattice_vectors, np.array(indices))
            for i in range(3):
                directions = np.eye(3)
                dst = src + directions[:, i]
                tmp = np.concatenate([src, dst]).reshape(2, 3).T
                ax.plot(tmp[0], tmp[1], tmp[2])

        superlattice_vectors = np.dot(self.lattice_vectors, self.hnf)
        for i in range(3):
            origin = [0, 0, 0]
            ax.quiver(*origin, *superlattice_vectors[:, i].tolist(), arrow_length_ratio=0)


if __name__ == '__main__':
    hnf = np.array([[4]])
    num_type = 2
    labelings = Labeling(hnf, num_type)

    list_labelings = labelings.generate_possible_labelings()
    print('generate possible labelings')
    print(len(list_labelings))
    print(list_labelings)

    list_labelings = labelings.remove_duplicates(list_labelings)
    print('remove duplicates')
    print(len(list_labelings))
    print(list_labelings)
