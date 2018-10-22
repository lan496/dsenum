import itertools

import numpy as np

from permutation import Permutation


class Labeling(object):

    def __init__(self, hnf, num_type, list_rotations=None):
        self.hnf = hnf
        self.num_type = num_type
        self.num_site = np.prod(self.hnf.diagonal())

        self.permutation = Permutation(self.hnf)
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
        flags = [True for _ in labelings]
        for i, lbl in enumerate(labelings):
            if not flags[i]:
                continue

            # remove superperiodic
            if any([(self.act_permutation(lbl, prm) == lbl) for prm in self.prm_t[1:]]):
                flags[i] = False
                continue

            # remove symmetry operation and exchanging duplicates
            for type_prm in itertools.permutations(range(self.num_type)):
                lbl_ex = [type_prm[e] for e in lbl]
                for prm in self.prm_all:
                    if (lbl_ex == lbl) and (prm == range(self.num_site)):
                        continue
                    operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                    idx = labelings.index(operated_lbl_ex)
                    flags[idx] = False

            flags[i] = True
            unique_labelings.append(lbl)

        return unique_labelings


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
