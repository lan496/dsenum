import itertools

import numpy as np

from permutation import Permutation


class Labeling(object):

    def __init__(self, hnf, num_type, list_rotations=None):
        self.hnf = hnf
        self.num_type = num_type
        self.num_site = np.prod(self.hnf.diagonal())

        self.permutation = Permutation(self.hnf)
        self.list_prm = self.permutation.get_permutation_list()

        self.list_rotations = list_rotations
        self.list_rotations_prm = self.per

    def act_permutation(self, labeling, prm):
        ret = tuple([labeling[prm[i]] for i in range(self.num_site)])
        return ret

    def generate_possible_labelings(self):
        list_labelings = []

        for labeling in itertools.product(range(self.num_type), repeat=self.num_site):
            if len(set(labeling)) == self.num_type:
                list_labelings.append(labeling)

        return list_labelings

    def remove_translation_symmetry(self, list_labelings):
        # remove translation and super-periodic duplicates
        flags = [True for _ in list_labelings]
        list_reduced_labelings = []
        for i, labeling in enumerate(list_labelings):
            if not flags[i]:
                continue

            is_superperiodic = False
            # assuming that self.list_prm[0] is identity operation
            for prm in self.list_prm[1:]:
                prm_labeling = self.act_permutation(labeling, prm)
                if prm_labeling == labeling:
                    is_superperiodic = True
                    break

                idx = list_labelings.index(prm_labeling)
                flags[idx] = False

            if is_superperiodic:
                continue

            list_reduced_labelings.append(labeling)

        return list_reduced_labelings

    def remove_label_exchange_symmetry(self, list_labelings):
        flags = [True for _ in list_labelings]
        list_reduced_labelings = []
        for i, labeling in enumerate(list_labelings):
            if not flags[i]:
                continue
            list_reduced_labelings.append(labeling)

            for type_prm in itertools.permutations(range(self.num_site)):
                if type_prm == range(self.num_site):
                    continue

                labeling_ex = [type_prm[e] for e in labeling]
                for prm in self.list_prm[1:]:
                    labeling_ex_prm = self.act_permutation(labeling_ex,
                                                           prm)
                    try:
                        idx = list_labelings.index(labeling_ex_prm)
                        flags[idx] = False
                    except ValueError:
                        continue

        return list_reduced_labelings

    def remove_rotation_symmetry(self, list_labelings):
        flags = [True for _ in list_labelings]
        list_reduced_labelings = []
        for i, labeling in enumerate(list_labelings):


if __name__ == '__main__':
    hnf = np.array([[4]])
    num_type = 2
    labelings = Labeling(hnf, num_type)
    list_labelings = labelings.generate_possible_labelings()
    print(len(list_labelings))
    print(list_labelings)
    list_labelings = labelings.remove_translation_symmetry(list_labelings)
    print(len(list_labelings))
    print(list_labelings)
    list_labelings = labelings.remove_label_exchange_symmetry(list_labelings)
    print(len(list_labelings))
    print(list_labelings)

