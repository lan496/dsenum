from functools import reduce
import itertools

import numpy as np

from derivative.permutation import Permutation
from derivative.utils import msp


class Labeling(object):
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    num_type: int
        # of kinds of atoms
    num_site_parent: int
        # of atoms in parent multilattice
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates of A in multilattice site
    rotations: array, (# of symmetry operations, dim, dim)
        rotations in fractinal coordinations of A
    translations: array, (# of symmetry operations, dim)
        translations in fractinal coordinations of A
    ignore_site_property: bool
        if True, treat label-exchange duplicates as distinct

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    num_site: int
        # of sites in unit cell of superlattice
    """
    def __init__(self, hnf, num_type, num_site_parent=1, displacement_set=None,
                 rotations=None, translations=None, ignore_site_property=False):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())
        self.num_type = num_type
        self.ignore_site_property = ignore_site_property

        self.num_site_parent = num_site_parent
        if self.num_site_parent == 1:
            self.displacement_set = np.array([[0, 0, 0]])
        else:
            self.displacement_set = displacement_set
            assert self.displacement_set.shape[0] == self.num_site_parent

        self.num_site = self.index * self.num_site_parent

        self.permutation = Permutation(self.hnf, self.num_site_parent, self.displacement_set,
                                       rotations, translations)
        # assuming that the 0-th element of permutaions is identity operation
        self.prm_all = self.permutation.get_symmetry_operation_permutaions()

    @property
    def prm_t(self):
        return self.permutation.prm_t

    @property
    def rotaions(self):
        return self.permutation.rotations

    @property
    def translations(self):
        return self.permutation.translations

    def act_permutation(self, labeling, prm):
        ret = tuple([labeling[prm[i]] for i in range(self.num_site)])
        return ret

    def convert_base(self, labeling):
        # labeling(i.e. list of int) -> int
        ret = reduce(lambda x, y: x * self.num_type + y, labeling)
        return ret

    def generate_possible_labelings(self):
        list_labelings = [lbl for lbl
                          in itertools.product(range(self.num_type), repeat=self.num_site)
                          if len(set(lbl)) == self.num_type]
        self.valid_flags = [True for _ in range(self.num_type ** self.num_site)]
        return list_labelings

    def remove_translation_and_exchanging_duplicates(self, labelings):
        unique_labelings = []
        for lbl in labelings:
            idx = self.convert_base(lbl)
            if not self.valid_flags[idx]:
                continue

            # remove superperiodic
            if any([(self.act_permutation(lbl, prm) == lbl) for prm in self.prm_t[1:]]):
                self.valid_flags[idx] = False
                continue

            # remove translation and exchanging duplicates
            if self.ignore_site_property:
                for prm in self.prm_t:
                    if prm == range(self.num_site):
                        continue
                    operated_lbl = self.act_permutation(lbl, prm)
                    idx_lbl = self.convert_base(operated_lbl)
                    self.valid_flags[idx_lbl] = False
            else:
                for type_prm in itertools.permutations(range(self.num_type)):
                    lbl_ex = [type_prm[e] for e in lbl]
                    for prm in self.prm_t:
                        if (lbl_ex == lbl) and (prm == range(self.num_site)):
                            continue
                        operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                        idx_lbl_ex = self.convert_base(operated_lbl_ex)
                        self.valid_flags[idx_lbl_ex] = False

            self.valid_flags[idx] = True
            unique_labelings.append(lbl)
        return unique_labelings

    def remove_rotation_duplicates(self, labelings):
        unique_labelings = []
        for lbl in labelings:
            idx = self.convert_base(lbl)
            if not self.valid_flags[idx]:
                continue

            if self.ignore_site_property:
                for prm in self.prm_all:
                    if prm == range(self.num_site):
                        continue
                    operated_lbl = self.act_permutation(lbl, prm)
                    idx_lbl = self.convert_base(operated_lbl)
                    self.valid_flags[idx_lbl] = False
            else:
                for type_prm in itertools.permutations(range(self.num_type)):
                    lbl_ex = [type_prm[e] for e in lbl]
                    for prm in self.prm_all:
                        if (lbl_ex == lbl) and (prm == range(self.num_site)):
                            continue
                        operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                        idx_lbl_ex = self.convert_base(operated_lbl_ex)
                        self.valid_flags[idx_lbl_ex] = False

            self.valid_flags[idx] = True
            unique_labelings.append(lbl)

        return unique_labelings

    def remove_duplicates(self, labelings):
        unique_labelings = self.remove_translation_and_exchanging_duplicates(labelings)
        unique_labelings = self.remove_rotation_duplicates(unique_labelings)
        return unique_labelings

    def get_inequivalent_labelings(self):
        labelings = self.generate_possible_labelings()
        labelings = self.remove_duplicates(labelings)
        return labelings

    def check_uniqueness(self, labelings):
        for i, lbl in enumerate(labelings):
            for type_prm in itertools.permutations(range(self.num_type)):
                lbl_ex = [type_prm[e] for e in lbl]
                for prm in self.prm_all:
                    operated_lbl_ex = self.act_permutation(lbl_ex, prm)
                    try:
                        idx = labelings.index(operated_lbl_ex)
                    except ValueError:
                        continue

                    if idx != i:
                        return False
        return True


class ConstraintedLabeling(Labeling):
    """
    Attributes
    ----------
    constraints: [[0, 1], [0, 2, 3], [0, 2, 3], ...]
    oxi_states: e.g. [0, -1, 2, 4]
    """
    def __init__(self, hnf, num_type, num_site_parent=1, displacement_set=None,
                 rotations=None, translations=None, ignore_site_property=False,
                 constraints=None, oxi_states=None):
        super().__init__(hnf, num_type, num_site_parent, displacement_set,
                         rotations, translations, ignore_site_property)
        self.constraints = constraints
        self.oxi_states = np.array(oxi_states)

    def get_valid_ratios(self):
        ratios = []

        def _dfs(ratio, cnt):
            if (len(ratio) == self.num_type) and (cnt == self.num_site):
                if self.oxi_states is None:
                    ratios.append(ratio)
                elif np.sum(np.array(ratio) * self.oxi_states) == 0:
                    ratios.append(ratio)
            else:
                for i in range(1, self.num_site - cnt + 1):
                    tmp = ratio + [i, ]
                    _dfs(tmp, cnt + i)

        _dfs([], 0)
        return ratios

    def is_valid_labeling(self, lbl):
        raveled = np.array(lbl).reshape(-1)
        for i in range(self.num_site_parent):
            if not set(raveled[i].reshape(-1)) < set(self.constraints[i]):
                return False
        return True

    def generate_possible_labelings(self):
        valid_ratios = self.get_valid_ratios()
        print("# of possible ratios: {}".format(len(valid_ratios)))

        list_labelings = []
        for ratio in valid_ratios:
            items = [i for i in range(self.num_type) for r in range(ratio[i])]
            lbls_tmp = [lbl for lbl in msp(items) if self.is_valid_labeling(lbl)]
            list_labelings.extend(lbls_tmp)
        print("# of possible labelings: {}".format(len(list_labelings)))

        self.valid_flags = {self.convert_base(lbl): True for lbl in list_labelings}

        """
        self.valid_flags = [True for _ in range(self.num_type ** self.num_site)]
        gen = itertools.product(*[itertools.product(cns, repeat=self.index)
                                  for cns in self.constraints])
        list_labelings = [np.array(e).reshape(-1).tolist() for e in gen]
        list_labelings = [lbl for lbl in list_labelings if len(set(lbl)) == self.num_type]
        """
        return list_labelings


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
