from functools import reduce
import itertools

from joblib import Parallel, delayed
import numpy as np

from dsenum.permutation import Permutation
from dsenum.utils import msp


class Labeling(object):
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    num_type: int
        # of kinds of atoms
    labelgen: LabelGenerator
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
    leave_superperiodic: bool
        if True, leave superperiodic labeling

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    num_site: int
        # of sites in unit cell of superlattice
    """
    def __init__(self, hnf, num_type, labelgen, num_site_parent=1, displacement_set=None,
                 rotations=None, translations=None,
                 ignore_site_property=False, leave_superperiodic=False):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())
        self.num_type = num_type
        self.ignore_site_property = ignore_site_property
        self.leave_superperiodic = leave_superperiodic

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

        self.labelgen = labelgen

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
        list_labelings = self.labelgen.list_labelings
        self.valid_flags = {self.convert_base(lbl): True for lbl in list_labelings}
        return list_labelings

    def remove_translation_and_exchanging_duplicates(self, labelings):
        unique_labelings = []
        for lbl in labelings:
            idx = self.convert_base(lbl)
            if not self.valid_flags[idx]:
                continue

            # remove superperiodic
            if not self.leave_superperiodic:
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


class LabelGenerator:

    def __init__(self, index, num_type, num_site_parent=1,
                 constraints=None, oxi_states=None,
                 force_unitcell_neutraliry=False,
                 n_jobs=-1):
        self.index = index
        self.num_type = num_type
        self.num_site_parent = num_site_parent
        self.num_site = self.num_site_parent * self.index
        self.force_unitcell_neutraliry = force_unitcell_neutraliry
        self.n_jobs = n_jobs

        self.constraints = constraints
        if oxi_states is None:
            self.oxi_states = None
        else:
            self.oxi_states = np.array(oxi_states)

        self.list_labelings = self.generate_possible_labelings()
        print("# of possible labelings: {}".format(len(self.list_labelings)))

    def _get_valid_ratios(self):
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

    def _is_valid_labeling(self, lbl):
        if len(set(lbl)) != self.num_type:
            return False

        if self.constraints is not None:
            # remove a labeling that does not satisfy the constraints for site preferences
            raveled = np.array(lbl).reshape(self.num_site_parent, self.index)
            for i in range(self.num_site_parent):
                if not set(raveled[i]) <= set(self.constraints[i]):
                    return False

        if self.force_unitcell_neutraliry and (self.oxi_states is not None):
            # only take labeling that satisfy neutrality in each primitive cell
            lbl_oxi_states = self.oxi_states[lbl].reshape(self.num_site_parent, self.index)
            if not np.allclose(np.sum(lbl_oxi_states, axis=0), np.zeros(self.index)):
                return False

        return True

    def _generate_with_ratio(self, ratio):
        items = [i for i in range(self.num_type) for r in range(ratio[i])]
        lbls_tmp = [lbl for lbl in msp(items) if self._is_valid_labeling(lbl)]
        print("# of possible labelings with ratio {}: {}".format(ratio, len(lbls_tmp)))
        return lbls_tmp

    def generate_possible_labelings(self):
        print("generate possible labelings")
        if self.oxi_states is not None:
            # restrict ratio by composition neutrality
            valid_ratios = self._get_valid_ratios()
            print("# of possible ratios: {}".format(len(valid_ratios)))
            list_labelings = Parallel(n_jobs=self.n_jobs, verbose=10)([
                delayed(self._generate_with_ratio)(ratio) for ratio in valid_ratios])
            list_labelings = list(itertools.chain.from_iterable(list_labelings))
        else:
            list_labelings = [lbl for lbl
                              in itertools.product(range(self.num_type), repeat=self.num_site)
                              if self._is_valid_labeling(lbl)]
        return list_labelings

    def convert_base(self, labeling):
        # labeling(i.e. list of int) -> int
        ret = reduce(lambda x, y: x * self.num_type + y, labeling)
        return ret


class ListBasedLabelGenerator:
    """
    generate possible labeling from a given list
    """
    def __init__(self, list_labelings):
        self.list_labelings = list_labelings
