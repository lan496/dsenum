from itertools import permutations

from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.coloring_generator import ColoringGenerator
from dsenum.core import hash_in_all_configuration, act_permutation


class DirectColoringEnumerator:
    """
    Parameters
    ----------
    permutation_group: list of permutation
    num_color: int
    cl_generator: ColoringGenerator
    color_exchange: bool
    """
    def __init__(self, permutation_group, num_color, cl_generator: ColoringGenerator,
                 color_exchange=True):
        self.permutation_group = permutation_group
        self.num_color = num_color
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange

        self.list_colorings, self.flags = self.cl_generator.generate_all_colorings()

    def _hash(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)

    def _walk_orbit(self, coloring, include_identity=False):
        offset = 0 if include_identity else 1
        # assume self.permutation_group[0] is identity
        for prm in self.permutation_group[offset:]:
            acted_cl = act_permutation(prm, coloring)
            acted_cl_hash = self._hash(acted_cl)
            self.flags[acted_cl_hash] = False

    def coset_enumerate(self):
        colorings = []

        for cl in self.list_colorings:
            cl_hash = self._hash(cl)
            # avoid already-visited coloring
            if not self.flags[cl_hash]:
                continue
            colorings.append(cl)

            self._walk_orbit(cl, include_identity=False)

            if self.color_exchange:
                for cl_prm in permutations(range(self.num_color)):
                    exchanged_cl = [cl_prm[c] for c in cl]
                    exchanged_cl_hash = self._hash(exchanged_cl)
                    if exchanged_cl_hash == cl_hash:
                        continue
                    self._walk_orbit(exchanged_cl, include_identity=True)

        return colorings


class SiteColoringEnumerator(object):
    """
    Parameters
    ----------
    ds_permutation: DerivativeStructurePermutation object
    num_color: int
    color_exchange: bool
    leave_superperiodic: bool
    """

    def __init__(self, num_color,
                 ds_permutation: DerivativeStructurePermutation,
                 cl_generator: ColoringGenerator,
                 color_exchange=True, leave_superperiodic=False, use_all_colors=True):
        self.num_color = num_color
        self.ds_permutation = ds_permutation
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        self.leave_superperiodic = leave_superperiodic
        self.use_all_colors = use_all_colors

        """
        if self.ds_permutation.num_site < self.num_color:
            raise ValueError('too many num_types')
        """

        self.permutation_group = self.ds_permutation.get_symmetry_operation_permutaions()

        self.clenum = DirectColoringEnumerator(self.permutation_group, self.num_color,
                                               self.cl_generator, color_exchange=color_exchange)

    @property
    def translation_permutations(self):
        return self.ds_permutation.prm_t

    def _hash(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)

    def unique_colorings(self):
        symmetry_uniqued_coloring = self.clenum.coset_enumerate()
        colorings = []

        for cl in symmetry_uniqued_coloring:
            if self.use_all_colors and (not self._has_all_colors(cl)):
                continue
            if (not self.leave_superperiodic) and self._is_superperiodic(cl):
                continue
            colorings.append(cl)

        return colorings

    def _has_all_colors(self, coloring):
        if len(set(coloring)) == self.num_color:
            return True
        else:
            return False

    def _is_superperiodic(self, coloring):
        # assume self.translation_permutations[0] is identity
        if any([self._hash(act_permutation(prm, coloring)) == self._hash(coloring)
               for prm in self.translation_permutations[1:]]):
            return True
        else:
            return False
