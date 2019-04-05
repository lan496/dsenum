from itertools import permutations

from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.coloring_generator import ColoringGenerator


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
                 color_exchange=True, leave_superperiodic=False):
        self.num_color = num_color
        self.ds_permutation = ds_permutation
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        self.leave_superperiodic = leave_superperiodic

        self.permutation_group = self.ds_permutation.get_symmetry_operation_permutaions()

    @property
    def translation_permutations(self):
        return self.ds_permutation.prm_t

    def unique_colorings(self):
        self.colorings, self.flags = self.cl_generator.generate_all_colorings()
        self.remove_symmetry_duplicates()
        if not self.leave_superperiodic:
            self.remove_superperiodic()

        colorings = [cl for cl_hash, cl in self.colorings.items() if self.flags[cl_hash]]
        return colorings

    def remove_superperiodic(self):
        for cl in self.colorings:
            cl_hash = self.cl_generator.hash_coloring(cl)
            # avoid already-visited coloring
            if not self.flags[cl_hash]:
                continue

            # assume self.translation_permutations[0] is identity
            if any([self.product_permutation(prm, cl) == cl]
                   for prm in self.translation_permutations[1:]):
                self.flags[cl_hash] = False

    def remove_rotation_duplicates(self):
        for cl in self.colorings:
            cl_hash = self.cl_generator.hash_coloring(cl)
            # avoid already-visited coloring
            if not self.flags[cl_hash]:
                continue

            # assume self.permutation_group[0] is identity
            for prm in self.permutation_group[1:]:
                acted_cl = product_permutations(prm, cl)
                acted_cl_hash = self.cl_generator.hash_coloring(acted_cl)
                self.flags[acted_cl_hash] = False

            if self.color_exchange:
                for cl_prm in permutations(range(self.num_color)):
                    exchanged_cl = [cl_prm[c] for c in cl]
                    exchanged_cl_hash = self.cl_generator.hash_coloring(exchanged_cl)
                    if exchanged_cl_hash == cl_hash:
                        continue
                    # assume self.permutation_group[0] is identity
                    for prm in self.permutation_group[1:]:
                        acted_cl = product_permutations(prm, exchanged_cl)
                        acted_cl_hash = self.cl_generator.hash_coloring(acted_cl)
                        self.flags[acted_cl_hash] = False


def product_permutations(p1, p2):
    """
    (p1 p2)(i) = p1(p2(i))
    """
    perm = [p1[p2[i]] for i in range(len(p1))]
    return perm
