import unittest

from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.utils import get_lattice

from dsenum.permutation_group import DerivativeStructurePermutation, is_permutation_group


class TestPermutationGroup(unittest.TestCase):

    def setUp(self):
        self.obj = {
            'fcc': {
                'structure': get_lattice('fcc'),
                'num_type': 2,
                'indices': range(1, 10 + 1),
            },
            'bcc': {
                'structure': get_lattice('bcc'),
                'indices': range(1, 10 + 1)
            },
            'sc': {
                'structure': get_lattice('sc'),
                'indices': range(1, 10 + 1),
            },
            'hex': {
                'structure': get_lattice('hex'),
                'indices': range(1, 10 + 1),
            },
            'tetragonal': {
                'structure': get_lattice('tet'),
                'indices': range(1, 10 + 1),
            },
            'hcp': {
                'structure': get_lattice('hcp'),
                'indices': range(1, 10 + 1)
            }
        }

    def test_permutations(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                list_reduced_HNF, rotations, translations = \
                    generate_symmetry_distinct_superlattices(index, structure, return_symops=True)
                for hnf in list_reduced_HNF:
                    dsperm = DerivativeStructurePermutation(hnf, frac_coords,
                                                            rotations, translations)
                    self.assertTrue(is_permutation_group(dsperm.prm_t))
                    prm_all = dsperm.get_symmetry_operation_permutations()
                    self.assertTrue(is_permutation_group(prm_all))


if __name__ == '__main__':
    unittest.main()
