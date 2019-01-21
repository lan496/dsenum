import unittest

from dsenum.superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)
from dsenum.permutation import Permutation
from dsenum.utils import get_lattice, get_symmetry_operations


class TestPermutation(unittest.TestCase):

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

    def test_translation_permutation(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = get_symmetry_operations(structure)
                rotations, translations = get_symmetry_operations(structure)

                list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                                              pl_rotations)
                for hnf in list_reduced_HNF:
                    permutation = Permutation(hnf, frac_coords.shape[0],
                                              frac_coords,
                                              rotations,
                                              translations)
                    self.assertTrue(self.validate_permutations(permutation.prm_t))

    def test_rigid_permutation(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = get_symmetry_operations(structure)
                rotations, translations = get_symmetry_operations(structure)

                list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                                              pl_rotations)
                for hnf in list_reduced_HNF:
                    permutation = Permutation(hnf, frac_coords.shape[0],
                                              frac_coords,
                                              rotations,
                                              translations)
                    self.assertTrue(self.validate_permutations(permutation.prm_rigid))

    def test_symmetry_permutation(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = get_symmetry_operations(structure)
                rotations, translations = get_symmetry_operations(structure)

                list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                                              pl_rotations)
                for hnf in list_reduced_HNF:
                    permutation = Permutation(hnf, frac_coords.shape[0],
                                              frac_coords,
                                              rotations,
                                              translations)
                    prm_all = permutation.get_symmetry_operation_permutaions()
                    self.assertTrue(self.validate_permutations(prm_all))

    def validate_permutations(self, permutations):
        # check if it is permutation
        for prm in permutations:
            if len(set(prm)) != len(prm):
                print(prm)
                print('not permutation')
                return False

        # check uniqueness
        if len(set(permutations)) != len(permutations):
            print('not unique')
            print(permutations)
            return False

        # TODO: add to check it is group

        return True


if __name__ == '__main__':
    unittest.main()
