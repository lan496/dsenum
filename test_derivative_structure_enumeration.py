import unittest

import numpy as np

from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)
from smith_normal_form import smith_normal_form
from permutation import Permutation
from derivative_structure import SuperMultilattice
from utils import get_lattice, get_symmetry_operations


class TestEnumerateSuperlattice(unittest.TestCase):

    def test_generate_all_superlattices(self):
        # https://oeis.org/A001001
        num_expected = [1, 7, 13, 35, 31, 91, 57, 155, 130, 217,
                        133, 455, 183, 399, 403, 651, 307, 910, 381, 1085,
                        741, 931, 553, 2015, 806, 1281, 1210, 1995, 871, 2821,
                        993, 2667, 1729, 2149, 1767, 4550, 1407, 2667, 2379, 4805,
                        1723, 5187, 1893, 4655, 4030, 3871, 2257, 8463, 2850, 5642,
                        3991, 6405, 2863]
        max_index = len(num_expected)

        for index, expected in zip(range(1, max_index + 1), num_expected):
            list_HNF = generate_all_superlattices(index)
            self.assertEqual(len(list_HNF), expected)

    def test_reduce_HNF_list_by_parent_lattice_symmetry(self):
        # confirm table 4
        obj = {
            'fcc': {
                'structure': get_lattice('fcc'),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
            },
            'bcc': {
                'structure': get_lattice('bcc'),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
            },
            'sc': {
                'structure': get_lattice('sc'),
                'num_expected': [1, 3, 3, 9, 5, 13, 7, 24, 14, 23]
            },
            'hex': {
                'structure': get_lattice('hex'),
                'num_expected': [1, 3, 5, 11, 7, 19, 11, 34, 23, 33]
            },
            'tetragonal': {
                'structure': get_lattice('tet'),
                'num_expected': [1, 5, 5, 17, 9, 29, 13, 51, 28, 53]
            },
            'hcp': {
                'structure': get_lattice('hcp'),
                'num_expected': [1, 3, 5, 11, 7, 19, 11, 34, 23, 33]
            }
        }

        for name, dct in obj.items():
            print('#' * 40)
            structure = dct['structure']
            for index, expected in zip(range(1, len(dct['num_expected']) + 1), dct['num_expected']):
                list_HNF = generate_all_superlattices(index)
                rotations, _ = get_symmetry_operations(structure)

                list_reduced_HNF = \
                    reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                               rotations)
                print('{}, index {}: superlattices {} {}'.format(name, index,
                                                                 len(list_reduced_HNF),
                                                                 expected))
                self.assertEqual(len(list_reduced_HNF), expected)


class TestSmithNormalForm(unittest.TestCase):

    def test_smf(self):
        list_matrix = [
            np.array([
                [2, 0],
                [1, 4]
            ]),
            np.array([
                [2, 4, 4],
                [-6, 6, 12],
                [10, -4, -16]
            ]),
            np.array([
                [8, 4, 8],
                [4, 8, 4]
            ]),
            np.array([
                [-6, 111, -36, 6],
                [5, -672, 210, 74],
                [0, -255, 81, 24],
                [-7, 255, -81, -10]
            ]),
            np.array([
                [3, -1, -1],
                [-1, 3, -1],
                [-1, -1, 3]
            ]),
            np.array([
                [1, 0, 0],
                [1, 2, 0],
                [0, 0, 2]
            ]),
        ]
        list_expected = [
            np.diag([1, 8]),
            np.diag([2, 6, 12]),
            np.array([
                [4, 0, 0],
                [0, 12, 0]
            ]),
            np.diag([1, 3, 21, 0]),
            np.diag([1, 4, 4]),
            np.diag([1, 2, 2])
        ]

        for M, expected in zip(list_matrix, list_expected):
            D, L, R = smith_normal_form(M)
            D_re = np.dot(L, np.dot(M, R))
            self.assertAlmostEqual(np.linalg.det(L) ** 2, 1)
            self.assertAlmostEqual(np.linalg.det(R) ** 2, 1)
            self.assertTrue(np.array_equal(D_re, D))

    def test_number_of_snf(self):
        # confirm table-3
        num_hnf_expected = [1, 7, 13, 35, 31, 91, 57, 155, 130, 217,
                            133, 455, 183, 399, 403, 651]
        num_snf_expected = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1,
                            1, 2, 1, 1, 1, 4]
        max_index = len(num_hnf_expected)

        for index, hnf_expected, snf_expected in zip(range(1, max_index + 1), num_hnf_expected, num_snf_expected):
            list_HNF = generate_all_superlattices(index)
            self.assertEqual(len(list_HNF), hnf_expected)

            list_SNF = set()
            for hnf in list_HNF:
                snf, _, _ = smith_normal_form(hnf)
                dag = tuple(snf.diagonal())
                list_SNF.add(dag)

            self.assertEqual(len(list_SNF), snf_expected)


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
            print('*' * 40)
            print(name)
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                print('    index={}'.format(index),)
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = \
                    get_symmetry_operations(structure, parent_lattice=True)
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
            print('*' * 40)
            print(name)
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                print('    index={}'.format(index),)
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = \
                    get_symmetry_operations(structure, parent_lattice=True)
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
            print('*' * 40)
            print(name)
            structure = dct['structure']
            frac_coords = structure.frac_coords
            for index in dct['indices']:
                print('    index={}'.format(index),)
                list_HNF = generate_all_superlattices(index)
                pl_rotations, pl_translations = \
                    get_symmetry_operations(structure, parent_lattice=True)
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


class TestSmall(unittest.TestCase):

    def test(self):
        parent_lattice = get_lattice('fcc')
        frac_coords = parent_lattice.frac_coords
        A = parent_lattice.lattice.matrix.T
        rotations, translations = get_symmetry_operations(parent_lattice)

        hnf = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 3]])

        permutation = Permutation(hnf, frac_coords.shape[0], frac_coords,
                                  rotations, translations)
        actual = permutation.rotations

        sl = SuperMultilattice(hnf, A, frac_coords.shape[0], frac_coords)
        sl_rotations, sl_translations = get_symmetry_operations(sl.get_structure())
        print(len(rotations), len(actual), len(sl_rotations))
        # import pdb; pdb.set_trace()
        # self.assertEqual(len(actual), len(sl_rotations))


if __name__ == '__main__':
    unittest.main()
