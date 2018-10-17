import unittest

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from derivative_structure_enumeration import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)

from smith_normal_form import smith_normal_form


class TestDerivativeStructureEnumeration(unittest.TestCase):

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

    def test_reduce_HNF_list_by_parent_lattice_symmetry_fcc_bcc(self):
        # https://oeis.org/A045790
        lst_num = [1, 2, 3, 7, 5, 10, 7, 20, 14, 18,
                   11, 41, 15, 28, 31, 58, 21, 60, 25, 77,
                   49, 54, 33, 144, 50, 72, 75, 123, 49, 158,
                   55, 177, 97, 112, 99, 268, 75, 136, 129, 286,
                   89, 268, 97, 249, 218, 190, 113, 496, 146, 280]
        obj = {
            'fcc': {
                'structure': self.get_face_centered_cubic(),
                'num_expected': lst_num
            },
            'bcc': {
                'structure': self.get_body_centered_cubic(),
                'num_expected': lst_num
            }
        }
        obj = {}

        for name, dct in obj.items():
            print('#' * 40)
            for index, expected in zip(range(1, len(dct['num_expected']) + 1), dct['num_expected']):
                list_HNF = generate_all_superlattices(index)
                sp = SpacegroupAnalyzer(dct['structure'])
                list_rotation_matrix = sp.get_symmetry_dataset()['rotations']

                list_reduced_HNF = \
                    reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                               list_rotation_matrix,
                                                               n_jobs=-1)
                print('{}, index {}: superlattices {} {}'.format(name, index,
                                                                 len(list_reduced_HNF),
                                                                 expected))
                self.assertEqual(len(list_reduced_HNF), expected)

    def test_reduce_HNF_list_by_parent_lattice_symmetry(self):
        # confirm table 4
        obj = {
            'fcc': {
                'structure': self.get_face_centered_cubic(),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
            },
            'bcc': {
                'structure': self.get_body_centered_cubic(),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]
            },
            'sc': {
                'structure': self.get_simple_cubic(),
                'num_expected': [1, 3, 3, 9, 5, 13, 7, 24, 14, 23]
            },
            'hex': {
                'structure': self.get_hexagonal(),
                'num_expected': [1, 3, 5, 11, 7, 19, 11, 34, 23, 33]
            },
            'tetragonal': {
                'structure': self.get_tetragonal(),
                'num_expected': [1, 5, 5, 17, 9, 29, 13, 51, 28, 53]
            }
        }

        for name, dct in obj.items():
            print('#' * 40)
            for index, expected in zip(range(1, len(dct['num_expected']) + 1), dct['num_expected']):
                list_HNF = generate_all_superlattices(index)
                sp = SpacegroupAnalyzer(dct['structure'])
                list_rotation_matrix = sp.get_symmetry_dataset()['rotations']

                from time import time
                start = time()
                list_reduced_HNF = \
                    reduce_HNF_list_by_parent_lattice_symmetry(list_HNF,
                                                               list_rotation_matrix)
                print('{}, index {}: superlattices {} {}'.format(name, index,
                                                                 len(list_reduced_HNF),
                                                                 expected))
                print(time() - start, 'sec')
                self.assertEqual(len(list_reduced_HNF), expected)

    def get_simple_cubic(self):
        latt = Lattice(np.eye(3))
        struct = Structure(latt, ['Po'], [[0, 0, 0]])
        return struct

    def get_face_centered_cubic(self):
        latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
        struct = Structure(latt, ['Al'], [[0, 0, 0]])
        return struct

    def get_body_centered_cubic(self):
        latt = Lattice(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]))
        struct = Structure(latt, ['Fe'], [[0, 0, 0]])
        return struct

    def get_hexagonal(self):
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
        struct = Structure(latt, ['Zn'], [[0, 0, 0]])
        return struct

    def get_tetragonal(self):
        latt = Lattice(np.diag([1, 1, 1.2]))
        struct = Structure(latt, ['Po'], [[0, 0, 0]])
        return struct


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


if __name__ == '__main__':
    unittest.main()
