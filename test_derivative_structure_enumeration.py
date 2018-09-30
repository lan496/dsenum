import unittest

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from derivative_structure_enumeration import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)


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

    def test_reduce_HNF_list_by_parent_lattice_symmetry(self):
        obj = {
            'fcc': {
                'structure': self.get_face_centered_cubic(),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18,
                                 11, 41, 15, 28, 31, 58, 21, 60, 25, 77,
                                 49, 54, 33, 144, 50, 72, 75, 123, 49, 158,
                                 55, 177, 97, 112, 99, 268, 75, 136, 129, 286,
                                 89, 268, 97, 249, 218, 190, 113, 496, 146, 280]
            },
            'bcc': {
                'structure': self.get_body_centered_cubic(),
                'num_expected': [1, 2, 3, 7, 5, 10, 7, 20, 14, 18,
                                 11, 41, 15, 28, 31, 58, 21, 60, 25, 77,
                                 49, 54, 33, 144, 50, 72, 75, 123, 49, 158,
                                 55, 177, 97, 112, 99, 268, 75, 136, 129, 286,
                                 89, 268, 97, 249, 218, 190, 113, 496, 146, 280]
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
            for index, expected in zip(range(1, len(dct['num_expected'])), dct['num_expected']):
                list_HNF = generate_all_superlattices(index)
                sp = SpacegroupAnalyzer(dct['structure'])
                list_rotation_matrix = sp.get_symmetry_dataset()['rotations']

                list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, list_rotation_matrix)
                print('{}, index {}: superlattices {} {}'.format(name, index,
                                                              len(list_reduced_HNF),
                                                                 expected))
                self.assertEqual(len(list_reduced_HNF), expected)

    def get_simple_cubic(self):
        latt = Lattice(np.eye(3))
        struct = Structure(latt, ['Po'], [[0, 0, 0]])
        return struct

    def get_face_centered_cubic(self):
        latt = Lattice(np.eye(3))
        struct = Structure(latt, ['Al'] * 4,
                           [[0, 0, 0], [0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]])
        return struct

    def get_body_centered_cubic(self):
        latt = Lattice(np.eye(3))
        struct = Structure(latt, ['Fe'] * 2, [[0, 0, 0], [0.5, 0.5, 0.5]])
        return struct

    def get_hexagonal(self):
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
        struct = Structure(latt, ['Zn'] * 2, [[0, 0, 0], [1 / 3, 1 / 3, 0.5]])
        return struct

    def get_tetragonal(self):
        latt = Lattice(np.diag([1, 1, 1.2]))
        struct = Structure(latt, ['Po'], [[0, 0, 0]])
        return struct


if __name__ == '__main__':
    unittest.main()
