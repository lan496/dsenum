import unittest

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)


class TestDerivativeStructureEnumeration(unittest.TestCase):

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

    def get_face_centered_cubic(self):
        latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
        struct = Structure(latt, ['Al'], [[0, 0, 0]])
        return struct

    def get_body_centered_cubic(self):
        latt = Lattice(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]))
        struct = Structure(latt, ['Fe'], [[0, 0, 0]])
        return struct


if __name__ == '__main__':
    unittest.main()
