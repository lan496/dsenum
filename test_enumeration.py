import unittest

from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry,
)
from labeling import Labeling
from derivative_structure import get_lattice


class TestUniqueLabeling(unittest.TestCase):

    def setUp(self):
        self.obj = {
            'fcc': {
                'structure': get_lattice('fcc'),
                'num_type': 2,
                'indices': range(1, 23 + 1),
                'num_expected': [0, 2, 3, 12, 14, 50, 52, 229, 252, 685,
                                 682, 3875, 2624, 9628, 16584,
                                 49764, 42135, 212612, 174104, 867893,
                                 1120708, 2628180, 3042732]
            },
            'sc': {
                'structure': get_lattice('sc'),
                'num_type': 2,
                'indices': range(1, 4 + 1),
                'num_expected': [0, 3, 3, 15]
            },
            'fcc_ternary': {
                'structure': get_lattice('fcc'),
                'num_type': 3,
                'indices': range(1, 10 + 1),
                'num_expected': [0, 0, 3, 13, 23, 130, 197, 1267, 2322, 9332]
            },
            'fcc_quaternary': {
                'structure': get_lattice('fcc'),
                'num_type': 4,
                'indices': range(1, 10 + 1),
                'num_expected': [0, 0, 0, 7, 9, 110, 211, 2110, 5471, 32362]
            },
            'hcp': {
                'structure': get_lattice('hcp'),
                'num_type': 2,
                'indices': range(1, 10 + 1),
                'num_expected': [0, 7, 30, 163, 366,
                                 2613, 5268, 42901, 119528, 662193]
            }
        }

    def test_labelings(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            num_type = dct['num_type']
            for index, expected in zip(dct['indices'], dct['num_expected']):
                if index > 8:
                    continue
                list_HNF = generate_all_superlattices(index)
                sym_dataset = SpacegroupAnalyzer(structure)\
                    .get_symmetry_dataset()
                rotations = sym_dataset['rotations']
                list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)

                lbls = []
                for hnf in list_reduced_HNF:
                    labeling = Labeling(hnf, num_type, rotations)
                    lbls_tmp = labeling.get_inequivalent_labelings()
                    lbls.extend(lbls_tmp)

                print('{}, index {}, labelings {} (expected {})'.format(name, index,
                                                                        len(lbls), expected))
                self.assertEqual(len(lbls), expected)


if __name__ == '__main__':
    unittest.main()
