import unittest

from dsenum.enumerate import enumerate_derivative_structures

from dsenum.enumerate import enumerate_derivatives
from dsenum.utils import get_lattice


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
                'num_expected': [1, 7, 30, 163, 366,
                                 2613, 5268, 42901, 119528, 662193]
            }
        }

    def test_labelings(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            # displacement_set = structure.frac_coords
            num_type = dct['num_type']
            for index, expected in zip(dct['indices'], dct['num_expected']):
                if index > 8:
                    continue
                actual = enumerate_derivative_structures(structure,
                                                         index,
                                                         num_type)
                self.assertEqual(len(actual), expected)


class TestUniqueColoring(unittest.TestCase):

    def setUp(self):
        self.obj = {
            'fcc': {
                'structure': get_lattice('fcc'),
                'num_type': 2,
                'indices': range(2, 23 + 1),
                'num_expected': [2, 3, 12, 14, 50, 52, 229, 252, 685,
                                 682, 3875, 2624, 9628, 16584,
                                 49764, 42135, 212612, 174104, 867893,
                                 1120708, 2628180, 3042732]
            },
            'sc': {
                'structure': get_lattice('sc'),
                'num_type': 2,
                'indices': range(2, 4 + 1),
                'num_expected': [3, 3, 15]
            },
            'fcc_ternary': {
                'structure': get_lattice('fcc'),
                'num_type': 3,
                'indices': range(3, 10 + 1),
                'num_expected': [3, 13, 23, 130, 197, 1267, 2322, 9332]
            },
            'fcc_quaternary': {
                'structure': get_lattice('fcc'),
                'num_type': 4,
                'indices': range(4, 10 + 1),
                'num_expected': [7, 9, 110, 211, 2110, 5471, 32362]
            },
            'hcp': {
                'structure': get_lattice('hcp'),
                'num_type': 2,
                'indices': range(2, 10 + 1),
                'num_expected': [7, 30, 163, 366,
                                 2613, 5268, 42901, 119528, 662193]
            }
        }

    def test_labelings(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            # displacement_set = structure.frac_coords
            num_type = dct['num_type']
            for index, expected in zip(dct['indices'], dct['num_expected']):
                if index > 8:
                    continue
                actual = enumerate_derivatives(structure, index, num_type,
                                               color_exchange=True,
                                               leave_superperiodic=False)
                self.assertEqual(len(actual), expected)


if __name__ == '__main__':
    unittest.main()
