import unittest

from tqdm import tqdm

from dsenum.enumerate import enumerate_derivative_structures

from dsenum.enumerate import enumerate_derivatives
from dsenum.coloring_generator import ColoringGenerator, FixedConcentrationColoringGenerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.utils import get_lattice
from dsenum.polya import polya_counting, polya_fixed_degrees_counting
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.coloring import SiteColoringEnumerator


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

    @unittest.skip
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

    @unittest.skip
    def test_colorings(self):
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

    @unittest.skip
    def test_colorings_with_polya(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            displacement_set = structure.frac_coords
            num_type = dct['num_type']
            num_sites_base = structure.num_sites

            for index, expected in zip(dct['indices'], dct['num_expected']):
                if index >= 6:
                    continue

                list_reduced_HNF, rotations, translations = \
                    generate_symmetry_distinct_superlattices(index, structure, return_symops=True)
                num_sites = num_sites_base * index
                cl_generator = ColoringGenerator(num_sites, num_type)

                for hnf in tqdm(list_reduced_HNF):
                    ds_permutaion = DerivativeStructurePermutation(hnf, displacement_set,
                                                                   rotations, translations)
                    sc_enum = SiteColoringEnumerator(num_type, ds_permutaion, cl_generator,
                                                     color_exchange=False,
                                                     leave_superperiodic=True,
                                                     use_all_colors=False)
                    colorings = sc_enum.unique_colorings()
                    cnt_polya = polya_counting(sc_enum.permutation_group, num_type)
                    self.assertEqual(len(colorings), cnt_polya)

    def test_fixed_colorings_with_polya(self):
        for name, dct in self.obj.items():
            structure = dct['structure']
            displacement_set = structure.frac_coords
            num_type = dct['num_type']
            num_sites_base = structure.num_sites

            for index, expected in zip(dct['indices'], dct['num_expected']):
                if index >= 6:
                    continue

                list_reduced_HNF, rotations, translations = \
                    generate_symmetry_distinct_superlattices(index, structure, return_symops=True)
                num_sites = num_sites_base * index
                color_ratio = [1] * (num_type - 1) + [num_sites - num_type + 1, ]
                cl_generator = FixedConcentrationColoringGenerator(num_sites, num_type, color_ratio)

                for hnf in tqdm(list_reduced_HNF):
                    ds_permutaion = DerivativeStructurePermutation(hnf, displacement_set,
                                                                   rotations, translations)
                    sc_enum = SiteColoringEnumerator(num_type, ds_permutaion, cl_generator,
                                                     color_exchange=False,
                                                     leave_superperiodic=True,
                                                     use_all_colors=False)
                    colorings = sc_enum.unique_colorings()
                    cnt_polya = polya_fixed_degrees_counting(sc_enum.permutation_group,
                                                             num_type, color_ratio)
                    self.assertEqual(len(colorings), cnt_polya)


if __name__ == '__main__':
    unittest.main()
