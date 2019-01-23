import unittest

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Specie

from dsenum.enumerate import enumerate_derivative_structures
from dsenum.utils import get_lattice


class TestDerivativeStructure(unittest.TestCase):

    def test_fcc(self):
        ps = get_lattice('fcc')
        max_index = 4
        num_type = 2
        species = [Specie('Cu'), Specie('Au')]
        list_num_exp = [0, 2, 3, 12]

        self.check(ps, max_index, num_type, species, list_num_exp)

    def test_sc(self):
        ps = get_lattice('sc')
        max_index = 4
        num_type = 2
        species = [Specie('Cu'), Specie('Au')]
        list_num_exp = [0, 3, 3, 15]

        self.check(ps, max_index, num_type, species, list_num_exp)

    def test_hcp(self):
        ps = get_lattice('hcp')
        max_index = 4
        num_type = 2
        species = [Specie('Cu'), Specie('Au')]
        list_num_exp = [1, 7, 30]
        # inconsistent with StructureMatcher in index=4
        # list_num_exp = [1, 7, 30, 163]

        self.check(ps, max_index, num_type, species, list_num_exp)

    def check(self, primitive_structure, max_index, num_type, species, list_num_exp):
        for index, num_exp in zip(range(1, max_index + 1), list_num_exp):
            list_ds = enumerate_derivative_structures(primitive_structure, index, num_type)
            self.assertEqual(len(list_ds), num_exp)

            dstructs = [ds.get_structure(species) for ds in list_ds]
            stm = StructureMatcher()
            grouped = stm.group_structures(dstructs)
            self.assertEqual(len(grouped), num_exp)


if __name__ == '__main__':
    unittest.main()
