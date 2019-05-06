import os
import unittest

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Specie

from dsenum.enumerate import enumerate_derivatives
from dsenum.utils import get_lattice, write_cif


class TestDerivativeStructure(unittest.TestCase):

    def test_fcc(self):
        base_structure = get_lattice('fcc')
        num_type = 2
        indices = [2, 3, 4]
        species = [Specie('Cu'), Specie('Au')]

        self.check(base_structure, num_type, indices, species, 'fcc')

    def test_sc(self):
        base_structure = get_lattice('sc')
        num_type = 2
        indices = [2, 3, 4]
        species = [Specie('Cu'), Specie('Au')]

        self.check(base_structure, num_type, indices, species, 'sc')

    def test_hcp(self):
        base_structure = get_lattice('hcp')
        num_type = 2
        indices = [2, 3, 4]
        species = [Specie('Cu'), Specie('Au')]
        # inconsistent with StructureMatcher in index=4

        self.check(base_structure, num_type, indices, species, 'hcp')

    def check(self, base_structure, num_type, indices, species, name):
        os.makedirs(os.path.join('tests', name), exist_ok=True)
        for index in indices:
            list_ds = enumerate_derivatives(base_structure, index, num_type, species,
                                            color_exchange=True,
                                            leave_superperiodic=False)

            stm = StructureMatcher(ltol=1e-4, stol=1e-4)
            grouped = stm.group_structures(list_ds)
            self.assertEqual(len(grouped), len(list_ds))

            """
            for i, dstruct in enumerate(list_ds):
                filename = os.path.join('tests', name, '{}_N={}_{}.cif'.format(name, index, i))
                write_cif(filename, dstruct, refine_cell=True, resize_volume=True)
            """


if __name__ == '__main__':
    unittest.main()
