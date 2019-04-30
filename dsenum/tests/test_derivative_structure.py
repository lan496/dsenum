import os
import unittest

from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Specie
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from dsenum.enumerate import enumerate_derivatives
from dsenum.utils import get_lattice


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


def write_cif(filename, struct, refine_cell=False, resize_volume=False):
    struct = refine_and_resize_structure(struct, refine_cell, resize_volume)
    if not struct.is_valid(1e-4):
        return
    cw = CifWriter(struct)
    cw.write_file(filename)


def refine_and_resize_structure(struct, refine_cell=True, resize_volume=True):
    if resize_volume:
        dls = DLSVolumePredictor()
        struct = dls.get_predicted_structure(struct)
        struct.apply_strain(0.5)

    if refine_cell:
        sga = SpacegroupAnalyzer(struct, symprec=1e-6, angle_tolerance=1e-2)
        struct = sga.get_primitive_standard_structure()

    return struct


if __name__ == '__main__':
    unittest.main()
