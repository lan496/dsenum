from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core.periodic_table import Specie

from dsenum.enumerate import enumerate_derivative_structures
from dsenum.utils import get_lattice, write_cif  # noqa


def test_fcc():
    base_structure = get_lattice("fcc")
    num_type = 2
    indices = [2, 3, 4]
    species = [Specie("Cu"), Specie("Au")]

    check(base_structure, num_type, indices, species, "fcc")


def test_sc():
    base_structure = get_lattice("sc")
    num_type = 2
    indices = [2, 3, 4]
    species = [Specie("Cu"), Specie("Au")]

    check(base_structure, num_type, indices, species, "sc")


def test_hcp():
    base_structure = get_lattice("hcp")
    num_type = 2
    indices = [2, 3, 4]
    species = [Specie("Cu"), Specie("Au")]
    # inconsistent with StructureMatcher in index=4

    check(base_structure, num_type, indices, species, "hcp")


def check(base_structure, num_type, indices, species, name):
    for index in indices:
        list_ds = enumerate_derivative_structures(
            base_structure,
            index,
            num_type,
            species,
            color_exchange=True,
            leave_superperiodic=False,
        )

        stm = StructureMatcher(ltol=1e-4, stol=1e-4)
        grouped = stm.group_structures(list_ds)
        assert len(grouped) == len(list_ds)
