import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure, Specie, DummySpecie
import pytest

from dsenum import StructureEnumerator
from dsenum.utils import get_lattice


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
        se = StructureEnumerator(
            base_structure,
            index,
            num_type,
            species,
            color_exchange=True,
            leave_superperiodic=False,
        )
        list_ds = se.generate()

        stm = StructureMatcher(ltol=1e-4, stol=1e-4)
        grouped = stm.group_structures(list_ds)
        assert len(grouped) == len(list_ds)


def test_coloring_with_fixed_species():
    lattice = Lattice(3.945 * np.eye(3))
    species = ["Sr", "Ti", "O", "O", "O"]
    frac_coords = np.array(
        [[0, 0, 0], [0.5, 0.5, 0.5], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]
    )
    aristo = Structure(lattice, species, frac_coords)
    base_structure = aristo.copy()
    base_structure.remove_species(["Sr", "Ti"])
    additional_species = species[:2]
    additional_frac_coords = frac_coords[:2]

    mapping_color_species = [DummySpecie("X"), "O"]
    num_types = len(mapping_color_species)
    index = 2

    se = StructureEnumerator(
        base_structure,
        index,
        num_types,
        mapping_color_species=mapping_color_species,
        color_exchange=False,
        leave_superperiodic=False,
        use_all_colors=False,
    )
    list_dstructs = se.generate(
        additional_species=additional_species, additional_frac_coords=additional_frac_coords
    )

    # with base_site_constraints
    mapping_color_species2 = [DummySpecie("X"), "O", "Sr", "Ti"]
    num_types2 = len(mapping_color_species2)
    base_site_constraints = [
        [2],  # Sr site
        [3],  # Cu site
        [0, 1],  # O or V
        [0, 1],  # O or V
        [0, 1],  # O or V
    ]
    se2 = StructureEnumerator(
        aristo,
        index,
        num_types2,
        mapping_color_species=mapping_color_species2,
        base_site_constraints=base_site_constraints,
        color_exchange=False,
        leave_superperiodic=False,
        use_all_colors=False,
    )
    list_dstructs2 = se2.generate()

    # check uniqueness by StructureMatcher
    stm = StructureMatcher(ltol=1e-4, stol=1e-4)
    grouped = stm.group_structures(list_dstructs + list_dstructs2)
    assert len(grouped) == len(list_dstructs)
    assert all([(len(matched) == 2) for matched in grouped])
