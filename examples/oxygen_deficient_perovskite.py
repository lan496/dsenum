import os

import numpy as np
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie, Specie

from dsenum import StructureEnumerator
from dsenum.utils import refine_and_resize_structure, write_cif

if __name__ == "__main__":
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
        remove_superperiodic=True,
        remove_incomplete=False,
    )
    list_dstructs = se.generate(
        additional_species=additional_species, additional_frac_coords=additional_frac_coords
    )  # type: ignore

    # leave only SrTiO_{3-x} (0 <= x <= 1)
    num_oxygen_lb = index * 2
    num_oxygen_ub = index * 3
    structures = []
    for ds in list_dstructs:
        num_oxygen = sum((str(sp) == "O") for sp in ds.species)  # type: ignore
        if num_oxygen_lb <= num_oxygen and num_oxygen <= num_oxygen_ub:
            structures.append(ds)
    print(index, len(structures))

    dirname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SrTiO3-x")
    os.makedirs(dirname, exist_ok=True)
    for i, dstruct in enumerate(structures):
        # remove void
        dstruct.remove_species([mapping_color_species[0]])  # type: ignore

        filename = os.path.join(dirname, f"{index}_{i}.cif")
        write_cif(filename, dstruct, refine_cell=True)
