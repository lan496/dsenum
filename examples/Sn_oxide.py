import os

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie, Specie

from dsenum import StructureEnumerator
from dsenum.utils import write_cif


def get_rutile_structure():
    # rutile structure taken from mp-856
    a = 4.832
    c = 3.243
    x_4f = 0.3066

    lattice = Lattice.from_parameters(a, a, c, 90, 90, 90)
    species = ["Sn", "Sn", "O", "O", "O", "O"]
    # fmt: off
    frac_coords = np.array([
        [0, 0, 0],                      # Sn(2a)
        [0.5, 0.5, 0.5],                # Sn(2a)
        [x_4f, x_4f, 0],                # O(4f)
        [1 - x_4f, 1 - x_4f, 0],        # O(4f)
        [0.5 - x_4f, 0.5 + x_4f, 0.5],  # O(4f)
        [0.5 + x_4f, 0.5 - x_4f, 0.5],  # O(4f)
    ])
    # fmt: on
    structure = Structure(lattice, species, frac_coords)
    return structure


if __name__ == "__main__":
    # enumerate rutile-like SnO1-x derivative structures
    rutile = get_rutile_structure()

    max_index = 3

    mapping_color_species = [DummySpecie("X"), Specie("O"), Specie("Sn")]
    num_types = len(mapping_color_species)
    composition_constraints = None
    # fmt: off
    base_site_constraints = [
        [2],  # 2a
        [2],  # 2a
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
    ]
    # fmt: on

    dirname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "SnO1-x")
    os.makedirs(dirname, exist_ok=True)

    for index in range(1, max_index + 1):
        se = StructureEnumerator(
            rutile,
            index,
            num_types,
            mapping_color_species=mapping_color_species,
            composition_constraints=composition_constraints,
            base_site_constraints=base_site_constraints,
            color_exchange=False,
            remove_superperiodic=True,
            remove_incomplete=False,
        )
        list_dstructs = se.generate()  # type: ignore

        for i, dstruct in enumerate(list_dstructs):
            # remove void
            dstruct.remove_species([mapping_color_species[0]])  # type: ignore

            filename = os.path.join(dirname, f"SnO1-x_{index}_{i}.cif")
            write_cif(filename, dstruct, refine_cell=True)
