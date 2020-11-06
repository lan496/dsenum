import os

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import Specie, DummySpecie
from tqdm import tqdm

from dsenum import ZddStructureEnumerator, StructureEnumerator
from dsenum.utils import write_cif


def get_fcc_structure():
    a = 2.856

    lattice = Lattice.from_parameters(a, a, a, 60, 60, 60)
    species = ["Cu"]
    frac_coords = np.array([[0, 0, 0]])
    structure = Structure(lattice, species, frac_coords)
    return structure


if __name__ == "__main__":
    # enumerate fcc derivative structures
    aristo = get_fcc_structure()

    max_index = 4

    mapping_color_species = [Specie("Cu"), Specie("Au")]
    num_types = len(mapping_color_species)
    composition_constraints = None
    base_site_constraints = None

    dirname = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fcc")
    os.makedirs(dirname, exist_ok=True)

    for index in range(2, max_index + 1):
        print(f"index={index}")
        zse = ZddStructureEnumerator(
            aristo,
            index,
            num_types,
            mapping_color_species=mapping_color_species,
            composition_constraints=composition_constraints,
            base_site_constraints=base_site_constraints,
            remove_superperiodic=True,
            remove_incomplete=True,
        )
        list_dstructs = zse.generate()

        for i, dstruct in enumerate(tqdm(list_dstructs)):
            filename = os.path.join(dirname, f"fcc_{index}_{i}.cif")
            write_cif(filename, dstruct, refine_cell=True)
