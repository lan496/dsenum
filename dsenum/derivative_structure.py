import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite

from dsenum.converter import DerivativeMultiLatticeHash, get_species_list


class ColoringToStructure:
    def __init__(
        self,
        base_structure: Structure,
        dshash: DerivativeMultiLatticeHash,
        mapping_color_to_species,
    ):
        self.base_structure = base_structure
        self.dshash = dshash
        self.mapping_color_to_species = mapping_color_to_species

        self.lattice_matrix = np.dot(self.base_matrix.T, self.dshash.hnf).T
        # lattice of derivative structure
        self.lattice = Lattice(self.lattice_matrix)

        self.canonical_derivative_sites = dshash.get_canonical_and_derivative_sites_list()

        list_coords = []
        for _, dsite in self.canonical_derivative_sites:
            coords = dshash.get_frac_coords(dsite)
            cart_coords = np.dot(coords, self.base_matrix)
            list_coords.append(cart_coords)
        self.list_coords = list_coords

        self.precomputed_psites = [
            [
                PeriodicSite(sp, coords, self.lattice, coords_are_cartesian=True)
                for sp in self.mapping_color_to_species
            ]
            for coords in self.list_coords
        ]

    @property
    def base_matrix(self):
        return self.base_structure.lattice.matrix

    def convert_to_structure(self, coloring) -> Structure:
        list_psites = [
            self.precomputed_psites[i][coloring[i]] for i in range(self.dshash.num_sites)
        ]
        dstruct = Structure.from_sites(list_psites)

        return dstruct


def get_supercell(structure: Structure, scaling_matrix: np.ndarray):
    """
    Parameters
    ----------
    structure: pymatgen.core.Structure
    scaling_matrix: array, lower triangular integer matrix
        the i-th new lattice vector is np.dot(scaling_matrix.T, structure.lattice.matrix)[i, :]

    Returns
    -------
    ds: pymatgen.core.Structure
    """
    dshash = DerivativeMultiLatticeHash(scaling_matrix, structure.frac_coords)

    base_lattice_matrix = structure.lattice.matrix
    lattice = Lattice(np.dot(scaling_matrix.T, base_lattice_matrix))

    species = get_species_list(dshash.index, [site.species for site in structure])

    list_dsites = dshash.get_distinct_derivative_sites_list()
    new_cart_coords = [
        np.dot(dshash.get_frac_coords(dsite), base_lattice_matrix) for dsite in list_dsites
    ]

    ds = Structure(lattice, species, new_cart_coords, coords_are_cartesian=True)
    return ds
