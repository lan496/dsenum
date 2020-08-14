from typing import List, Union

import numpy as np
from pymatgen.core import Lattice, Structure, Specie, DummySpecie
from pymatgen.core.sites import PeriodicSite

from dsenum.converter import DerivativeMultiLatticeHash, get_species_list


class ColoringToStructure:
    def __init__(
        self,
        base_structure: Structure,
        dshash: DerivativeMultiLatticeHash,
        mapping_color_to_species: List[Union[Specie, DummySpecie]],
        additional_species=None,
        additional_frac_coords=None,
    ):
        """
        Parameters
        ----------
        base_structure:
            Structure with only ordering species
        dshash
        mapping_color_to_species
        additional_species: list of pymatgen.core.Species, optional
            species which are nothing to do with ordering
        additional_frac_coords: np.ndarray, optional
            fractional coordinates of species which are nothing to do with ordering
        """
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

        # additional fixed sites
        self.additional_species = additional_species
        self.additional_frac_coords = additional_frac_coords
        if self.additional_species is not None:
            assert len(self.additional_species) == len(self.additional_frac_coords)

        self.additional_psites = []
        if self.additional_species is not None:
            lattice_points = self.dshash.get_lattice_points()
            for sp, disp in zip(self.additional_species, self.additional_frac_coords):
                cart_coords = [
                    np.dot(np.array(disp) + np.array(lp), self.base_matrix)
                    for lp in lattice_points
                ]
                self.additional_psites.extend(
                    [
                        PeriodicSite(sp, coords, self.lattice, coords_are_cartesian=True)
                        for coords in cart_coords
                    ]
                )

    @property
    def base_matrix(self):
        return self.base_structure.lattice.matrix

    def convert_to_structure(self, coloring) -> Structure:
        list_psites = [
            self.precomputed_psites[i][coloring[i]] for i in range(self.dshash.num_sites)
        ]
        if self.additional_psites is not None:
            list_psites.extend(self.additional_psites)

        dstruct = Structure.from_sites(list_psites)

        return dstruct
