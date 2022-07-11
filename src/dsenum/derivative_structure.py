from typing import List, Union

import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie, Element, Specie
from pymatgen.core.sites import PeriodicSite

import dsenum
from dsenum.converter import DerivativeMultiLatticeHash


class ColoringToStructure:
    def __init__(
        self,
        base_structure: Structure,
        dshash: DerivativeMultiLatticeHash,
        mapping_color_to_species: List[Union[str, Element, Specie, DummySpecie]],
        additional_species=None,
        additional_frac_coords=None,
    ):
        """
        Parameters
        ----------
        base_structure:
            Structure with only ordering species
        dshash:
        mapping_color_to_species:
        additional_species: list of pymatgen.core.Species, optional
            species which are nothing to do with ordering
        additional_frac_coords: np.ndarray, optional
            fractional coordinates of species which are nothing to do with ordering
        """
        self.base_structure = base_structure
        self.dshash = dshash
        self.mapping_color_to_species = mapping_color_to_species
        self.num_sites = self.dshash.num_sites

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

        # map color to specie str (e.g. 0 -> "Cu")
        self.mapping_color_to_species_str = [
            str(specie) for specie in self.mapping_color_to_species
        ]

        self._precompute_for_poscar_string()

    def _precompute_for_poscar_string(self):
        # precompute comment and lattice vectors for POSCAR output
        self.head_lines = []
        version = dsenum.__version__
        comment = "generated by dsenum " + version
        self.head_lines.append(comment)

        # scale
        self.head_lines.append("1.0")

        # lattice vectors
        self.head_lines.extend(
            [
                str(self.lattice_matrix[i, 0])
                + " "
                + str(self.lattice_matrix[i, 1])
                + " "
                + str(self.lattice_matrix[i, 2])
                for i in range(3)
            ]
        )

        self.precompute_coords_str = [
            [
                str(psite.frac_coords[0])
                + " "
                + str(psite.frac_coords[1])
                + " "
                + str(psite.frac_coords[2])
                for psite in list_psites
            ]
            for list_psites in self.precomputed_psites
        ]

        if self.additional_species is not None:
            assert len(self.additional_species) == len(self.additional_frac_coords)
            self.additional_species_str = [str(specie) for specie in self.additional_species]
        else:
            self.additional_species_str = None

        if self.additional_psites is not None:
            self.additional_coords_str = [
                str(psite.frac_coords[0])
                + " "
                + str(psite.frac_coords[1])
                + " "
                + str(psite.frac_coords[2])
                for psite in self.additional_psites
            ]
        else:
            self.additional_coords_str = None

    @property
    def base_matrix(self):
        return self.base_structure.lattice.matrix

    def convert_to_structure(self, coloring) -> Structure:
        list_psites = [self.precomputed_psites[i][coloring[i]] for i in range(self.num_sites)]
        if self.additional_psites is not None:
            list_psites.extend(self.additional_psites)

        dstruct = Structure.from_sites(list_psites)

        return dstruct

    def convert_to_poscar_string(self, coloring) -> str:
        list_coords_str = [
            self.precompute_coords_str[i][coloring[i]] for i in range(self.num_sites)
        ]
        if self.additional_psites is not None:
            list_coords_str.extend([s for s in self.additional_coords_str])

        list_species = [
            self.mapping_color_to_species_str[coloring[i]] for i in range(self.num_sites)
        ]

        head_lines = self.head_lines[::]
        poscar_str = prepare_poscar_string(head_lines, list_species, list_coords_str)
        return poscar_str


def prepare_poscar_string(
    head_lines: List[str], list_species: List[str], list_coords_str: List[str]
) -> str:
    # ref: https://www.vasp.at/wiki/index.php/POSCAR
    lines = head_lines

    # species
    # TODO: this part is dominant in runtime
    counter = []
    element = ""
    count = 0
    for species_str in list_species:
        if element == "":
            element = species_str
            count += 1
        elif species_str == element:
            count += 1
        else:
            counter.append((element, count))
            element = species_str
            count = 1
    if element != "":
        counter.append((element, count))
    lines.append(" ".join([e for e, _ in counter]))
    lines.append(" ".join([str(c) for _, c in counter]))

    # fractional coords
    lines.append("Direct")
    lines.extend(list_coords_str)

    poscar_str = "\n".join(lines)
    return poscar_str