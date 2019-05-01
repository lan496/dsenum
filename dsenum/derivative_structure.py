import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
# from pymatgen.core.sites import PeriodicSite
from pymatgen.core.periodic_table import DummySpecie

from dsenum.permutation import DerivativeStructureHash
from dsenum.permutation_group import DerivativeMultiLatticeHash


class DerivativeStructure:
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    num_type: int
        # of kinds of atoms
    lattice_vectors: array, (dim, dim)
        lattice_vectors[:, i] is the i-th lattice vector
    labeling: list
    num_site_parent: int
        # of atoms in parent multilattice
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates of A in multilattice site
    """

    def __init__(self, hnf, num_type, lattice_vectors, labeling,
                 num_site_parent=1, displacement_set=None):
        self.hnf = hnf
        self.num_site_parent = num_site_parent
        self.dshash = DerivativeStructureHash(hnf, num_site_parent)

        self.lattice_vectors = lattice_vectors
        self.num_type = num_type
        self.labeling = labeling

        if self.num_site_parent == 1:
            self.displacement_set = np.array([[0, 0, 0]])
        else:
            self.displacement_set = displacement_set
            assert self.displacement_set.shape[0] == self.num_site_parent

    @property
    def left_inv(self):
        return self.dshash.left_inv

    def get_structure(self, species=None):
        if species is None:
            self.species = [DummySpecie(str(i)) for i in range(1, self.num_type + 1)]
        else:
            self.species = species
        self.list_species = [self.species[idx] for idx in self.labeling]

        # (1 + dim, num_site)
        factors_e = self.dshash.get_distinct_factors().T
        # (dim, num_site)
        parent_frac_coords = self.displacement_set[factors_e[0, :]].T \
            + np.dot(self.left_inv, factors_e[1:, :])

        frac_coords = np.linalg.solve(self.hnf, parent_frac_coords).T

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)

        struct = Structure(lattice, self.list_species, frac_coords)
        return struct


class ColoringToStructure:

    def __init__(self, base_structure: Structure,
                 dshash: DerivativeMultiLatticeHash,
                 mapping_color_to_species):
        self.base_structure = base_structure
        self.dshash = dshash
        self.mapping_color_to_species = mapping_color_to_species

        self.lattice_matrix = np.dot(self.base_matrix.T, self.dshash.hnf).T
        self.lattice = Lattice(self.lattice_matrix)

        self.canonical_derivative_sites = dshash.get_canonical_and_derivative_sites_list()
        self.csite_indices = [self.dshash.hash_canonical_site(csite)
                              for csite, _ in self.canonical_derivative_sites]

        list_coords = []
        for _, dsite in self.canonical_derivative_sites:
            coords = dshash.get_frac_coords(dsite)
            cart_coords = np.dot(coords, self.base_matrix)
            list_coords.append(cart_coords)

        self.list_coords = list_coords

    @property
    def base_matrix(self):
        return self.base_structure.lattice.matrix

    def convert_to_structure(self, coloring) -> Structure:
        list_species = []
        for (csite, dsite), indices in zip(self.canonical_derivative_sites, self.csite_indices):
            site_index, dimage = dsite
            species = self.mapping_color_to_species[coloring[indices]]
            list_species.append(species)

        dstruct = Structure(self.lattice, list_species, self.list_coords,
                            coords_are_cartesian=True)

        return dstruct
