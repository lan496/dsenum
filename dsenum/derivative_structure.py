import numpy as np
from pymatgen.core import Lattice
from pymatgen.core.structure import Structure
from pymatgen.core.sites import PeriodicSite
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


def coloring_to_derivative_structure(base_structure: Structure, dshash: DerivativeMultiLatticeHash,
                                     mapping_color_to_species, coloring) -> Structure:
    sites = []
    for dsite in dshash.get_distinct_derivative_sites_list():
        site_index, dimage = dsite

        csite = dshash.hash_derivative_site(dsite)
        indices = dshash.hash_canonical_site(csite)
        species = mapping_color_to_species[coloring[indices]]

        coords = dshash.get_frac_coords(dsite)
        base_matrix = base_structure.lattice.matrix
        cart_coords = np.dot(coords, base_matrix)

        lattice_matrix = np.dot(base_matrix.T, dshash.hnf).T
        lattice = Lattice(lattice_matrix)
        psite = PeriodicSite(species, cart_coords, lattice, coords_are_cartesian=True)
        sites.append(psite)

    dstruct = Structure.from_sites(sites)
    return dstruct
