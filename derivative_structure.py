import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie

from derivative.smith_normal_form import smith_normal_form


class SuperMultilattice:
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    lattice_vectors: array, (dim, dim)
        lattice_vectors[:, i] is the i-th lattice vector
    num_site_parent: int
        # of atoms in parent multilattice
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates of A in multilattice site

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    num_site: int
        # of sites in unit cell of superlattice
    shape: tuple of int
        (1 + dim, num_site)
    snf: array, (dim, dim)
        Smith normal form of hnf
    left: array, (dim, dim)
        left unimodular matrix
    left_inv: array, (dim, dim)
        inverse of left matrix
    right: array, (dim, dim)
        right unimodular matrix
    struct: pymatgen.core.Structure
    """

    def __init__(self, hnf, lattice_vectors,
                 num_site_parent=1, displacement_set=None):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())
        self.lattice_vectors = lattice_vectors

        self.num_site_parent = num_site_parent
        if self.num_site_parent == 1:
            self.displacement_set = np.array([[0, 0, 0]])
        else:
            self.displacement_set = displacement_set
            assert self.displacement_set.shape[0] == self.num_site_parent
        self.num_site = self.num_site_parent * self.index

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)

        self.shape = tuple([self.num_site_parent]
                           + self.snf.diagonal().tolist())

        self.list_species = [DummySpecie('X')] * self.num_site

    def get_structure(self):
        # (1 + dim, num_site)
        factors_e = np.array([np.unravel_index(indices, self.shape)
                              for indices in range(self.num_site)]).T

        # (dim, num_site)
        parent_frac_coords = self.displacement_set[factors_e[0, :]].T \
            + np.dot(self.left_inv, factors_e[1:, :])

        frac_coords = np.linalg.solve(self.hnf, parent_frac_coords).T

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)

        struct = Structure(lattice, self.list_species, frac_coords)
        return struct


class DerivativeStructure(SuperMultilattice):

    def __init__(self, hnf, num_type, lattice_vectors, labeling,
                 num_site_parent=1, displacement_set=None):
        super().__init__(hnf, lattice_vectors, num_site_parent, displacement_set)
        self.num_type = num_type
        self.labeling = labeling

    def get_structure(self, species=None):
        if species is None:
            self.species = [DummySpecie(str(i)) for i in range(1, self.num_type + 1)]
        else:
            self.species = species
        self.list_species = [self.species[idx] for idx in self.labeling]
        return super().get_structure()

    def draw(self, ax=None):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        for indices in itertools.product(*[range(e + 1) for e in self.hnf.diagonal().tolist()]):
            src = np.dot(self.lattice_vectors, np.array(indices))
            for i in range(3):
                directions = np.eye(3)
                dst = src + directions[:, i]
                tmp = np.concatenate([src, dst]).reshape(2, 3).T
                ax.plot(tmp[0], tmp[1], tmp[2])

        superlattice_vectors = np.dot(self.lattice_vectors, self.hnf)
        for i in range(3):
            origin = [0, 0, 0]
            ax.quiver(*origin, *superlattice_vectors[:, i].tolist(), arrow_length_ratio=0)
