import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer

from smith_normal_form import smith_normal_form


class DerivativeStructure(object):

    def __init__(self, hnf, num_type, lattice_vectors, labeling):
        self.hnf = hnf
        self.num_sites = np.prod(self.hnf.diagonal())
        self.num_type = num_type
        self.lattice_vectors = lattice_vectors
        self.labeling = labeling

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)

        self.list_species = [DummySpecie(str(i)) for i in range(1, self.num_type + 1)]

        self.struct = self._get_structure()

    def _get_structure(self):
        # (dims, num)
        factors_e = np.array([
            np.unravel_index(indices, tuple(self.snf.diagonal()))
            for indices in range(self.num_sites)]).T
        # (dims, num)
        self.points = np.dot(self.lattice_vectors, np.dot(self.left_inv, factors_e))
        frac_coords = np.linalg.solve(np.dot(self.lattice_vectors, self.hnf), self.points).T
        self.frac_coords = np.mod(frac_coords, np.ones(self.hnf.shape[0]))

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)
        species = [self.list_species[idx] for idx in self.labeling]
        struct = Structure(lattice, species, self.frac_coords)
        return struct

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


class SuperMultilattice(object):
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

        self.struct = self._get_structure()

    def _get_structure(self):
        # (1 + dim, num_site)
        factors_e = np.array([np.unravel_index(indices, self.shape)
                              for indices in range(self.num_site)]).T

        # (dim, num_site)
        parent_frac_coords = self.displacement_set[factors_e[0, :]].T \
            + np.dot(self.left_inv, factors_e[1:, :])

        frac_coords = np.linalg.solve(self.hnf, parent_frac_coords).T

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)
        list_species = [DummySpecie('X')] * self.num_site

        struct = Structure(lattice, list_species, frac_coords)
        return struct


def get_lattice(kind):
    if kind == 'hcp':
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
        coords = [[0, 0, 0], [0.5, np.sqrt(3) / 2, np.sqrt(6) / 3]]
    else:
        coords = [[0, 0, 0]]

    if kind == 'sc':
        latt = Lattice(np.eye(3))
    elif kind == 'fcc':
        latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
    elif kind == 'bcc':
        latt = Lattice(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]))
    elif kind == 'hex':
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
    elif kind == 'tet':
        latt = Lattice(np.diag([1, 1, 1.2]))

    struct = Structure(latt, [DummySpecie('X')] * len(coords), coords)
    return struct


def get_symmetry_operations(structure, parent_lattice=False):
    """
    Parameters
    ----------
    structure: pymatgen.core.Structure
    parent_lattice: bool
        if True, return symmetry operations of parent lattice of a given structure

    Returns
    -------
    rotations: array, (# of symmetry operations, 3, 3)
    translations: array, (# of symmetry operations, 3)
    """
    if parent_lattice:
        struct = Structure(structure.lattice, [DummySpecie('X')],
                           coords=[[0, 0, 0]])
    else:
        struct = structure

    sym_dataset = SpacegroupAnalyzer(struct).get_symmetry_dataset()
    rotations = sym_dataset['rotations']
    translations = sym_dataset['translations']
    return rotations, translations
