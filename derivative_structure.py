import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.analysis.structure_matcher import StructureMatcher
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


class Superlattice(object):
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    lattice_vectors: array, (dim, dim)
        lattice_vectors[:, i] is the i-th lattice vector

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    list_species: list of DummySpecie
    snf: array, (dim, dim)
        Smith normal form of hnf
    left: array, (dim, dim)
        left unimodular matrix
    left_inv: array, (dim, dim)
        inverse of left matrix
    right: array, (dim, dim)
        right unimodular matrix
    """

    def __init__(self, hnf, lattice_vectors):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())
        self.lattice_vectors = lattice_vectors
        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)

        self.struct = self._get_structure()

    def _get_structure(self):
        # (dim, index)
        lattice_factors_e = np.array([
            np.unravel_index(indices, tuple(self.snf.diagonal()))
            for indices in range(self.index)]).T
        # (dim, index)
        self.lattice_points = np.dot(self.lattice_vectors,
                                     np.dot(self.left_inv,
                                            lattice_factors_e))
        frac_coords = np.linalg.solve(self.hnf,
                                      np.dot(self.left_inv,
                                             lattice_factors_e)).T
        frac_coords = np.mod(frac_coords, np.ones(self.dim))

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)
        list_species = [DummySpecie('X')] * self.index

        struct = Structure(lattice, list_species, frac_coords)
        return struct


def unique_structures(structures):
    uniqued = []
    stm = StructureMatcher()

    for struct in structures:
        sp = SpacegroupAnalyzer(struct)
        sym_dataset = sp.get_symmetry_dataset()

        is_unique = True
        for r, t in zip(sym_dataset['rotations'], sym_dataset['translations']):
            new_frac_coords = np.dot(r, struct.frac_coords.T) + t[:, np.newaxis]
            struct_tmp = Structure(struct.lattice, struct.species, new_frac_coords.T)
            if any([stm.fit(struct_tmp, unique_struct) for unique_struct in uniqued]):
                is_unique = False
                break

        if is_unique:
            uniqued.append(struct)

    return uniqued


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
    parent_lattice = Structure(structure.lattice, [DummySpecie('X')],
                               coords=[[0, 0, 0]])
    sym_dataset = SpacegroupAnalyzer(parent_lattice).get_symmetry_dataset()
    rotations = sym_dataset['rotations']
    translations = sym_dataset['translations']
    return rotations, translations


if __name__ == '__main__':
    hcp = get_lattice('hcp')
    print(hcp)
