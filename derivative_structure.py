import itertools

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer

from smith_normal_form import smith_normal_form
from permutation import is_same_lattice


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


class SuperLattice(object):

    def __init__(self, hnf, lattice_vectors):
        self.hnf = hnf
        self.num_sites = np.prod(self.hnf.diagonal())
        self.lattice_vectors = lattice_vectors
        # self.list_species = [DummySpecie(str(i))
        #                      for i in range(1, self.num_sites + 1)]
        self.list_species = [DummySpecie('X')
                             for _ in range(self.num_sites)]

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)

        self.struct = self._get_structure()

    def _get_structure(self):
        # (dims, num)
        factors_e = np.array([
            np.unravel_index(indices, tuple(self.snf.diagonal()))
            for indices in range(self.num_sites)]).T
        # (dims, num)
        self.points = np.dot(self.lattice_vectors,
                             np.dot(self.left_inv, factors_e))
        frac_coords = np.linalg.solve(
            np.dot(self.lattice_vectors, self.hnf), self.points).T
        self.frac_coords = np.mod(frac_coords, np.ones(self.hnf.shape[0]))

        lattice = Lattice(np.dot(self.lattice_vectors, self.hnf).T)
        struct = Structure(lattice, self.list_species, self.frac_coords)
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


def check_valid_rotations(dstruct, rotations):
    rotations_sl = [r.astype(np.int) for r in rotations
                    if is_same_lattice(np.dot(r, dstruct.hnf), dstruct.hnf)]
    rotations_sl = dstruct.labeling.permutation.rotations
    rotations_spglib = SpacegroupAnalyzer(dstruct.struct).get_symmetry_dataset()['rotations']
    print(len(rotations_sl), len(rotations_spglib))
    if len(rotations_sl) == len(rotations_spglib):
        return True
    else:
        print('implementation')
        print(rotations_sl)
        print('spglib')
        print(rotations_spglib)
        return False
