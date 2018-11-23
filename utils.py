import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie


def get_lattice(kind):
    if kind == 'hcp':
        latt = Lattice(np.array([[1, 0, 0],
                                 [0.5, np.sqrt(3) / 2, 0],
                                 [0, 0, 2 * np.sqrt(6) / 3]]))
        coords = [[0, 0, 0], [1 / 3, 1 / 3, 0.5]]
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


def get_fcc_with_vacancy():
    latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
    displacement_set = [[0, 0, 0],  # lattice point
                        [0.25, 0.25, 0.25],  # tetrahedral site
                        [0.5, 0.5, 0.5],  # octahedral site
                        [0.75, 0.75, 0.75],  # tetrahedral site
                        ]
    struct = Structure(latt, [DummySpecie('X')] * len(displacement_set), displacement_set)
    return struct


def get_symmetry_operations(structure):
    """
    Parameters
    ----------
    structure: pymatgen.core.Structure

    Returns
    -------
    rotations: array, (# of symmetry operations, 3, 3)
    translations: array, (# of symmetry operations, 3)
    """
    sym_dataset = SpacegroupAnalyzer(structure).get_symmetry_dataset()
    rotations = sym_dataset['rotations']
    translations = sym_dataset['translations']
    return rotations, translations


def unique_structures(structures):
    uniqued = []
    stm = StructureMatcher()

    for struct in structures:
        rotations, translations = get_symmetry_operations(struct)

        is_unique = True
        for r, t in zip(rotations, translations):
            new_frac_coords = np.dot(r, struct.frac_coords.T) + t[:, np.newaxis]
            struct_tmp = Structure(struct.lattice, struct.species, new_frac_coords.T)
            if any([stm.fit(struct_tmp, unique_struct) for unique_struct in uniqued]):
                is_unique = False
                break

        if is_unique:
            uniqued.append(struct)

    return uniqued
