import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor


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
    find symmetry operations for a given structure

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


def write_cif(filename, struct, refine_cell=False, resize_volume=False):
    struct = refine_and_resize_structure(struct, refine_cell, resize_volume)
    if not struct.is_valid(1e-4):
        return
    cw = CifWriter(struct)
    cw.write_file(filename)


def refine_and_resize_structure(struct, refine_cell=True, resize_volume=True):
    if resize_volume:
        dls = DLSVolumePredictor()
        struct = dls.get_predicted_structure(struct)
        struct.apply_strain(0.5)

    if refine_cell:
        sga = SpacegroupAnalyzer(struct, symprec=1e-6, angle_tolerance=1e-2)
        struct = sga.get_primitive_standard_structure()

    return struct
