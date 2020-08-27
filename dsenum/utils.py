import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor

from ._version import get_versions  # type: ignore


def get_lattice(kind):
    if kind == "hcp":
        latt = Lattice(np.array([[1, 0, 0], [0.5, np.sqrt(3) / 2, 0], [0, 0, 2 * np.sqrt(6) / 3]]))
        coords = [[0, 0, 0], [1 / 3, 1 / 3, 0.5]]
    else:
        coords = [[0, 0, 0]]

    if kind == "sc":
        latt = Lattice(np.eye(3))
    elif kind == "fcc":
        latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
    elif kind == "bcc":
        latt = Lattice(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]))
    elif kind == "hex":
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
    elif kind == "tet":
        latt = Lattice(np.diag([1, 1, 1.2]))

    struct = Structure(latt, [DummySpecie("X")] * len(coords), coords)
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
    rotations = sym_dataset["rotations"]
    translations = sym_dataset["translations"]
    return rotations, translations


def write_cif(
    filename: str,
    struct: Structure,
    refine_cell=False,
    resize_volume=False,
    symprec=1e-2,
    comment="",
):
    """
    dump structure in CIF format after resizing to feasible volume and refing cell by symmetry.

    Parameters
    ----------
    filename: str
    struct: pymatgen.core.Structure
    refine_cell: bool, optional
        if true, refine cell setting by spglib
    resize_volume: bool, optional
        if true, resize lattice by DLSVolumePredictor in pymatgen
    symprec: float, optional
        symprec in spglib
    """
    struct = refine_and_resize_structure(struct, refine_cell, resize_volume)
    cw = CifWriter(struct, symprec=symprec)
    version = get_versions()["version"]
    comment = f"# generated by {version}\n" + comment
    with open(filename, "w") as f:
        if comment:
            f.write(comment + "\n" + cw.__str__())
        else:
            f.write(cw.__str__())


def refine_and_resize_structure(struct, refine_cell=True, resize_volume=True):
    if resize_volume:
        dls = DLSVolumePredictor()
        struct = dls.get_predicted_structure(struct)
        struct.apply_strain(0.5)

    if refine_cell:
        sga = SpacegroupAnalyzer(struct)
        # SpacegroupAnalyzer.get_primitive_standard_structure may return incorrect structures
        # So we avoid to use SpacegroupAnalyzer.get_primitive_standard_structure
        struct = sga.get_refined_structure()

    return struct


def cast_integer_matrix(arr: np.ndarray) -> np.ndarray:
    arr_int = np.around(arr).astype(np.int)
    return arr_int
