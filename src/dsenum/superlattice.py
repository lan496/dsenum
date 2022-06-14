from typing import List

import numpy as np
from pymatgen.core import Structure

from dsenum.utils import get_symmetry_operations


def generate_all_superlattices(index: int) -> List[np.ndarray]:
    """
    enumerate 3-by-3 Hermite normal form which determinant is equal to `index`

    Parameters
    ----------
    index: positive integer

    Returns
    -------
    list_HNF: list of matrices, each element is Hermite normal form
    """

    def make_HNF(a, b, c, d, e, f):
        arr = np.zeros((3, 3), dtype=int)
        arr[np.tril_indices(3)] = np.array([a, b, c, d, e, f])
        return arr

    list_HNF = []
    for a in range(1, index + 1):
        if index % a != 0:
            continue
        for c in range(1, index // a + 1):
            if (index % c != 0) or (index % (a * c) != 0):
                continue
            f = index // (a * c)
            list_HNF.extend(
                [make_HNF(a, b, c, d, e, f) for b in range(c) for d in range(f) for e in range(f)]
            )

    return list_HNF


def reduce_HNF_list_by_parent_lattice_symmetry(
    list_HNF: List[np.ndarray], list_rotation_matrix: np.ndarray
) -> List[np.ndarray]:
    """
    reduce equivalent HNF with parent lattice symmetry

    Parameters
    ----------
    list_HNF: list of matrices, each element is Hermite normal form
    list_rotation_matrix: list of matrices
        each element represents the symmetry of parent lattice
    lattice_vector: array, (3, 3)
        lattice_vector[:, i] is the i-th lattice vector.

    Returns
    -------
    list_reduced_HNF: list of matrices, unique by symmetry
    """

    def is_equivalent(Bi, list_RBj_inv):
        if list_RBj_inv is None:
            return False
        list_H = np.dot(list_RBj_inv, Bi)
        close = np.isclose(list_H, np.around(list_H))
        if np.any(np.all(close, axis=(1, 2)), axis=-1):
            return True
        return False

    list_reduced_HNF = []
    list_RBj_inv = None

    sgn = np.linalg.det(list_rotation_matrix) == 1
    rotations = list_rotation_matrix[sgn, ...]

    for Bi in list_HNF:
        if not is_equivalent(Bi, list_RBj_inv):
            list_reduced_HNF.append(Bi)
            preinv = np.linalg.solve(np.dot(rotations, Bi), np.identity(Bi.shape[0])[None, ...])
            if list_RBj_inv is None:
                list_RBj_inv = preinv
            else:
                list_RBj_inv = np.concatenate([list_RBj_inv, preinv])

    return list_reduced_HNF


def generate_symmetry_distinct_superlattices(
    index: int,
    structure: Structure,
    return_symops=False,
    symprec=1e-2,
):
    """
    generate symmetry distict HNF

    Parameters
    ----------
    index: positive integer
    structure: pymatgen.core.Structure
    return_symops: bool

    Returns
    -------
    list_reduced_HNF: list of matrices, unique by symmetry
    (Optional) rotations: array, (# of symmetry operations, 3, 3)
    (Optional) translations: array, (# of symmetry operations, 3)
    """
    rotations, translations = get_symmetry_operations(structure, symprec=symprec)
    list_HNF = generate_all_superlattices(index)
    list_reduced_HNF = reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)
    if return_symops:
        return list_reduced_HNF, rotations, translations
    else:
        return list_reduced_HNF
