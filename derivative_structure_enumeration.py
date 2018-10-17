import numpy as np


def generate_all_superlattices(index):
    """
    Parameters
    ----------
    index: positive integer

    Returns
    -------
    list_HNF: list of matrices, each element is Hermite normal form
    """
    def make_HNF(a, b, c, d, e, f):
        arr = np.zeros((3, 3), dtype=np.int)
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
            list_HNF.extend([make_HNF(a, b, c, d, e, f)
                             for b in range(c) for d in range(f) for e in range(f)])

    return list_HNF


def reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, list_rotation_matrix, n_jobs=1):
    """
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
        for RBj_inv in list_RBj_inv:
            H = np.dot(RBj_inv, Bi)
            if np.allclose(H, np.around(H)):
                return True
        return False

    list_reduced_HNF = []
    list_RBj_inv = []

    for Bi in list_HNF:
        if not is_equivalent(Bi, list_RBj_inv):
            list_reduced_HNF.append(Bi)
            preinv = [np.linalg.solve(np.dot(R, Bi), np.identity(R.shape[0])) for R
                      in list_rotation_matrix]
            list_RBj_inv.extend(preinv)

    return list_reduced_HNF
