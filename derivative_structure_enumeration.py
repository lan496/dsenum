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
