cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
def hash_in_all_configuration(coloring, num_color):
    ret = 0
    cdef int e
    for e in coloring:
        ret = ret * num_color + e
    return ret


@cython.boundscheck(False)
@cython.wraparound(False)
def act_permutation(perm, coloring):
    new_coloring = [coloring[perm[i]] for i in range(len(coloring))]
    return new_coloring


@cython.boundscheck(False)
@cython.wraparound(False)
def get_composition(coloring, num_color):
    composition = [0 for _ in range(num_color)]
    for c in coloring:
        composition[c] += 1
    return composition
