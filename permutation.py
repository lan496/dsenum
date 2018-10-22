import itertools
import numpy as np

from smith_normal_form import smith_normal_form


class Permutation(object):
    """
    Parameters
    ----------
    hnf: array, (dims, dims)
        Hermite normal form
    rotations: list
        list of symmetry rotation operations in fractional coordinates

    Attributes
    ----------
    num: int
        the number of atoms in unit cell of super lattice
    dims: int
        dimentions of snf
    shapes: tuple of int
        diagonal of snf
    snf: array, (dims, dims)
        Smith normal form
    left: array, (dims, dims)
        left unimodular matrix
    left_inv: array, (dims, dims)
        inverse of left matrix
    right: array, (dims, dims)
        right unimodular matrix
    factors_e: array, (dims, num)
    images_e: array, (dims, num)
    """
    def __init__(self, hnf, rotations=None):
        self.hnf = hnf
        self.rotations = rotations

        self.num = np.prod(self.hnf.diagonal())
        self.dims = self.hnf.shape[0]

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.left_inv = np.around(np.linalg.inv(self.left))
        self.right = R

        self.shapes = tuple(self.snf.diagonal())

        self._init_unique_images()

    def _init_unique_images(self):
        # (dims, num)
        self.factors_e = np.array([
            np.unravel_index(indices, self.shapes)
            for indices in range(self.num)]).T
        # (dims, num)
        self.images_e = np.dot(self.left_inv, self.factors_e)

    def get_translation_permutations(self):
        list_permutations = []
        for i in range(self.num):
            di = self.factors_e[:, i]
            factors = np.mod(self.factors_e + di[:, np.newaxis],
                             np.array(self.shapes)[:, np.newaxis])
            raveled_factors = np.ravel_multi_index(factors, self.shapes)
            list_permutations.append(raveled_factors.tolist())
        return list_permutations

    def get_rotation_permutations(self):
        list_permutations = []
        if not self.rotations:
            return list_permutations

        for R in self.rotations:
            R = np.around(R)
            r_tmp = np.dot(self.left, np.dot(R, self.left_inv))
            factors = np.mod(np.dot(r_tmp, self.factors_e),
                             np.array(self.shapes)[:, np.newaxis])
            raveled_factors = np.ravel_multi_index(factors, self.shapes)
            list_permutations.append(raveled_factors.tolist())
        return list_permutations

    def get_symmetry_operation_permutaions(self):
        prm_t = self.get_translation_permutations()
        if not self.rotations:
            return prm_t

        prm_r = self.get_rotation_permutations()
        list_permutations = [self.product_permutations(pr, pt)
                             for pr, pt
                             in itertools.product(prm_r, prm_t)]
        return list_permutations

    def product_permutations(self, prm1, prm2):
        ret = [0 for _ in prm1]
        for i in range(self.num):
            ret[i] = prm2[prm1[i]]
        return ret


# https://repl.it/@smichr/msp
def msp(items):
    '''
    Yield the permutations of `items` where items is either a list
    of integers representing the actual items or a list of hashable items.
    The output are the unique permutations of the items given as a list
    of integers 0, ..., n-1 that represent the n unique elements in
    `items`.

    Examples
    ========

    >>> for i in msp('xoxox'):
    ...   print(i)

    [1, 1, 1, 0, 0]
    [0, 1, 1, 1, 0]
    [1, 0, 1, 1, 0]
    [1, 1, 0, 1, 0]
    [0, 1, 1, 0, 1]
    [1, 0, 1, 0, 1]
    [0, 1, 0, 1, 1]
    [0, 0, 1, 1, 1]
    [1, 0, 0, 1, 1]
    [1, 1, 0, 0, 1]

    Reference: "An O(1) Time Algorithm for Generating Multiset Permutations", Tadao Takaoka
    https://pdfs.semanticscholar.org/83b2/6f222e8648a7a0599309a40af21837a0264b.pdf
    '''

    def visit(head):
        (rv, j) = ([], head)
        for i in range(N):
            (dat, j) = E[j]
            rv.append(dat)
        return rv

    u = list(set(items))
    E = list(reversed(sorted([u.index(i) for i in items])))
    N = len(E)
    # put E into linked-list format
    (val, nxt) = (0, 1)
    for i in range(N):
        E[i] = [E[i], i + 1]
    E[-1][nxt] = None
    head = 0
    afteri = N - 1
    i = afteri - 1
    yield visit(head)
    while E[afteri][nxt] is not None or E[afteri][val] < E[head][val]:
        j = E[afteri][nxt]  # added to algorithm for clarity
        if j is not None and E[i][val] >= E[j][val]:
            beforek = afteri
        else:
            beforek = i
        k = E[beforek][nxt]
        E[beforek][nxt] = E[k][nxt]
        E[k][nxt] = head
        if E[k][val] < E[head][val]:
            i = k
        afteri = E[i][nxt]
        head = k
        yield visit(head)


if __name__ == '__main__':
    hnf = np.array([
        [2, 0],
        [1, 4]
    ])

    permutation = Permutation(hnf)
    prm_t = permutation.get_translation_permutations()
    print(permutation.shapes)
    print(permutation.snf)
    print(permutation.left)
    print(permutation.right)
    print(prm_t)
