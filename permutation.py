from itertools import zip_longest

import numpy as np

from smith_normal_form import smith_normal_form


class Permutation(object):
    """
    Parameters
    ----------
    hnf: array, (dims, dims)
        Hermite normal form
    num_site_parent: int
        # of atoms in parent multilattice
    displacement_set: array, (num_site_parent, 3)
        fractinal coordinates of A in multilattice site
    rotations: array, (# of symmetry operations, 3, 3)
        rotations in fractinal coordinations of A
    translations: array, (# of symmetry operations, 3)
        translations in fractinal coordinations of A

    Attributes
    ----------
    num: int
        the number of atoms in unit cell of super lattice
    shape: tuple of int
        (num_site_parent, D11, D22, D33)
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
    def __init__(self, hnf, num_site_parent=1, displacement_set=None,
                 rotations=None, translations=None):
        self.hnf = hnf
        self.hnf_inv = np.linalg.inv(self.hnf)
        self.index = np.prod(self.hnf.diagonal())

        self.num_site_parent = num_site_parent
        self.num_site = self.index * self.num_site_parent

        if self.num_site_parent == 1:
            displacement_set = np.array([[0, 0, 0]])
        else:
            self.displacement_set = displacement_set
            assert self.displacement_set.shape[0] == self.num_site_parent

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)
        self.right = R

        self.shape = tuple([self.num_site_parent]
                           + self.snf.diagonal().tolist())

        # (num_site_parent, D11, D22, D33)
        self.factors_e = np.array([np.unravel_index(indices, self.shape)
                                   for indices in range(self.num_site)]).T
        self.rotations, self.translations = \
            self._get_superlattice_symmetry_operations(rotations, translations)

    def _get_superlattice_symmetry_operations(self,
                                              rotations, translations):
        """
        return symmetry operations of superlattice
        """
        valid_rotations = []
        valid_translations = []
        self.rotation_factors = []

        for R, tau in zip_longest(rotations, translations):
            if not is_same_lattice(np.dot(R, self.hnf), self.hnf):
                continue

            """
            r_tmp = np.dot(self.left, np.dot(R, self.left_inv))
            factors = np.dot(r_tmp, self.factors_e)
            factors = np.mod(factors,
                             np.array(self.shape)[:, np.newaxis])
            self.rotation_factors.append(factors)
            """
            valid_rotations.append(R)
            valid_translations.append(tau)

        return valid_rotations, valid_translations

    def get_translation_permutations(self):
        list_permutations = []
        for i in range(self.num):
            di = np.dot(self.left, self.factors_e[:, i])
            factors = np.mod(self.factors_e + di[:, np.newaxis],
                             np.array(self.shape)[:, np.newaxis])
            raveled_factors = np.ravel_multi_index(factors, self.shape)
            list_permutations.append(tuple(raveled_factors.tolist()))
        assert self.validate_permutations(list_permutations)
        return list_permutations

    def get_symmetry_operation_permutaions(self):
        if self.rotations is None:
            return prm_t

        list_permutations = []

        for i in range(self.num):
            for factors_r in self.rotation_factors:
                di = self.factors_e[:, i]
                factors = np.mod(factors_r + di[:, np.newaxis],
                                 np.array(self.shape)[:, np.newaxis])
                raveled_factors = tuple(np.ravel_multi_index(factors, self.shape).tolist())
                if raveled_factors in list_permutations:
                    continue
                if len(set(raveled_factors)) != len(raveled_factors):
                    continue
                list_permutations.append(raveled_factors)

        assert self.validate_permutations(list_permutations)
        return list_permutations

    def validate_permutations(self, permutations):
        for prm in permutations:
            if len(set(prm)) != self.num:
                print(prm)
                print('not permutation')
                return False

        if len(set(permutations)) != len(permutations):
            print('not unique')
            print(permutations)
            return False

        return True


def is_same_lattice(H1, H2):
    M = np.linalg.solve(H1, H2)
    M_re = np.around(M).astype(np.int)

    if np.array_equal(np.dot(H1, M_re), H2):
        assert np.around(np.linalg.det(M_re) ** 2) == 1
        return True
    else:
        return False


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
        [2, 4]
    ])

    rotations = np.array([
        [[1, 0], [0, 1]],
        [[0, -1], [1, 0]],
        [[-1, 0], [0, -1]],
        [[0, 1], [-1, 0]],
    ])

    permutation = Permutation(hnf, rotations)
    prm_t = permutation.get_translation_permutations()

    print(permutation.snf)
    print(permutation.left)
    print(permutation.right)
    print()

    print('factors_e')
    print(permutation.factors_e)
    print('images_e')
    print(permutation.images_e)
