import numpy as np

from dsenum.smith_normal_form import smith_normal_form


class Permutation(object):
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form
    num_site_parent: int
        # of atoms in parent multilattice
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates of A in multilattice site
    rotations: array, (# of symmetry operations, dim, dim)
        rotations in fractinal coordinations of A
    translations: array, (# of symmetry operations, dim)
        translations in fractinal coordinations of A

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    num_site: int
        # of sites in unit cell of superlattice
    shape: tuple of int
        (1 + dim, num_site)
    snf: array, (dim, dim)
        Smith normal form
    left: array, (dim, dim)
        left unimodular matrix
    left_inv: array, (dim, dim)
        inverse of left matrix
    right: array, (dim, dim)
        right unimodular matrix
    factors_e: array, (1 + dim, num_site)
    parent_frac_coords_e: array (dim, num_site)
    """
    def __init__(self, hnf, num_site_parent=1, displacement_set=None,
                 rotations=None, translations=None):
        self.hnf = hnf
        self.num_site_parent = num_site_parent
        self.dshash = DerivativeStructureHash(hnf, num_site_parent)

        if self.num_site_parent == 1:
            self.displacement_set = np.array([[0, 0, 0]])
        else:
            self.displacement_set = displacement_set
            assert self.displacement_set.shape[0] == self.num_site_parent

        """
        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)
        self.right = R
        """

        # (1 + dim, num_site)
        self.factors_e = self.dshash.get_distinct_factors().T
        # self.factors_e = np.array([np.unravel_index(indices, self.shape)
        #                            for indices in range(self.num_site)]).T

        # (dim, num_site)
        self.parent_frac_coords_e = self.displacement_set[self.factors_e[0, :]].T \
            + np.dot(self.left_inv, self.factors_e[1:, :])

        self.rotations, self.translations = \
            self._get_superlattice_symmetry_operations(rotations, translations)

        self.prm_t = self.get_translation_permutations()
        self.prm_rigid = self.get_rigid_permutations()

    @property
    def dim(self):
        return self.dshash.dim

    @property
    def index(self):
        return self.dshash.index

    @property
    def num_site(self):
        return self.dshash.num_site

    @property
    def shape(self):
        return self.dshash.shape

    @property
    def left(self):
        return self.dshash.left

    @property
    def left_inv(self):
        return self.dshash.left_inv

    def _get_superlattice_symmetry_operations(self, rotations, translations):
        """
        return symmetry operations of superlattice
        """
        valid_rotations = []
        valid_translations = []

        if rotations is not None:
            for i, R in enumerate(rotations):
                if not is_same_lattice(np.dot(R, self.hnf), self.hnf):
                    continue

                valid_rotations.append(R)
                if translations is not None:
                        valid_translations.append(translations[i])

        return np.array(valid_rotations), np.array(valid_translations)

    def get_rigid_permutations(self, eps=1e-10):
        prm_rigid = []
        self.rigid_factors_ = []

        for R, tau in zip(self.rotations, self.translations):
            # (dim, num_site)
            parent_frac_coords = np.dot(R, self.parent_frac_coords_e) \
                + tau[:, np.newaxis]
            dset = np.remainder(parent_frac_coords + eps, 1)
            lattice_factors = np.dot(self.left,
                                     np.around(parent_frac_coords - dset).astype(np.int))
            lattice_factors = np.mod(lattice_factors,
                                     np.array(self.shape[1:])[:, np.newaxis])

            factors = -np.ones_like(self.factors_e)
            factors[1:, :] = lattice_factors
            # TODO: awkward
            for i in range(self.num_site):
                for j, ds in enumerate(self.displacement_set):
                    if np.allclose(dset[:, i], ds, rtol=1e-5):
                        factors[0, i] = j
                        break

            if np.any(factors[0] == -1):
                continue
            # assert (np.sum(factors == -1) == 0)

            raveled_factors = tuple(np.ravel_multi_index(factors, self.shape).tolist())
            if raveled_factors not in prm_rigid:
                prm_rigid.append(raveled_factors)
                self.rigid_factors_.append(factors)

        return prm_rigid

    def get_translation_permutations(self):
        list_permutations = []
        for i in range(self.index):
            di = self.factors_e[1:, i]
            factors = np.copy(self.factors_e)
            factors[1:, :] = np.mod(self.factors_e[1:, :] + di[:, np.newaxis],
                                    np.array(self.shape[1:])[:, np.newaxis])
            raveled_factors = np.ravel_multi_index(factors, self.shape)
            list_permutations.append(tuple(raveled_factors.tolist()))
        return list_permutations

    def get_symmetry_operation_permutaions(self):
        if self.rotations is None:
            return self.prm_t

        list_permutations = []

        for i in range(self.index):
            for factors_r in self.rigid_factors_:
                di = self.factors_e[1:, i]
                factors = np.copy(factors_r)
                factors[1:, :] = np.mod(factors[1:, :] + di[:, np.newaxis],
                                        np.array(self.shape[1:])[:, np.newaxis])
                raveled_factors = tuple(np.ravel_multi_index(factors, self.shape).tolist())

                if raveled_factors in list_permutations:
                    continue
                # assert len(set(raveled_factors)) == len(raveled_factors)
                list_permutations.append(raveled_factors)

        return list_permutations


class DerivativeStructureHash:

    def __init__(self, hnf, num_site_parents=1):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())

        self.num_site_parents = num_site_parents
        self.num_site = self.num_site_parents * self.index

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = np.around(np.linalg.inv(self.left)).astype(np.int)

        self.shape = tuple([self.num_site_parents]
                           + self.snf.diagonal().tolist())

    def hash_fractional_coordinates(self, indexes):
        """
        (d, n) -> (d, Ln mod D)

        Parameters
        ----------
        indexes: array, ([N,] 1 + self.dim)

        Returns
        -------
        factors: array, ([N,] 1 + self.dim)
        """
        indexes = np.atleast_2d(indexes)
        factors = np.empty_like(indexes)
        factors[:, 0] = indexes[:, 0][:]
        factors[:, 1:] = np.sum(self.left[np.newaxis, :, :] * indexes[:, np.newaxis, 1:],
                                axis=2)
        factors[:, 1:] = np.mod(factors[:, 1:], np.array(self.shape[1:])[np.newaxis, :])
        if factors.shape[0] == 1:
            factors = np.squeeze(factors, axis=0)
        return factors

    def unhash_factors(self, factors):
        """
        f -> L^{-1}f

        Parameters
        ----------
        factors: array, ([N,] 1 + self.dim)

        Returns
        -------
        indexes: array, ([N,] 1 + self.dim)
        """
        factors = np.atleast_2d(factors)
        indexes = np.empty_like(factors)
        indexes[:, 0] = factors[:, 0][:]
        indexes[:, 1:] = np.sum(self.left_inv[np.newaxis, :, :]
                                * factors[:, np.newaxis, 1:], axis=2)
        if indexes.shape[0] == 1:
            indexes = np.squeeze(indexes, axis=0)
        return indexes

    def get_distinct_factors(self):
        """
        return a list of factors corresponding to distinct lattice points

        Parameters
        ----------

        Returns
        -------
        factors_e: array, (self.num_site, 1 + self.dim)
        """
        factors_e = np.array([np.unravel_index(indices, self.shape)
                              for indices in range(self.num_site)])
        return factors_e

    def get_distinct_indexes(self):
        """
        return a list of lattice indexes corresponding to distinct lattice points

        Parameters
        ----------

        Returns
        -------
        indexes_e: array, (self.num_site, 1 + self.dim)
        """
        indexes_e = self.unhash_factors(self.get_distinct_factors())
        return indexes_e


def is_same_lattice(H1, H2):
    M = np.linalg.solve(H1, H2)
    M_re = np.around(M).astype(np.int)

    if np.array_equal(np.dot(H1, M_re), H2):
        assert np.around(np.linalg.det(M_re) ** 2) == 1
        return True
    else:
        return False
