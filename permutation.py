import numpy as np

from smith_normal_form import smith_normal_form


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
        self.hnf_inv = np.linalg.inv(self.hnf)
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())

        self.num_site_parent = num_site_parent
        self.num_site = self.index * self.num_site_parent

        if self.num_site_parent == 1:
            self.displacement_set = np.array([[0, 0, 0]])
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

        # (1 + dim, num_site)
        self.factors_e = np.array([np.unravel_index(indices, self.shape)
                                   for indices in range(self.num_site)]).T
        # (dim, num_site)
        self.parent_frac_coords_e = self.displacement_set[self.factors_e[0, :]].T \
            + np.dot(self.left_inv, self.factors_e[1:, :])

        self.rotations, self.translations = \
            self._get_superlattice_symmetry_operations(rotations, translations)

        self.prm_t = self.get_translation_permutations()
        self.prm_rigid = self.get_rigid_permutations()

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
                assert len(set(raveled_factors)) == len(raveled_factors)
                list_permutations.append(raveled_factors)

        return list_permutations


def is_same_lattice(H1, H2):
    M = np.linalg.solve(H1, H2)
    M_re = np.around(M).astype(np.int)

    if np.array_equal(np.dot(H1, M_re), H2):
        assert np.around(np.linalg.det(M_re) ** 2) == 1
        return True
    else:
        return False
