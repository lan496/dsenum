from itertools import product

import numpy as np

from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.site import CanonicalSite
from dsenum.utils import cast_integer_matrix


class DerivativeStructurePermutation:
    """
    Permutation Representation of space group of superlattice

    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form(lower triangular)
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates in primitive cell of base structure
    rotations: array, (# of symmetry operations, dim, dim)
        rotations with primitive basis for base structure
    translations: array, (# of symmetry operations, dim)
        translations with primitive basis for base structure

    Attributes
    ----------
    dim : int
        dimention of lattice
    index: int
        # of parent multilattice in super lattice
    num_site: int
        # of sites in unit cell of superlattice
    """

    def __init__(
        self,
        hnf: np.ndarray,
        displacement_set: np.ndarray,
        rotations: np.ndarray,
        translations: np.ndarray,
    ):
        self.hnf = hnf
        self.num_sites_base = len(displacement_set)
        # TODO: check each site in displacement_set is in [0, 1)^dim
        self.displacement_set = displacement_set

        self.dhash = DerivativeMultiLatticeHash(self.hnf, self.displacement_set)

        self.rotations, self.translations = self._get_superlattice_invariant_subgroup(
            rotations, translations
        )

        self.list_dsites = self.dhash.get_distinct_derivative_sites_list()
        self.list_csites = self.dhash.get_canonical_sites_list()

        self._prm_t = self._get_translation_permutations()
        self.prm_rigid = self._get_rigid_permutations()

    @property
    def dim(self):
        return self.dhash.dim

    @property
    def index(self):
        return self.dhash.index

    @property
    def num_sites(self):
        return self.dhash.num_sites

    @property
    def prm_t(self):
        return self._prm_t

    def _get_superlattice_invariant_subgroup(
        self, rotations: np.ndarray, translations: np.ndarray
    ):
        valid_rotations = []
        valid_translations = []

        for R, tau in zip(rotations, translations):
            if not is_same_lattice(np.dot(R, self.hnf), self.hnf):
                continue

            valid_rotations.append(R)
            valid_translations.append(tau)

        assert len(rotations) % len(valid_rotations) == 0
        return np.array(valid_rotations), np.array(valid_translations)

    def _get_translation_permutations(self):
        list_permutations = []
        for add_factor in self.dhash.get_all_factors():
            new_list_csites = []
            for csite in self.list_csites:
                new_factor = self.dhash.modulus_factor(
                    np.array(csite.factor) + np.array(add_factor)
                )
                new_csite = CanonicalSite(csite.site_index, new_factor)
                new_list_csites.append(new_csite)

            # permutation represenation
            perm = [self.dhash.ravel_canonical_site(csite) for csite in new_list_csites]
            assert is_permutation(perm)
            list_permutations.append(perm)

        # assume list_permutations[0] is identity
        assert is_identity_permutation(list_permutations[0])

        return list_permutations

    def _get_rigid_permutations(self):
        identity = list(range(self.num_sites))
        list_permutations = [
            identity,
        ]

        for R, tau in zip(self.rotations, self.translations):
            new_list_csites = []
            for dsite in self.list_dsites:
                frac_coord = self.displacement_set[dsite.site_index] + np.array(dsite.jimage)
                acted_frac_coord = np.dot(R, frac_coord) + tau
                new_csite = self.dhash.hash_frac_coords(acted_frac_coord)
                assert new_csite is not None
                new_list_csites.append(new_csite)

            perm = [self.dhash.ravel_canonical_site(csite) for csite in new_list_csites]
            assert is_permutation(perm)
            if perm not in list_permutations:
                list_permutations.append(perm)

        # this set of permutations is not group!
        return list_permutations

    def get_symmetry_operation_permutations(self):
        list_permutations = []

        for p1 in self.prm_t:
            for p2 in self.prm_rigid:
                perm = product_permutations(p1, p2)
                assert perm not in list_permutations
                list_permutations.append(perm)

        # assume list_permutations[0] is identity
        assert is_identity_permutation(list_permutations[0])

        return list_permutations


def is_unimodular(M: np.ndarray) -> bool:
    if np.abs(np.around(np.linalg.det(M))) == 1:
        return True
    else:
        return False


def is_same_lattice(H1: np.ndarray, H2: np.ndarray) -> bool:
    M = np.linalg.solve(H1, H2)
    M_int = cast_integer_matrix(M)

    if np.array_equal(np.dot(H1, M_int), H2):
        assert is_unimodular(M_int)
        return True
    else:
        return False


def is_permutation(perm):
    return len(set(perm)) == len(perm)


def is_identity_permutation(perm):
    if all([index == i for i, index in enumerate(perm)]):
        return True
    else:
        return False


def is_permutation_group(list_permutations):
    list_permutations_tuple = [tuple(perm) for perm in list_permutations]
    if len(set(list_permutations_tuple)) != len(list_permutations_tuple):
        print(list_permutations)
        raise ValueError("not unique permutations")

    # identity
    if not any([is_identity_permutation(perm) for perm in list_permutations]):
        raise ValueError("There is not identity")

    # close
    for p1, p2 in product(list_permutations, repeat=2):
        p1p2 = product_permutations(p1, p2)
        if p1p2 not in list_permutations:
            raise ValueError("not closed in product")

    # inverse
    for perm in list_permutations:
        perm_inv = [-1 for _ in range(len(perm))]
        for i, idx in enumerate(perm):
            perm_inv[idx] = i

        if perm_inv not in list_permutations:
            raise ValueError("not contains inverse")

    return True


def product_permutations(p1, p2):
    """
    (p1 p2)(i) = p1(p2(i))
    """
    perm = [p1[p2[i]] for i in range(len(p1))]
    return perm
