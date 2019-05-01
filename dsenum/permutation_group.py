from collections import namedtuple
from itertools import product
from typing import Optional, List

import numpy as np

from dsenum.smith_normal_form import smith_normal_form


DerivativeSite = namedtuple('DerivativeSite', ('site_index', 'jimage'))
CanonicalSite = namedtuple('CanonicalSite', ('site_index', 'factor'))


class DerivativeStructurePermutation(object):
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
    def __init__(self, hnf, displacement_set, rotations, translations):
        self.hnf = hnf
        self.num_sites_base = len(displacement_set)
        # TODO: check each site in displacement_set is in [0, 1)^dim
        self.displacement_set = displacement_set

        self.dhash = DerivativeMultiLatticeHash(self.hnf, self.displacement_set)

        self.rotations, self.translations = \
            self._get_superlattice_invariant_subgroup(rotations, translations)

        self.list_dsites = self.dhash.get_distinct_derivative_sites_list()
        self.list_csites = self.dhash.get_canonical_sites_list()

        self.prm_t = self._get_translation_permutations()
        self.prm_rigid = self._get_rigid_permutations()

    @property
    def dim(self):
        return self.dhash.dim

    @property
    def index(self):
        return self.dhash.index

    @property
    def num_site(self):
        return self.dhash.num_site

    def _get_superlattice_invariant_subgroup(self, rotations, translations):
        valid_rotations = []
        valid_translations = []

        for R, tau in zip(rotations, translations):
            if not is_same_lattice(np.dot(R, self.hnf), self.hnf):
                continue

            valid_rotations.append(R)
            valid_translations.append(tau)

        assert(len(rotations) % len(valid_rotations) == 0)
        return np.array(valid_rotations), np.array(valid_translations)

    def _get_translation_permutations(self):
        list_permutations = []
        for add_factor in self.dhash.get_all_factors():
            new_list_csites = []
            for site_index, factor in self.list_csites:
                new_factor = self.dhash.modulus_factor(np.array(factor) + np.array(add_factor))
                new_csite = CanonicalSite(site_index, new_factor)
                new_list_csites.append(new_csite)

            # permutation represenation
            perm = [self.dhash.hash_canonical_site(csite) for csite in new_list_csites]
            assert(is_permutation(perm))
            list_permutations.append(perm)

        # assume list_permutations[0] is identity
        assert(is_identity_permutation(list_permutations[0]))

        return list_permutations

    def _get_rigid_permutations(self):
        identity = list(range(self.num_site))
        list_permutations = [identity, ]

        for R, tau in zip(self.rotations, self.translations):
            new_list_csites = []
            for site_index, dimage in self.list_dsites:
                frac_coord = self.displacement_set[site_index] + np.array(dimage)
                acted_frac_coord = np.dot(R, frac_coord) + tau
                new_csite = self.dhash.hash_frac_coords(acted_frac_coord)
                assert(new_csite is not None)
                new_list_csites.append(new_csite)

            perm = [self.dhash.hash_canonical_site(csite) for csite in new_list_csites]
            assert(is_permutation(perm))
            if perm not in list_permutations:
                list_permutations.append(perm)

        # this set of permutations is not group!
        return list_permutations

    def get_symmetry_operation_permutaions(self):
        list_permutations = []

        for p1 in self.prm_t:
            for p2 in self.prm_rigid:
                perm = product_permutations(p1, p2)
                assert(perm not in list_permutations)
                list_permutations.append(perm)

        # assume list_permutations[0] is identity
        assert(is_identity_permutation(list_permutations[0]))

        return list_permutations


class DerivativeMultiLatticeHash(object):
    """
    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form(lower triangular)
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates in primitive cell of base structure
    """
    def __init__(self, hnf, displacement_set):
        self.hnf = hnf
        self.dim = self.hnf.shape[0]
        self.index = np.prod(self.hnf.diagonal())

        self.displacement_set = displacement_set

        self.num_site_base = len(displacement_set)
        self.num_site = self.num_site_base * self.index

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = cast_integer_matrix(np.linalg.inv(self.left))

        self.invariant_factors = tuple(self.snf.diagonal())
        self.shape = (self.num_site_base, ) + self.invariant_factors

    def hash_derivative_site(self, dsite: DerivativeSite) -> CanonicalSite:
        site_index, jimage = dsite.site_index, dsite.jimage

        factor_tmp = np.dot(self.left, np.array(jimage, dtype=int))
        factor = np.mod(factor_tmp, np.array(self.invariant_factors))

        derivative_jimage = cast_integer_matrix(factor_tmp - factor) \
            / np.array(self.invariant_factors)
        derivative_jimage = cast_integer_matrix(derivative_jimage)
        derivative_jimage = np.dot(self.right, derivative_jimage)

        csite = CanonicalSite(site_index, tuple(factor))
        return csite

    def hash_frac_coords(self, frac_coord) -> Optional[CanonicalSite]:
        for site_index, fc in enumerate(self.displacement_set):
            jimage = cast_integer_matrix(frac_coord - fc)
            if np.allclose(fc + jimage, frac_coord):
                dsite = DerivativeSite(site_index, jimage)
                csite = self.hash_derivative_site(dsite)
                return csite
        return None

    def get_canonical_sites_list(self) -> List[CanonicalSite]:
        list_csites = []
        for site_index in range(len(self.displacement_set)):
            for factor in product(*[range(f) for f in self.invariant_factors]):
                csite = CanonicalSite(site_index, tuple(factor))
                list_csites.append(csite)

        assert(all([self.hash_canonical_site(csite) == i
                    for i, csite in enumerate(list_csites)]))
        return list_csites

    def get_distinct_derivative_sites_list(self) -> List[DerivativeSite]:
        list_csites = self.get_canonical_sites_list()
        list_dsites = [self.embed_to_derivative_site(csite) for csite in list_csites]
        return list_dsites

    def get_canonical_and_derivative_sites_list(self):
        list_csites_dsites = []
        for site_index in range(len(self.displacement_set)):
            for factor in product(*[range(f) for f in self.invariant_factors]):
                csite = CanonicalSite(site_index, tuple(factor))
                dsite = self.embed_to_derivative_site(csite)
                list_csites_dsites.append((csite, dsite))

        return list_csites_dsites

    def get_all_factors(self):
        list_factors = list(product(*[range(f) for f in self.invariant_factors]))
        return list_factors

    def modulus_factor(self, factor):
        modded_factor = np.mod(factor, np.array(self.invariant_factors))
        return modded_factor

    def get_frac_coords(self, dsite: DerivativeSite):
        return self.displacement_set[dsite.site_index] + np.array(dsite.jimage)

    def hash_canonical_site(self, csite: CanonicalSite) -> int:
        # hash canonical site s.t. self.get_canonical_sites_list <=> identity
        multi_index = (csite.site_index, ) + tuple(csite.factor)
        raveled = np.ravel_multi_index(multi_index, self.shape)
        return raveled

    def unhash_indices_to_canonical_site(self, indices: int) -> CanonicalSite:
        unraveled = np.unravel_index(indices, self.shape)
        site_index, factor = unraveled[0], unraveled[1:]
        csite = CanonicalSite(site_index, factor)
        return csite

    def embed_to_derivative_site(self, csite: CanonicalSite) -> DerivativeSite:
        site_index, factor = csite
        jimage_base = cast_integer_matrix(np.dot(self.left_inv, factor))
        dsite = DerivativeSite(site_index, jimage_base)
        return dsite


def cast_integer_matrix(arr: np.ndarray) -> np.ndarray:
    arr_int = np.around(arr).astype(np.int)
    return arr_int


def is_unimodular(M: np.ndarray) -> bool:
    if np.abs(np.around(np.linalg.det(M))) == 1:
        return True
    else:
        return False


def is_same_lattice(H1: np.ndarray, H2: np.ndarray) -> bool:
    M = np.linalg.solve(H1, H2)
    M_int = cast_integer_matrix(M)

    if np.array_equal(np.dot(H1, M_int), H2):
        assert(is_unimodular(M_int))
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
        raise ValueError('not unique permutations')

    # identity
    if not any([is_identity_permutation(perm) for perm in list_permutations]):
        raise ValueError('There is not identity')

    # close
    for p1, p2 in product(list_permutations, repeat=2):
        p1p2 = product_permutations(p1, p2)
        if p1p2 not in list_permutations:
            raise ValueError('not closed in product')

    # inverse
    for perm in list_permutations:
        perm_inv = [-1 for _ in range(len(perm))]
        for i, idx in enumerate(perm):
            perm_inv[idx] = i

        if perm_inv not in list_permutations:
            raise ValueError('not contains inverse')

    return True


def product_permutations(p1, p2):
    """
    (p1 p2)(i) = p1(p2(i))
    """
    perm = [p1[p2[i]] for i in range(len(p1))]
    return perm
