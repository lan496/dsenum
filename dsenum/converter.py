from itertools import product
from typing import Optional, List, Tuple

import numpy as np

from dsenum.site import DerivativeSite, CanonicalSite
from dsenum.smith_normal_form import smith_normal_form
from dsenum.utils import cast_integer_matrix


class DerivativeMultiLatticeHash:
    """
                hash_frac_coords               ravel_canonical_site
    frac_coords ---------------> CanonicalSite -------------------------> Int
       |                            |^         <-------------------------
       |                            ||          unravel_to_canonical_site
       |                            ||
       |    embed_to_derivative_site||hash_derivative_site
       |                            ||
       |                            V|
       |-----------------------> DerivativeSite
            get_frac_coords

    Parameters
    ----------
    hnf: array, (dim, dim)
        Hermite normal form(lower triangular)
    displacement_set: array, (num_site_parent, dim)
        fractinal coordinates in primitive cell of base structure
    """

    def __init__(self, hnf: np.ndarray, displacement_set: np.ndarray):
        self._hnf = hnf
        self._index = np.around(np.abs(np.linalg.det(self.hnf))).astype(int)
        assert self.index != 0

        self.displacement_set = displacement_set
        assert self.dim == self.displacement_set.shape[1]

        self.num_site_base = len(displacement_set)

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = cast_integer_matrix(np.linalg.inv(self.left))

        self.invariant_factors = tuple(self.snf.diagonal())
        self.shape = (self.num_site_base,) + self.invariant_factors

    @property
    def hnf(self):
        return self._hnf

    @property
    def dim(self):
        return self.hnf.shape[0]

    @property
    def index(self):
        return self._index

    @property
    def num_sites(self):
        return self.num_site_base * self.index

    def hash_derivative_site(self, dsite: DerivativeSite, return_image=False) -> CanonicalSite:
        """
        Returns
        -------
        csite: CanoncalSite
        derivative_jimage: (Optional), array

        self.get_frac_coords(dsite) == self.displacement_set[dsite.site_index]
                                        + np.dot(self.left_inv, csite.factor)
                                        + np.dot(self.hnf, derivative_jimage)
        """
        site_index, jimage = dsite.site_index, dsite.jimage

        factor_tmp = np.dot(self.left, np.array(jimage, dtype=int))
        factor = np.mod(factor_tmp, np.array(self.invariant_factors))

        csite = CanonicalSite(site_index, tuple(factor))

        if return_image:
            derivative_jimage = cast_integer_matrix(factor_tmp - factor) / np.array(
                self.invariant_factors
            )
            derivative_jimage = cast_integer_matrix(derivative_jimage)
            derivative_jimage = np.dot(self.right, derivative_jimage)
            return csite, derivative_jimage
        else:
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

        assert all([self.ravel_canonical_site(csite) == i for i, csite in enumerate(list_csites)])
        return list_csites

    def get_distinct_derivative_sites_list(self) -> List[DerivativeSite]:
        list_csites = self.get_canonical_sites_list()
        list_dsites = [self.embed_to_derivative_site(csite) for csite in list_csites]
        return list_dsites

    def get_canonical_and_derivative_sites_list(
        self,
    ) -> List[Tuple[CanonicalSite, DerivativeSite]]:
        list_csites_dsites = []
        for site_index in range(len(self.displacement_set)):
            for factor in product(*[range(f) for f in self.invariant_factors]):
                csite = CanonicalSite(site_index, tuple(factor))
                dsite = self.embed_to_derivative_site(csite)
                list_csites_dsites.append((csite, dsite))

        return list_csites_dsites

    def get_all_factors(self) -> List[Tuple[int, ...]]:
        list_factors = list(product(*[range(f) for f in self.invariant_factors]))
        return list_factors

    def modulus_factor(self, factor: np.ndarray) -> np.ndarray:
        modded_factor = np.mod(factor, np.array(self.invariant_factors))
        return modded_factor

    def get_frac_coords(self, dsite: DerivativeSite) -> np.ndarray:
        return self.displacement_set[dsite.site_index] + np.array(dsite.jimage)

    def ravel_canonical_site(self, csite: CanonicalSite) -> int:
        # hash canonical site s.t. self.get_canonical_sites_list <=> identity
        multi_index = (csite.site_index,) + tuple(csite.factor)
        raveled = np.ravel_multi_index(multi_index, self.shape)
        return raveled

    def unravel_to_canonical_site(self, indices: int) -> CanonicalSite:
        unraveled = np.unravel_index(indices, self.shape)
        site_index, factor = unraveled[0], unraveled[1:]
        csite = CanonicalSite(site_index, factor)
        return csite

    def embed_to_derivative_site(self, csite: CanonicalSite) -> DerivativeSite:
        jimage_base = cast_integer_matrix(np.dot(self.left_inv, csite.factor))
        dsite = DerivativeSite(csite.site_index, tuple(jimage_base))
        return dsite


def get_species_list(index, list_species):
    """
    tile list of species for derivative structure
    """
    species = []
    for sp in list_species:
        species.extend([sp] * index)
    return species


def convert_site_constraints(base_site_constraints, index):
    site_constraints = []
    for sc in base_site_constraints:
        for _ in range(index):
            site_constraints.append(sc)
    return site_constraints
