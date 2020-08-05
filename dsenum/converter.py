from itertools import product
from typing import Optional, List
from collections import namedtuple

import numpy as np

from dsenum.smith_normal_form import smith_normal_form
from dsenum.utils import cast_integer_matrix


DerivativeSite = namedtuple("DerivativeSite", ("site_index", "jimage"))
CanonicalSite = namedtuple("CanonicalSite", ("site_index", "factor"))


class DerivativeMultiLatticeHash:
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
        self.index = np.around(np.linalg.det(hnf)).astype(int)
        assert self.index != 0

        self.displacement_set = displacement_set

        self.num_site_base = len(displacement_set)
        self.num_site = self.num_site_base * self.index

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R
        self.left_inv = cast_integer_matrix(np.linalg.inv(self.left))

        self.invariant_factors = tuple(self.snf.diagonal())
        self.shape = (self.num_site_base,) + self.invariant_factors

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

        assert all([self.hash_canonical_site(csite) == i for i, csite in enumerate(list_csites)])
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
        multi_index = (csite.site_index,) + tuple(csite.factor)
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

    def get_species_list(self, list_species):
        """
        tile list of species for derivative structure
        """
        species = []
        for sp in list_species:
            species.extend([sp] * self.index)
        return species

    @classmethod
    def convert_site_constraints(cls, base_site_constraints, index):
        site_constraints = []
        for sc in base_site_constraints:
            for _ in range(index):
                site_constraints.append(sc)
        return site_constraints
