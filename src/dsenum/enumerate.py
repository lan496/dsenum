from __future__ import annotations

from abc import ABCMeta, abstractmethod
from time import time
from typing import Literal, cast

import numpy as np
from numpy.typing import NDArray
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.util.typing import SpeciesLike
from pyzdd import Permutation, Universe
from pyzdd.structure import construct_derivative_structures, enumerate_labelings
from tqdm import tqdm

from dsenum.coloring import SiteColoringEnumerator
from dsenum.coloring_generator import (
    BaseColoringGenerator,
    ColoringGenerator,
    FixedConcentrationColoringGenerator,
    ListBasedColoringGenerator,
)
from dsenum.converter import convert_site_constraints
from dsenum.derivative_structure import ColoringToStructure
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.utils import get_symmetry_operations


class AbstractStructureEnumerator(metaclass=ABCMeta):
    """
    Abstract class for enumerating derivative structures.
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        mapping_color_species: list[SpeciesLike] = None,
        composition_constraints=None,
        base_site_constraints=None,
        color_exchange=True,
        remove_superperiodic=True,
        remove_incomplete=True,
        verbose=True,
    ):
        self.base_structure = base_structure
        self.index = index
        self.num_types = num_types
        self.composition_constraints = composition_constraints

        # settings
        self.color_exchange = color_exchange
        self.remove_superperiodic = remove_superperiodic
        self.remove_incomplete = remove_incomplete

        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, base_structure, return_symops=True
        )
        self.list_reduced_HNF = list_reduced_HNF
        self.rotations = rotations
        self.translations = translations

        # self.site_constraints[i] is a list of allowed species at the i-th site in supercell
        self.site_constraints = None
        if base_site_constraints:
            assert len(base_site_constraints) == self.num_sites_base
            self.site_constraints = convert_site_constraints(base_site_constraints, self.index)

        if mapping_color_species and len(mapping_color_species) != self.num_types:
            raise ValueError("mapping_color_species must have num_type species.")
        if mapping_color_species is None:
            mapping_color_species = [DummySpecie(str(i)) for i in range(1, self.num_types + 1)]
        self.mapping_color_species = mapping_color_species

        self.verbose = verbose

    @property
    def num_sites_base(self):
        return self.base_structure.num_sites

    @property
    def num_sites(self):
        return self.num_sites_base * self.index

    def generate(
        self,
        return_colorings: bool = False,
        additional_species: list[SpeciesLike] | None = None,
        additional_frac_coords: NDArray | None = None,
        output: Literal["poscar", "pymatgen"] = "pymatgen",
    ) -> list[Structure | str] | tuple[list[Structure | str], list[NDArray], list[list[int]]]:
        """
        Generate derivative structures

        Parameters
        ----------
        return_colorings: bool, optional
            If true, return transformation matrices and colorings in addition
        additional_species: list[SpeciesLike] | None, optional
            species which are nothing to do with ordering. If specified, append these species in
            returned derivative structures.
        additional_frac_coords: NDArray | None, optional
            fractional coordinates of species which are nothing to do with ordering
        output: Literal['poscar', 'pymatgen']
            In default, `pymatgen.core.Structure` objects are returned. If specified
            `output='poscar'`, POSCAR strings are returned instead.

        Returns
        -------
        list_ds: list[Structure | str]
            List of derivative structures
        list[list[int]]

        list_transformations: list[NDArray], optional
            List of transformation matrices. Returned if `return_colorings=True`
        list_colorings: list[list[int]], optional
            List of colorings. Returned if `return_colorings=True`
        """
        assert (output == "pymatgen") or (output == "poscar")
        start = time()

        displacement_set = self.base_structure.frac_coords
        list_ds = []
        list_transformations = []
        list_colorings = []
        for hnf in tqdm(self.list_reduced_HNF, disable=not self.verbose):
            ds_permutation = DerivativeStructurePermutation(
                hnf, displacement_set, self.rotations, self.translations
            )
            # enumerate colorings
            list_colorings_hnf = self._generate_coloring_with_hnf(
                hnf, ds_permutation, additional_species, additional_frac_coords
            )

            # convert to Structure object
            cts = ColoringToStructure(
                self.base_structure,
                ds_permutation.dhash,
                self.mapping_color_species,
                additional_species=additional_species,
                additional_frac_coords=additional_frac_coords,
            )
            if output == "pymatgen":
                list_ds_hnf = [cts.convert_to_structure(cl) for cl in list_colorings_hnf]
            elif output == "poscar":
                list_ds_hnf = [cts.convert_to_poscar_string(cl) for cl in list_colorings_hnf]

            list_ds.extend(list_ds_hnf)
            if return_colorings:
                list_transformations.extend([hnf] * len(list_ds_hnf))
                list_colorings.extend(list_colorings_hnf)

        end = time()
        if self.verbose:
            print("total: {} (Time: {:.4}sec)".format(len(list_ds), end - start))

        if return_colorings:
            return list_ds, list_transformations, list_colorings
        else:
            return list_ds

    @abstractmethod
    def _generate_coloring_with_hnf(
        self,
        hnf: np.ndarray,
        ds_permutation: DerivativeStructurePermutation,
        additional_species,
        additional_frac_coords,
    ) -> list[list[int]]:
        raise NotImplementedError


class StructureEnumerator(AbstractStructureEnumerator):
    """
    Enumerate derivative structures.

    Parameters
    ----------
    base_structure: pymatgen.core.Structure
        Aristotype for derivative structures
    index: int
        How many times to expand unit cell
    num_types: int
        The number of species in derivative structures.
        `num_types` may be larger than the number of the kinds of species in
        `base_structure`: for example, you consider vacancies in derivative structures.
    mapping_color_species: list[SpecieLike] | None, optional
        If specified, use these species in derivative structures.
        The length of this list should be equal to `num_types`
    composition_constraints: list[int] | None, optional
        composition_constraints[i] is the ratio of the i-th species in mapping_color_species.
        For example, set `composition_constraints = [1, 2]` for enumerating TiO2 structures with
        `mapping_color_species = ["Ti", "O"]`.
    base_site_constraints: list[list[int]] | None, optional
        (num_elements, num_color) e.g. site_constraints[2] = [0, 3, 4] means color of site-2
        in base_structure must be 0, 3, or 4.
    color_exchange: bool, optional
        Iff true, identify color-exchanging
    remove_superperiodic: bool, optional
        Iff true, discard superperiodic coloring
    remove_incomplete: bool, optional
        Iff true, discard structures whose number of types are less then `num_types`.
    method: str, optional
        "direct" or "lexicographic", so far
    n_jobs: int, optional
        core in lexicographic coset enumeration(only used when method='lexicographic')
    verbose: bool, optional
        If true, print progress and number of enumerated structures.
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        mapping_color_species: list[SpeciesLike] | None = None,
        composition_constraints: list[int] | None = None,
        base_site_constraints: list[list[int]] | None = None,
        color_exchange: bool = True,
        remove_superperiodic: bool = True,
        remove_incomplete: bool = True,
        method: Literal["direct", "lexicographic"] = "direct",
        n_jobs: int = 1,
        verbose: bool = True,
    ):
        super().__init__(
            base_structure=base_structure,
            index=index,
            num_types=num_types,
            mapping_color_species=mapping_color_species,
            composition_constraints=composition_constraints,
            base_site_constraints=base_site_constraints,
            color_exchange=color_exchange,
            remove_superperiodic=remove_superperiodic,
            remove_incomplete=remove_incomplete,
            verbose=verbose,
        )
        self.method = method
        self.n_jobs = n_jobs

        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, base_structure, return_symops=True
        )
        self.list_reduced_HNF = list_reduced_HNF
        self.rotations = rotations
        self.translations = translations

        # composition constraints
        # typing.cast causes no runtime effect
        if self.composition_constraints is None:
            cl_generator = cast(
                BaseColoringGenerator,
                ColoringGenerator(
                    self.num_sites, self.num_types, site_constraints=self.site_constraints
                ),
            )
        else:
            cl_generator = cast(
                BaseColoringGenerator,
                FixedConcentrationColoringGenerator(
                    self.num_sites,
                    self.num_types,
                    self.composition_constraints,
                    site_constraints=self.site_constraints,
                ),
            )
        self.cl_generator = cl_generator

    def _generate_coloring_with_hnf(
        self,
        hnf: np.ndarray,
        ds_permutation: DerivativeStructurePermutation,
        additional_species,
        additional_frac_coords,
    ) -> list[list[int]]:
        sc_enum = SiteColoringEnumerator(
            self.num_types,
            ds_permutation,
            self.cl_generator,
            self.color_exchange,
            self.remove_superperiodic,
            self.remove_incomplete,
            method=self.method,
            n_jobs=self.n_jobs,
        )
        colorings = sc_enum.unique_colorings()

        return colorings


class ZddStructureEnumerator(AbstractStructureEnumerator):
    """
    Enumerate derivative structures with ZDD acceleration.

    Parameters
    ----------
    base_structure: pymatgen.core.Structure
        Aristotype for derivative structures
    index: int
        How many times to expand unit cell
    num_types: int
        The number of species in derivative structures.
        `num_types` may be larger than the number of the kinds of species
        in `base_structure`: for example, you consider vacancies in derivative structures.
    mapping_color_species: list[int] | None, optional
        If specified, use these species in derivative structures.
        The length of this list should be equal to `num_types`
    composition_constraints: list[int] | None, optional
        composition_constraints[i] is the ratio of the i-th species in mapping_color_species.
        For example, set `composition_constraints = [1, 2]` for enumerating TiO2 structures with
        `mapping_color_species = ["Ti", "O"]`.
    base_site_constraints: list[list[int]] | None, optional
        (num_elements, num_color) e.g. site_constraints[2] = [0, 3, 4] means color of site-2
        in base_structure must be 0, 3, or 4.
    remove_superperiodic: bool, optional
        Iff true, discard superperiodic coloring
    remove_incomplete: bool, optional
        Iff true, discard structures whose number of types are less then `num_types`.
    verbose: bool, optional
        If true, print progress and number of enumerated structures.
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        mapping_color_species: list[SpeciesLike] | None = None,
        composition_constraints: list[int] | None = None,
        base_site_constraints: list[list[int]] | None = None,
        remove_superperiodic: bool = True,
        remove_incomplete: bool = True,
        verbose: bool = True,
    ):
        super().__init__(
            base_structure=base_structure,
            index=index,
            num_types=num_types,
            mapping_color_species=mapping_color_species,
            composition_constraints=composition_constraints,
            base_site_constraints=base_site_constraints,
            remove_superperiodic=remove_superperiodic,
            remove_incomplete=remove_incomplete,
            verbose=verbose,
        )

        self.prohibited_site_constraints = []
        if self.site_constraints is not None:
            all_species: set[int] = set()
            all_species.update(range(self.num_types))
            for allowed_species in self.site_constraints:
                prohibited = list(all_species.difference(set(allowed_species)))
                self.prohibited_site_constraints.append(prohibited)
            assert len(self.prohibited_site_constraints) == self.num_sites

    def _generate_coloring_with_hnf(
        self,
        hnf: np.ndarray,
        ds_permutation: DerivativeStructurePermutation,
        additional_species,
        additional_frac_coords,
    ) -> list[list[int]]:
        dd = Universe()

        num_sites = ds_permutation.num_sites
        automorphism = [
            Permutation(sigma) for sigma in ds_permutation.get_symmetry_operation_permutations()
        ]
        translations = [Permutation(sigma) for sigma in ds_permutation._prm_t]

        composition_constraints_dd: list[int] = []
        if self.composition_constraints is not None:
            ratio_sum = np.sum(self.composition_constraints)
            if num_sites % ratio_sum != 0:
                # impossible to satisfy composition constraints
                return []
            else:
                composition_constraints_dd = [
                    ratio * (num_sites // ratio_sum) for ratio in self.composition_constraints
                ]

        construct_derivative_structures(
            dd,
            num_sites=num_sites,
            num_types=self.num_types,
            automorphism=automorphism,
            translations=translations,
            composition_constraints=composition_constraints_dd,
            site_constraints=self.prohibited_site_constraints,
            remove_incomplete=self.remove_incomplete,
            remove_superperiodic=self.remove_superperiodic,
        )

        colorings = list(enumerate_labelings(dd, num_sites, self.num_types))
        return colorings

    def count(self) -> int:
        """
        Count the number of derivative structures
        """
        start = time()

        displacement_set = self.base_structure.frac_coords
        count = 0
        for hnf in tqdm(self.list_reduced_HNF, disable=not self.verbose):
            ds_permutation = DerivativeStructurePermutation(
                hnf, displacement_set, self.rotations, self.translations
            )
            # enumerate colorings
            count_hnf = self._count_with_hnf(hnf, ds_permutation)

            count += count_hnf

        end = time()
        if self.verbose:
            print("total: {} (Time: {:.4}sec)".format(count, end - start))

        return count

    def _count_with_hnf(
        self,
        hnf: np.ndarray,
        ds_permutation: DerivativeStructurePermutation,
    ) -> int:
        dd = Universe()

        num_sites = ds_permutation.num_sites
        automorphism = [
            Permutation(sigma) for sigma in ds_permutation.get_symmetry_operation_permutations()
        ]
        translations = [Permutation(sigma) for sigma in ds_permutation._prm_t]

        composition_constraints_dd: list[int] = []
        if self.composition_constraints is not None:
            ratio_sum = np.sum(self.composition_constraints)
            if num_sites % ratio_sum != 0:
                # impossible to satisfy composition constraints
                return 0
            else:
                composition_constraints_dd = [
                    ratio * (num_sites // ratio_sum) for ratio in self.composition_constraints
                ]

        construct_derivative_structures(
            dd,
            num_sites=num_sites,
            num_types=self.num_types,
            automorphism=automorphism,
            translations=translations,
            composition_constraints=composition_constraints_dd,
            site_constraints=self.prohibited_site_constraints,
            remove_incomplete=self.remove_incomplete,
            remove_superperiodic=self.remove_superperiodic,
        )

        return int(dd.cardinality())


def remove_symmetry_duplicates(
    base_structure,
    hnf,
    num_type,
    list_colorings,
    color_exchange: bool,
    leave_superperiodic: bool,
    use_all_colors: bool,
):
    displacement_set = base_structure.frac_coords
    rotations, translations = get_symmetry_operations(base_structure)
    cl_generator = ListBasedColoringGenerator(num_type, list_colorings)

    ds_permutation = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
    sc_enum = SiteColoringEnumerator(
        num_type, ds_permutation, cl_generator, color_exchange, leave_superperiodic, use_all_colors
    )
    colorings = sc_enum.unique_colorings()
    return colorings


def remove_symmetry_duplicates_from_generator(
    base_structure,
    hnf,
    num_type,
    list_colorings,
    color_exchange: bool,
    leave_superperiodic: bool,
    use_all_colors: bool,
    method="direct",
    n_jobs=1,
):
    displacement_set = base_structure.frac_coords
    rotations, translations = get_symmetry_operations(base_structure)
    cl_generator = ListBasedColoringGenerator(num_type, list_colorings)

    ds_permutation = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
    sc_enum = SiteColoringEnumerator(
        num_type,
        ds_permutation,
        cl_generator,
        color_exchange,
        leave_superperiodic,
        use_all_colors,
        method=method,
        n_jobs=n_jobs,
    )

    colorings = sc_enum.unique_colorings()
    return colorings
