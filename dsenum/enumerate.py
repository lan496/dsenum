from abc import ABCMeta, abstractmethod
from time import time
from typing import List, Tuple, Union, cast
from warnings import warn

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie, Element, Specie
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

    Parameters
    ----------
    base_structure: pymatgen.core.Structure
        Aristotype for derivative structures
    index: int
        How many times to expand unit cell
    num_types: int
        The number of species in derivative structures.
        `num_types` may be larger than the number of the kinds of species in `base_structure`: for example, you consider vacancies in derivative structures.
    mapping_color_species: optional
        If specified, use these species in derivative structures.
        The length of this list should be equal to `num_types`
    composition_constraints: (Optional) List[int]
        composition_constraints[i] is the ratio of the i-th species in mapping_color_species.
        For example, set `composition_constraints = [1, 2]` for enumerating TiO2 structures with
        `mapping_color_species = ["Ti", "O"]`.
    base_site_constraints: (Optional) List[List[int]], (num_elements, num_color)
        e.g. site_constraints[2] = [0, 3, 4] means color of site-2 in base_structure must be 0, 3, or 4.
    color_exchange: (Optional) bool
        identify color-exchanging
    remove_superperiodic: (Optional) bool
        iff true, discard superperiodic coloring
    remove_incomplete: (Optional) bool

    Arguments
    ---------
    list_reduced_HNF
    rotations
    translations
    site_constraints
    mapping_color_species
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        mapping_color_species: List[Union[str, Element, Specie, DummySpecie]] = None,
        composition_constraints=None,
        base_site_constraints=None,
        color_exchange=True,
        remove_superperiodic=True,
        remove_incomplete=True,
    ):
        self.base_structure = base_structure
        self.index = index
        self.num_types = num_types

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

        # site constraints
        self.site_constraints = None
        if base_site_constraints:
            assert len(base_site_constraints) == self.num_sites_base
            self.site_constraints = convert_site_constraints(base_site_constraints, self.index)

        if mapping_color_species and len(mapping_color_species) != self.num_types:
            raise ValueError("mapping_color_species must have num_type species.")
        if mapping_color_species is None:
            mapping_color_species = [DummySpecie(str(i)) for i in range(1, self.num_types + 1)]
        self.mapping_color_species = mapping_color_species

    @property
    def num_sites_base(self):
        return self.base_structure.num_sites

    @property
    def num_sites(self):
        return self.num_sites_base * self.index

    def generate(
        self,
        return_colorings=False,
        additional_species=None,
        additional_frac_coords=None,
        output="pymatgen",
    ) -> Union[
        List[Union[Structure, str]],
        Tuple[List[Union[Structure, str]], List[np.ndarray], List[List[int]]],
    ]:
        """
        Parameters
        ----------
        return_colorings: bool, optional
            if true, return transformation matrices and colorings in addition
        additional_species: list of pymatgen.core.Species, optional
            species which are nothing to do with ordering
        additional_frac_coords: np.ndarray, optional
            fractional coordinates of species which are nothing to do with ordering

        Returns
        -------
        list_ds: list of derivative structure
        list_transformations: list of transformation matrices, optional
        list_colorings: list of colorings, optional
        """
        assert (output == "pymatgen") or (output == "poscar")
        start = time()

        displacement_set = self.base_structure.frac_coords
        list_ds = []
        list_transformations = []
        list_colorings = []
        for hnf in tqdm(self.list_reduced_HNF):
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
    ) -> List[List[int]]:
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
        `num_types` may be larger than the number of the kinds of species in `base_structure`: for example, you consider vacancies in derivative structures.
    mapping_color_species: optional
        If specified, use these species in derivative structures.
        The length of this list should be equal to `num_types`
    composition_constraints: (Optional) List[int]
        composition_constraints[i] is the ratio of the i-th species in mapping_color_species
    base_site_constraints: (Optional) List[List[int]], (num_elements, num_color)
        e.g. site_constraints[2] = [0, 3, 4] means color of site-2 in base_structure must be 0, 3, or 4.
    color_exchange: (Optional) bool
        identify color-exchanging
    remove_superperiodic: (Optional) bool
        iff true, discard superperiodic coloring
    remove_incomplete: (Optional) bool
    method: (Optional) str
        "direct" or "lexicographic", so far
    n_jobs: (Optional) int
        core in lexicographic coset enumeration(only used when method='lexicographic')

    Arguments
    ---------
    list_reduced_HNF
    rotations
    translations
    site_constraints
    mapping_color_species
    cl_generator
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        mapping_color_species: List[Union[str, Element, Specie, DummySpecie]] = None,
        composition_constraints=None,
        base_site_constraints=None,
        color_exchange=True,
        remove_superperiodic=True,
        remove_incomplete=True,
        method="direct",
        n_jobs=1,
    ):
        super().__init__(
            base_structure,
            index,
            num_types,
            mapping_color_species,
            composition_constraints,
            base_site_constraints,
            color_exchange,
            remove_superperiodic,
            remove_incomplete,
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
        if composition_constraints is None:
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
                    composition_constraints,
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
    ) -> List[List[int]]:
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


def enumerate_derivative_structures(
    base_structure,
    index,
    num_type,
    mapping_color_species=None,
    composition_constraints=None,
    base_site_constraints=None,
    color_exchange=True,
    leave_superperiodic=False,
    use_all_colors=True,
    method="direct",
    n_jobs=1,
):
    """
    Parameters
    ----------
    base_structure: Structure
    index: int
    num_type: int
    mapping_color_species: (Optional) if specified, use these species in derivative structures
    composition_constraints: (Optional) None or list of int
    base_site_constraints: (Optional) list (num_elements, num_color)
        e.g. site_constraints[2] = [0, 3, 4] means color of site-2 in base_structure must be 0, 3, or 4.
    color_exchange: identify color-exchanging
    leave_superperiodic: do not discard superperiodic coloring
    use_all_colors: bool
    method: "direct" or "lexicographic", so far
    n_jobs: core in lexicographic coset enumeration(only used when method='lexicographic')

    Returns
    -------
    list_ds: list of derivative structure
    """
    warn("Deprecated. Use StructureEnumerator instead", DeprecationWarning)
    se = StructureEnumerator(
        base_structure,
        index,
        num_type,
        mapping_color_species,
        composition_constraints,
        base_site_constraints,
        color_exchange,
        leave_superperiodic,
        use_all_colors,
        method,
        n_jobs,
    )
    return se.generate()


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
