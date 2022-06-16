import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher
from tqdm import tqdm

from dsenum.coloring import SiteColoringEnumerator
from dsenum.coloring_generator import (
    ColoringGenerator,
    FixedConcentrationColoringGenerator,
)
from dsenum.enumerate import StructureEnumerator, ZddStructureEnumerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.polya import polya_counting, polya_fixed_degrees_counting
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.utils import get_lattice

FCC_BINARY_COUNTS = [
    2,
    3,
    12,
    14,
    50,
    52,
    229,
    252,
    685,
    682,
    3875,
    2624,
    9628,
    16584,
    49764,
    42135,
    212612,
    174104,
    867893,
    1120708,
    2628180,
    3042732,
]


@pytest.mark.parametrize(
    "kind,num_types,index_start,num_expected",
    [
        ("fcc", 2, 2, FCC_BINARY_COUNTS),
        ("sc", 2, 2, [3, 3, 15]),
        ("fcc", 3, 3, [3, 13, 23, 130, 197, 1267, 2322, 9332]),
        ("fcc", 4, 4, [7, 9, 110, 211, 2110, 5471, 32362]),
        ("hcp", 2, 2, [7, 30, 163, 366, 2613, 5268, 42901, 119528, 662193]),
    ],
)
@pytest.mark.parametrize("method", ["direct", "lexicographic"])
def test_colorings_small(kind, num_types, index_start, num_expected, method):
    # Skip too time-consuming cases...
    max_index = 7

    structure = get_lattice(kind)
    index_end = index_start + len(num_expected)
    for index, expected in zip(range(index_start, index_end), num_expected):
        if index > max_index:
            continue
        se = StructureEnumerator(
            structure,
            index,
            num_types,
            color_exchange=True,
            remove_superperiodic=True,
            method=method,
        )
        actual = se.generate()
        assert len(actual) == expected


@pytest.mark.parametrize(
    "kind,num_types,composition_constraints,site_constraints",
    [
        ("fcc", 2, None, None),
        ("sc", 2, None, None),
        ("fcc", 3, None, None),
        ("fcc", 4, None, None),
        ("hcp", 2, None, None),
        ("fcc", 2, [2, 1], None),
        ("fcc", 3, [1, 2, 1], None),
        ("hcp", 2, None, [[0], [1]]),
    ],
)
def test_zdd_with_naive(kind, num_types, composition_constraints, site_constraints):
    structure = get_lattice(kind)

    for index in range(1, 6 + 1):
        if composition_constraints:
            ratio_sum = sum(composition_constraints)
            if structure.num_sites * index % ratio_sum != 0:
                # Incorrect composition ratio
                continue

        se = StructureEnumerator(
            structure,
            index,
            num_types,
            composition_constraints=composition_constraints,
            base_site_constraints=site_constraints,
            color_exchange=False,
            remove_incomplete=False,
            remove_superperiodic=False,
        )
        zse = ZddStructureEnumerator(
            structure,
            index,
            num_types,
            composition_constraints=composition_constraints,
            base_site_constraints=site_constraints,
            remove_superperiodic=False,
            remove_incomplete=False,
        )

        count_naive = len(se.generate())
        count_zdd = len(zse.generate())
        assert count_zdd == count_naive


@pytest.mark.parametrize(
    "kind,num_types,composition_constraints,site_constraints",
    [
        ("fcc", 2, None, None),
        ("sc", 2, None, None),
        ("fcc", 3, None, None),
        ("fcc", 4, None, None),
        ("hcp", 2, None, None),
        ("fcc", 2, [2, 1], None),
        ("fcc", 3, [1, 2, 1], None),
        ("hcp", 2, None, [[0], [1]]),
    ],
)
def test_enumerated_structures(kind, num_types, composition_constraints, site_constraints):
    structure = get_lattice(kind)

    for index in range(1, 5):
        # TODO: when remove_superperiodic = remove_incomplete = False, this test is failed.
        zse = ZddStructureEnumerator(
            structure,
            index,
            num_types,
            composition_constraints=composition_constraints,
            base_site_constraints=site_constraints,
            remove_superperiodic=True,
            remove_incomplete=True,
        )
        list_dstructs = zse.generate()

        stm = StructureMatcher(ltol=1e-4, stol=1e-4)
        grouped = stm.group_structures(list_dstructs)
        assert len(grouped) == len(list_dstructs)


@pytest.mark.parametrize(
    "kind,num_types",
    [
        ("fcc", 2),
        ("sc", 2),
        ("fcc", 3),
        ("fcc", 4),
        ("hcp", 2),
    ],
)
def test_colorings_with_polya(kind, num_types):
    structure = get_lattice(kind)
    displacement_set = structure.frac_coords
    num_sites_base = structure.num_sites

    for index in range(num_types, 5):
        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, structure, return_symops=True
        )
        num_sites = num_sites_base * index
        cl_generator = ColoringGenerator(num_sites, num_types)

        for hnf in tqdm(list_reduced_HNF):
            ds_permutation = DerivativeStructurePermutation(
                hnf, displacement_set, rotations, translations
            )
            sc_enum = SiteColoringEnumerator(
                num_types,
                ds_permutation,
                cl_generator,
                color_exchange=False,
                remove_superperiodic=False,
                remove_incomplete=False,
            )
            colorings = sc_enum.unique_colorings()
            cnt_polya = polya_counting(sc_enum.permutation_group, num_types)
            assert len(colorings) == cnt_polya


@pytest.mark.parametrize(
    "kind,num_types",
    [
        ("fcc", 2),
        ("sc", 2),
        ("fcc", 3),
        ("fcc", 4),
        ("hcp", 2),
    ],
)
def test_fixed_colorings_with_polya(kind, num_types):
    structure = get_lattice(kind)
    displacement_set = structure.frac_coords
    num_sites_base = structure.num_sites

    for index in range(num_types, 5):
        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, structure, return_symops=True
        )
        num_sites = num_sites_base * index
        # TODO: test at more color_ratio cases
        color_ratio = [1] * (num_types - 1) + [
            num_sites - num_types + 1,
        ]
        cl_generator = FixedConcentrationColoringGenerator(num_sites, num_types, color_ratio)

        for hnf in tqdm(list_reduced_HNF):
            ds_permutaion = DerivativeStructurePermutation(
                hnf, displacement_set, rotations, translations
            )
            sc_enum = SiteColoringEnumerator(
                num_types,
                ds_permutaion,
                cl_generator,
                color_exchange=False,
                remove_superperiodic=False,
                remove_incomplete=False,
            )
            colorings = sc_enum.unique_colorings()
            cnt_polya = polya_fixed_degrees_counting(
                sc_enum.permutation_group, num_types, color_ratio
            )
            assert len(colorings) == cnt_polya
