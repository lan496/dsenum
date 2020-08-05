from tqdm import tqdm

from dsenum.enumerate import StructureEnumerator
from dsenum.coloring_generator import ColoringGenerator, FixedConcentrationColoringGenerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.utils import get_lattice
from dsenum.polya import polya_counting, polya_fixed_degrees_counting
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.coloring import SiteColoringEnumerator


obj = {
    "fcc": {
        "structure": get_lattice("fcc"),
        "num_type": 2,
        "indices": range(2, 23 + 1),
        "num_expected": [
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
        ],
    },
    "sc": {
        "structure": get_lattice("sc"),
        "num_type": 2,
        "indices": range(2, 4 + 1),
        "num_expected": [3, 3, 15],
    },
    "fcc_ternary": {
        "structure": get_lattice("fcc"),
        "num_type": 3,
        "indices": range(3, 10 + 1),
        "num_expected": [3, 13, 23, 130, 197, 1267, 2322, 9332],
    },
    "fcc_quaternary": {
        "structure": get_lattice("fcc"),
        "num_type": 4,
        "indices": range(4, 10 + 1),
        "num_expected": [7, 9, 110, 211, 2110, 5471, 32362],
    },
    "hcp": {
        "structure": get_lattice("hcp"),
        "num_type": 2,
        "indices": range(2, 10 + 1),
        "num_expected": [7, 30, 163, 366, 2613, 5268, 42901, 119528, 662193],
    },
}


def test_colorings():
    for name, dct in obj.items():
        structure = dct["structure"]
        # displacement_set = structure.frac_coords
        num_type = dct["num_type"]
        for index, expected in zip(dct["indices"], dct["num_expected"]):
            if index >= 8:
                continue
            for method in ["direct", "lexicographic"]:
                se = StructureEnumerator(
                    structure,
                    index,
                    num_type,
                    color_exchange=True,
                    leave_superperiodic=False,
                    method=method,
                )
                actual = se.generate()
                assert len(actual) == expected


def test_colorings_with_polya():
    for name, dct in obj.items():
        structure = dct["structure"]
        displacement_set = structure.frac_coords
        num_type = dct["num_type"]
        num_sites_base = structure.num_sites

        for index, expected in zip(dct["indices"], dct["num_expected"]):
            if index >= 5:
                continue

            list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
                index, structure, return_symops=True
            )
            num_sites = num_sites_base * index
            cl_generator = ColoringGenerator(num_sites, num_type)

            for hnf in tqdm(list_reduced_HNF):
                ds_permutation = DerivativeStructurePermutation(
                    hnf, displacement_set, rotations, translations
                )
                sc_enum = SiteColoringEnumerator(
                    num_type,
                    ds_permutation,
                    cl_generator,
                    color_exchange=False,
                    leave_superperiodic=True,
                    use_all_colors=False,
                )
                colorings = sc_enum.unique_colorings()
                cnt_polya = polya_counting(sc_enum.permutation_group, num_type)
                assert len(colorings) == cnt_polya


def test_fixed_colorings_with_polya():
    for name, dct in obj.items():
        structure = dct["structure"]
        displacement_set = structure.frac_coords
        num_type = dct["num_type"]
        num_sites_base = structure.num_sites

        for index, expected in zip(dct["indices"], dct["num_expected"]):
            if index >= 5:
                continue

            list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
                index, structure, return_symops=True
            )
            num_sites = num_sites_base * index
            # TODO: test at more color_ratio cases
            color_ratio = [1] * (num_type - 1) + [
                num_sites - num_type + 1,
            ]
            cl_generator = FixedConcentrationColoringGenerator(num_sites, num_type, color_ratio)

            for hnf in tqdm(list_reduced_HNF):
                ds_permutaion = DerivativeStructurePermutation(
                    hnf, displacement_set, rotations, translations
                )
                sc_enum = SiteColoringEnumerator(
                    num_type,
                    ds_permutaion,
                    cl_generator,
                    color_exchange=False,
                    leave_superperiodic=True,
                    use_all_colors=False,
                )
                colorings = sc_enum.unique_colorings()
                cnt_polya = polya_fixed_degrees_counting(
                    sc_enum.permutation_group, num_type, color_ratio
                )
                assert len(colorings) == cnt_polya
