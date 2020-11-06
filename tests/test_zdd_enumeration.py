from tqdm import tqdm
import pytest
from pymatgen.analysis.structure_matcher import StructureMatcher

from dsenum.enumerate import StructureEnumerator, ZddStructureEnumerator
from dsenum.coloring_generator import ColoringGenerator, FixedConcentrationColoringGenerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.utils import get_lattice
from dsenum.polya import polya_counting, polya_fixed_degrees_counting
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.coloring import SiteColoringEnumerator


obj = {
    "fcc": {
        "structure": get_lattice("fcc"),
        "num_types": 2,
    },
    "sc": {
        "structure": get_lattice("sc"),
        "num_types": 2,
    },
    "fcc_ternary": {
        "structure": get_lattice("fcc"),
        "num_types": 3,
    },
    "fcc_quaternary": {
        "structure": get_lattice("fcc"),
        "num_types": 4,
    },
    "hcp": {
        "structure": get_lattice("hcp"),
        "num_types": 2,
    },
}


def test_with_naive_method():
    max_index = 6
    for name, dct in obj.items():
        structure = dct["structure"]
        num_types = dct["num_types"]
        for index in range(1, max_index + 1):
            se = StructureEnumerator(
                structure,
                index,
                num_types,
                color_exchange=False,
                remove_incomplete=True,
                remove_superperiodic=True,
            )
            zse = ZddStructureEnumerator(
                structure,
                index,
                num_types,
                remove_superperiodic=True,
                remove_incomplete=True,
            )

            count_naive = len(se.generate())
            count_zdd = len(zse.generate())
            assert count_zdd == count_naive


def test_enumerated_structures():
    max_index = 4

    for name, dct in obj.items():
        structure = dct["structure"]
        num_types = dct["num_types"]
        for index in range(1, max_index + 1):
            zse = ZddStructureEnumerator(
                structure,
                index,
                num_types,
                remove_superperiodic=True,
                remove_incomplete=True,
            )
            list_dstructs = zse.generate()

            stm = StructureMatcher(ltol=1e-4, stol=1e-4)
            grouped = stm.group_structures(list_dstructs)
            assert len(grouped) == len(list_dstructs)


@pytest.mark.skip
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
                    remove_superperiodic=False,
                    remove_incomplete=False,
                )
                colorings = sc_enum.unique_colorings()
                cnt_polya = polya_counting(sc_enum.permutation_group, num_type)
                assert len(colorings) == cnt_polya


@pytest.mark.skip
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
                    remove_superperiodic=False,
                    remove_incomplete=False,
                )
                colorings = sc_enum.unique_colorings()
                cnt_polya = polya_fixed_degrees_counting(
                    sc_enum.permutation_group, num_type, color_ratio
                )
                assert len(colorings) == cnt_polya
