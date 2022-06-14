import sys

import numpy as np
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
    "fcc_fixed-composition": {
        "structure": get_lattice("fcc"),
        "num_types": 2,
        "composition_constraints": [2, 1],
    },
    "fcc_ternary_fixed-composition": {
        "structure": get_lattice("fcc"),
        "num_types": 3,
        "composition_constraints": [1, 2, 1],
    },
    "hcp_site-constraint": {
        "structure": get_lattice("hcp"),
        "num_types": 2,
        "site_constraints": [
            [0],
            [1],
        ],
    },
}


@pytest.mark.skipif("pyzdd" not in sys.modules, reason="requires pyzdd")
def test_with_naive_method():
    max_index = 6
    for name, dct in obj.items():
        print(name)
        structure = dct["structure"]
        num_types = dct["num_types"]
        composition_constraints = dct.get("composition_constraints", None)
        base_site_constraints = dct.get("site_constraints", None)
        for index in range(1, max_index + 1):
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
                base_site_constraints=base_site_constraints,
                color_exchange=False,
                remove_incomplete=False,
                remove_superperiodic=False,
            )
            zse = ZddStructureEnumerator(
                structure,
                index,
                num_types,
                composition_constraints=composition_constraints,
                base_site_constraints=base_site_constraints,
                remove_superperiodic=False,
                remove_incomplete=False,
            )

            count_naive = len(se.generate())
            count_zdd = len(zse.generate())
            print(count_naive, count_zdd)
            assert count_zdd == count_naive


def test_enumerated_structures():
    max_index = 4

    for name, dct in obj.items():
        print("strcuture check", name)
        structure = dct["structure"]
        num_types = dct["num_types"]
        composition_constraints = dct.get("composition_constraints", None)
        base_site_constraints = dct.get("site_constraints", None)
        for index in range(1, max_index + 1):
            # TODO: when remove_superperiodic = remove_incomplete = False, this test is failed.
            zse = ZddStructureEnumerator(
                structure,
                index,
                num_types,
                composition_constraints=composition_constraints,
                base_site_constraints=base_site_constraints,
                remove_superperiodic=True,
                remove_incomplete=True,
            )
            list_dstructs = zse.generate()

            stm = StructureMatcher(ltol=1e-4, stol=1e-4)
            grouped = stm.group_structures(list_dstructs)
            assert len(grouped) == len(list_dstructs)
