import numpy as np
import pytest
from pymatgen.core import Structure, Lattice

from dsenum import StructureEnumerator, ZddStructureEnumerator
from dsenum.utils import get_lattice


def get_fcc_structure():
    a = 2.856

    lattice = Lattice.from_parameters(a, a, a, 60, 60, 60)
    species = ["Cu"]
    frac_coords = np.array([[0, 0, 0]])
    structure = Structure(lattice, species, frac_coords)
    return structure


def get_rutile_structure():
    # rutile structure taken from mp-856
    a = 4.832
    c = 3.243
    x_4f = 0.3066

    lattice = Lattice.from_parameters(a, a, c, 90, 90, 90)
    species = ["Sn", "Sn", "O", "O", "O", "O"]
    # fmt: off
    frac_coords = np.array([
        [0, 0, 0],                      # Sn(2a)
        [0.5, 0.5, 0.5],                # Sn(2a)
        [x_4f, x_4f, 0],                # O(4f)
        [1 - x_4f, 1 - x_4f, 0],        # O(4f)
        [0.5 - x_4f, 0.5 + x_4f, 0.5],  # O(4f)
        [0.5 + x_4f, 0.5 - x_4f, 0.5],  # O(4f)
    ])
    # fmt: on
    structure = Structure(lattice, species, frac_coords)
    return structure


def get_common_settings():
    data = {
        # "base_structure": get_rutile_structure(),
        "base_structure": get_fcc_structure(),
        "index": 16,
        "num_types": 2,
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }
    return data


@pytest.mark.benchmark(group="basic")
def test_direct(benchmark):
    setting = get_common_settings()
    setting["method"] = "direct"
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)


@pytest.mark.benchmark(group="basic")
def test_zdd(benchmark):
    setting = get_common_settings()
    zse = ZddStructureEnumerator(**setting)

    benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)


@pytest.mark.benchmark(group="basic")
def test_zdd_counting(benchmark):
    setting = get_common_settings()
    zse = ZddStructureEnumerator(**setting)

    benchmark.pedantic(zse.count, iterations=1)
