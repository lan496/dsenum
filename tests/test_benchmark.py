import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from dsenum import StructureEnumerator
from dsenum.utils import get_lattice


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
        "base_structure": get_rutile_structure(),
        "index": 2,
        "num_types": 3,
    }
    return data


@pytest.mark.benchmark(group="basic")
def test_direct(benchmark):
    setting = get_common_settings()
    setting["method"] = "direct"
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="basic")
def test_lexicographic(benchmark):
    setting = get_common_settings()
    setting["method"] = "lexicographic"
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="composition")
def test_direct_with_composition(benchmark):
    setting = get_common_settings()
    setting["method"] = "direct"
    setting["composition_constraints"] = [1, 3, 2]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="composition")
def test_lexicographic_with_composition(benchmark):
    setting = get_common_settings()
    setting["method"] = "lexicographic"
    setting["composition_constraints"] = [1, 3, 2]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="site")
def test_direct_with_site(benchmark):
    setting = get_common_settings()
    setting["method"] = "direct"
    setting["base_site_constraints"] = [
        [2],  # 2a
        [2],  # 2a
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
    ]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="site")
def test_lexicographic_with_site(benchmark):
    setting = get_common_settings()
    setting["method"] = "lexicographic"
    setting["base_site_constraints"] = [
        [2],  # 2a
        [2],  # 2a
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
    ]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="composition-site")
def test_direct_with_composition_and_site(benchmark):
    setting = get_common_settings()
    setting["method"] = "direct"
    setting["composition_constraints"] = [1, 3, 2]
    setting["base_site_constraints"] = [
        [2],  # 2a
        [2],  # 2a
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
    ]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)


@pytest.mark.benchmark(group="composition-site")
def test_lexicographic_with_composition_and_site(benchmark):
    setting = get_common_settings()
    setting["method"] = "lexicographic"
    setting["composition_constraints"] = [1, 3, 2]
    setting["base_site_constraints"] = [
        [2],  # 2a
        [2],  # 2a
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
        [0, 1],  # 4f
    ]
    se = StructureEnumerator(**setting)

    benchmark.pedantic(se.generate, iterations=1)
