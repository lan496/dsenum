import numpy as np
import pytest
from pymatgen.core import Lattice, Structure

from dsenum import StructureEnumerator, ZddStructureEnumerator
from dsenum.utils import get_lattice


@pytest.mark.parametrize("method", ["direct", "zdd", "zdd_count"])
@pytest.mark.benchmark(group="basic-binary")
def test_basic_binary(method, benchmark):
    setting = {
        "base_structure": get_lattice("fcc"),
        "index": 16,
        "num_types": 2,
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }

    if method == "direct":
        se = StructureEnumerator(**setting)
        benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)
    else:
        zse = ZddStructureEnumerator(**setting)
        if method == "zdd":
            benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)
        elif method == "zdd_count":
            benchmark.pedantic(zse.count, iterations=1)


@pytest.mark.parametrize("method", ["direct", "zdd", "zdd_count"])
@pytest.mark.benchmark(group="basic-ternary")
def test_basic_ternary(method, benchmark):
    setting = {
        "base_structure": get_lattice("fcc"),
        "index": 10,
        "num_types": 3,
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }

    if method == "direct":
        se = StructureEnumerator(**setting)
        benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)
    else:
        zse = ZddStructureEnumerator(**setting)
        if method == "zdd":
            benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)
        elif method == "zdd_count":
            benchmark.pedantic(zse.count, iterations=1)


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


@pytest.mark.parametrize("method", ["direct", "zdd", "zdd_count"])
@pytest.mark.benchmark(group="composition")
def test_composition_constraints(method, benchmark):
    setting = {
        "base_structure": get_rutile_structure(),
        "index": 2,
        "num_types": 3,
        "composition_constraints": [1, 3, 2],
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }

    if method == "direct":
        se = StructureEnumerator(**setting)
        benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)
    else:
        zse = ZddStructureEnumerator(**setting)
        if method == "zdd":
            benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)
        elif method == "zdd_count":
            benchmark.pedantic(zse.count, iterations=1)


@pytest.mark.parametrize("method", ["direct", "zdd", "zdd_count"])
@pytest.mark.benchmark(group="site")
def test_site_constraints(method, benchmark):
    setting = {
        "base_structure": get_rutile_structure(),
        "index": 4,
        "num_types": 3,
        "base_site_constraints": [
            [2],  # 2a
            [2],  # 2a
            [0, 1],  # 4f
            [0, 1],  # 4f
            [0, 1],  # 4f
            [0, 1],  # 4f
        ],
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }

    if method == "direct":
        se = StructureEnumerator(**setting)
        benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)
    else:
        zse = ZddStructureEnumerator(**setting)
        if method == "zdd":
            benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)
        elif method == "zdd_count":
            benchmark.pedantic(zse.count, iterations=1)


@pytest.mark.parametrize("method", ["direct", "zdd", "zdd_count"])
@pytest.mark.benchmark(group="composition-site")
def test_composition_site_constraints(method, benchmark):
    setting = {
        "base_structure": get_rutile_structure(),
        "index": 4,
        "num_types": 3,
        "composition_constraints": [1, 3, 2],
        "base_site_constraints": [
            [2],  # 2a
            [2],  # 2a
            [0, 1],  # 4f
            [0, 1],  # 4f
            [0, 1],  # 4f
            [0, 1],  # 4f
        ],
        "remove_incomplete": True,
        "remove_superperiodic": True,
    }

    if method == "direct":
        se = StructureEnumerator(**setting)
        benchmark.pedantic(se.generate, kwargs={"output": "poscar"}, iterations=1)
    else:
        zse = ZddStructureEnumerator(**setting)
        if method == "zdd":
            benchmark.pedantic(zse.generate, kwargs={"output": "poscar"}, iterations=1)
        elif method == "zdd_count":
            benchmark.pedantic(zse.count, iterations=1)
