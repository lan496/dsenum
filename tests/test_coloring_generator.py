import pytest

from dsenum.coloring_generator import (
    ColoringGenerator,
    FixedConcentrationColoringGenerator,
    satisfy_site_constraints,
)
from dsenum.utils import get_lattice


def test_coloring_generator():
    num_elements = 8
    num_color = 3
    cg = ColoringGenerator(num_elements, num_color)

    list_colorings, flags = cg.generate_all_colorings()
    assert len(set(list_colorings)) == num_color**num_elements
    assert len(set(flags.keys())) == num_color**num_elements

    list_colorings2 = list(cg.yield_coloring())
    assert len({tuple(e) for e in list_colorings2}) == num_color**num_elements


def test_coloring_generator_with_site_constraints():
    num_color = 3
    site_constraints = [
        [0],
        [1],
        [2],
        [1, 2],
        [0, 2],
        [0, 1, 2],
    ]
    num_elements = len(site_constraints)
    cg = ColoringGenerator(num_elements, num_color, site_constraints)

    colorings_expect = {
        (0, 1, 2, 1, 0, 0),
        (0, 1, 2, 1, 0, 1),
        (0, 1, 2, 1, 0, 2),
        (0, 1, 2, 1, 2, 0),
        (0, 1, 2, 1, 2, 1),
        (0, 1, 2, 1, 2, 2),
        (0, 1, 2, 2, 0, 0),
        (0, 1, 2, 2, 0, 1),
        (0, 1, 2, 2, 0, 2),
        (0, 1, 2, 2, 2, 0),
        (0, 1, 2, 2, 2, 1),
        (0, 1, 2, 2, 2, 2),
    }
    list_colorings, flags = cg.generate_all_colorings()
    assert {tuple(e) for e in list_colorings} == colorings_expect
    assert len(set(flags.keys())) == len(colorings_expect)

    list_colorings2 = list(cg.yield_coloring())
    assert {tuple(e) for e in list_colorings2} == colorings_expect


def test_fixed_concentration_coloring_generator():
    num_elements = 4
    num_color = 3
    color_ratio = [1, 1, 2]
    fcg = FixedConcentrationColoringGenerator(num_elements, num_color, color_ratio)

    colorings_expect = {
        (0, 1, 2, 2),
        (0, 2, 1, 2),
        (0, 2, 2, 1),
        (1, 0, 2, 2),
        (1, 2, 0, 2),
        (1, 2, 2, 0),
        (2, 0, 1, 2),
        (2, 0, 2, 1),
        (2, 1, 0, 2),
        (2, 1, 2, 0),
        (2, 2, 0, 1),
        (2, 2, 1, 0),
    }

    list_colorings, flags = fcg.generate_all_colorings()
    assert {tuple(e) for e in list_colorings} == colorings_expect
    assert len(set(flags.keys())) == len(colorings_expect)

    list_colorings2 = list(fcg.yield_coloring())
    assert {tuple(e) for e in list_colorings2} == colorings_expect


def test_fixed_concentration_coloring_generator_with_site_constraints():
    num_elements = 4
    num_color = 3
    color_ratio = [1, 1, 2]
    site_constraints = [[2], [1, 2], [0, 1, 2], [0, 1, 2]]
    fcg = FixedConcentrationColoringGenerator(
        num_elements, num_color, color_ratio, site_constraints
    )

    colorings_expect = {(2, 1, 0, 2), (2, 1, 2, 0), (2, 2, 0, 1), (2, 2, 1, 0)}

    list_colorings, flags = fcg.generate_all_colorings()
    assert {tuple(e) for e in list_colorings} == colorings_expect
    assert len(set(flags.keys())) == len(colorings_expect)

    list_colorings2 = list(fcg.yield_coloring())
    assert {tuple(e) for e in list_colorings2} == colorings_expect


def test_satisfy_site_constraints():
    site_constraints = [
        [0],
        [1],
        [0, 1],
        [0, 1, 2],
    ]
    coloring1 = [0, 1, 1, 0]
    coloring2 = [1, 1, 1, 1]
    assert satisfy_site_constraints(site_constraints, coloring1)
    assert not satisfy_site_constraints(site_constraints, coloring2)
