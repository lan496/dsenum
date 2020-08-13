import pytest

from dsenum.coloring_generator import ColoringGenerator
from dsenum.utils import get_lattice


def test_coloring_generator():
    num_elements = 8
    num_color = 3
    cg = ColoringGenerator(num_elements, num_color)

    list_colorings, flags = cg.generate_all_colorings()
    assert len(set(list_colorings)) == num_color ** num_elements
    assert len(set(flags.keys())) == num_color ** num_elements

    list_colorings2 = list(cg.yield_coloring())
    assert len(set([tuple(e) for e in list_colorings2])) == num_color ** num_elements


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

    colorings_expect = set(
        [
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
        ]
    )
    list_colorings, flags = cg.generate_all_colorings()
    assert set([tuple(e) for e in list_colorings]) == colorings_expect
    assert len(set(flags.keys())) == len(colorings_expect)

    list_colorings2 = list(cg.yield_coloring())
    assert set([tuple(e) for e in list_colorings2]) == colorings_expect


def test_composition_constraints():
    pass


def test_composition_site_constraints():
    pass
