import pytest

from dsenum.superlattice import (
    generate_all_superlattices,
    generate_symmetry_distinct_superlattices,
)
from dsenum.utils import get_lattice


def test_generate_all_superlattices():
    # https://oeis.org/A001001
    num_expected = [
        1,
        7,
        13,
        35,
        31,
        91,
        57,
        155,
        130,
        217,
        133,
        455,
        183,
        399,
        403,
        651,
        307,
        910,
        381,
        1085,
        741,
        931,
        553,
        2015,
        806,
        1281,
        1210,
        1995,
        871,
        2821,
        993,
        2667,
        1729,
        2149,
        1767,
        4550,
        1407,
        2667,
        2379,
        4805,
        1723,
        5187,
        1893,
        4655,
        4030,
        3871,
        2257,
        8463,
        2850,
        5642,
        3991,
        6405,
        2863,
    ]
    max_index = len(num_expected)

    for index, expected in zip(range(1, max_index + 1), num_expected):
        list_HNF = generate_all_superlattices(index)
        assert len(list_HNF) == expected


@pytest.mark.parametrize(
    "kind,list_expected",
    [
        ("fcc", [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]),
        ("bcc", [1, 2, 3, 7, 5, 10, 7, 20, 14, 18]),
        ("sc", [1, 3, 3, 9, 5, 13, 7, 24, 14, 23]),
        ("hex", [1, 3, 5, 11, 7, 19, 11, 34, 23, 33]),
        ("tet", [1, 5, 5, 17, 9, 29, 13, 51, 28, 53]),
        ("hcp", [1, 3, 5, 11, 7, 19, 11, 34, 23, 33]),
    ],
)
def test_reduce_HNF_list_by_parent_lattice_symmetry(kind, list_expected):
    # Confirm table 4
    structure = get_lattice(kind)
    for index, expected in zip(range(1, len(list_expected) + 1), list_expected):
        list_reduced_HNF = generate_symmetry_distinct_superlattices(index, structure)
        assert len(list_reduced_HNF) == expected


@pytest.mark.parametrize("kind", ["fcc", "bcc"])
def test_reduce_HNF_list_by_parent_lattice_symmetry_fcc_bcc(kind):
    # https://oeis.org/A045790
    lst_num = [
        1,
        2,
        3,
        7,
        5,
        10,
        7,
        20,
        14,
        18,
        11,
        41,
        15,
        28,
        31,
        58,
        21,
        60,
        25,
        77,
        49,
        54,
        33,
        144,
        50,
        72,
        75,
        123,
        49,
        158,
        55,
        177,
        97,
        112,
        99,
        268,
        75,
        136,
        129,
        286,
        89,
        268,
        97,
        249,
        218,
        190,
        113,
        496,
        146,
        280,
    ]

    structure = get_lattice(kind)
    for index, expected in zip(range(1, len(lst_num) + 1), lst_num):
        if index >= 20:
            continue
        list_reduced_HNF = generate_symmetry_distinct_superlattices(index, structure)
        assert len(list_reduced_HNF) == expected
