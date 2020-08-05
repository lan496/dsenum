import pytest

from dsenum.polya import polya_counting, polya_fixed_degrees_counting


@pytest.fixture
def dihedral4():
    Dihedral4 = [
        [0, 1, 2, 3],  # E
        [3, 0, 1, 2],  # C4
        [2, 3, 0, 1],
        [1, 2, 3, 0],
        [1, 0, 3, 2],  # mirror
        [2, 1, 0, 3],
        [3, 2, 1, 0],
        [0, 3, 2, 1],
    ]
    return Dihedral4


def test_polya_counting(dihedral4):
    num_color = 2
    expected = 6
    actual = polya_counting(dihedral4, num_color)
    assert actual == expected


def test_fixed_degrees(dihedral4):
    num_color = 2

    num_element_of_each_color = [4, 0]
    expected = 1
    actual = polya_fixed_degrees_counting(dihedral4, num_color, num_element_of_each_color)
    assert actual == expected

    num_element_of_each_color = [3, 1]
    expected = 1
    actual = polya_fixed_degrees_counting(dihedral4, num_color, num_element_of_each_color)
    assert actual == expected

    num_element_of_each_color = [2, 2]
    expected = 2
    actual = polya_fixed_degrees_counting(dihedral4, num_color, num_element_of_each_color)
    assert actual == expected

    num_element_of_each_color = [1, 3]
    expected = 1
    actual = polya_fixed_degrees_counting(dihedral4, num_color, num_element_of_each_color)
    assert actual == expected

    num_element_of_each_color = [0, 4]
    expected = 1
    actual = polya_fixed_degrees_counting(dihedral4, num_color, num_element_of_each_color)
    assert actual == expected
