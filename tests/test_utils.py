import numpy as np

from dsenum.utils import gcd_list, cast_integer_matrix


def test_gcd_list():
    assert gcd_list([2]) == 2
    assert gcd_list([2, 4]) == 2
    assert gcd_list([2, 3]) == 1
    assert gcd_list([2, 4, 6]) == 2
    assert gcd_list([2, 3, 4]) == 1


def test_cast_integer_matrix():
    arr = np.array([1.0, 2.0001, 2.9999, -0.9999, -2.0001])
    arr_rounded = np.array([1.0, 2.0, 3.0, -1.0, -2.0])
    assert np.allclose(cast_integer_matrix(arr), arr_rounded)
