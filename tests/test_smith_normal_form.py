import numpy as np

from dsenum.smith_normal_form import smith_normal_form
from dsenum.superlattice import generate_all_superlattices


def test_smf():
    list_matrix = [
        np.array([[2, 0], [1, 4]]),
        np.array([[2, 4, 4], [-6, 6, 12], [10, -4, -16]]),
        np.array([[8, 4, 8], [4, 8, 4]]),
        np.array([[-6, 111, -36, 6], [5, -672, 210, 74], [0, -255, 81, 24], [-7, 255, -81, -10]]),
        np.array([[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]),
        np.array([[1, 0, 0], [1, 2, 0], [0, 0, 2]]),
    ]
    list_expected = [
        np.diag([1, 8]),
        np.diag([2, 6, 12]),
        np.array([[4, 0, 0], [0, 12, 0]]),
        np.diag([1, 3, 21, 0]),
        np.diag([1, 4, 4]),
        np.diag([1, 2, 2]),
    ]

    for M, expected in zip(list_matrix, list_expected):
        D, L, R = smith_normal_form(M)
        D_re = np.dot(L, np.dot(M, R))
        assert np.isclose(np.linalg.det(L) ** 2, 1)
        assert np.isclose(np.linalg.det(R) ** 2, 1)
        assert np.array_equal(D_re, D)


def test_number_of_snf():
    # confirm table-3
    num_hnf_expected = [1, 7, 13, 35, 31, 91, 57, 155, 130, 217, 133, 455, 183, 399, 403, 651]
    num_snf_expected = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 4]
    max_index = len(num_hnf_expected)

    for index, hnf_expected, snf_expected in zip(
        range(1, max_index + 1), num_hnf_expected, num_snf_expected
    ):
        list_HNF = generate_all_superlattices(index)
        assert len(list_HNF) == hnf_expected

        list_SNF = set()
        for hnf in list_HNF:
            snf, _, _ = smith_normal_form(hnf)
            dag = tuple(snf.diagonal())
            list_SNF.add(dag)

        assert len(list_SNF) == snf_expected
