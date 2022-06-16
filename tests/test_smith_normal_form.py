from hsnf import smith_normal_form

from dsenum.superlattice import generate_all_superlattices


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
