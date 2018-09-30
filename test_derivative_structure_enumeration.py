import unittest

from derivative_structure_enumeration import (
    generate_all_superlattices,
)


class TestDerivativeStructureEnumeration(unittest.TestCase):

    def test_generate_all_superlattices(self):
        # https://oeis.org/A001001
        num_expected = [1, 7, 13, 35, 31, 91, 57, 155, 130, 217,
                        133, 455, 183, 399, 403, 651, 307, 910, 381, 1085,
                        741, 931, 553, 2015, 806, 1281, 1210, 1995, 871, 2821,
                        993, 2667, 1729, 2149, 1767, 4550, 1407, 2667, 2379, 4805,
                        1723, 5187, 1893, 4655, 4030, 3871, 2257, 8463, 2850, 5642,
                        3991, 6405, 2863]
        max_index = len(num_expected)

        for index, expected in zip(range(1, max_index + 1), num_expected):
            list_HNF = generate_all_superlattices(index)
            self.assertEqual(len(list_HNF), expected)


if __name__ == '__main__':
    unittest.main()
