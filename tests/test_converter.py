import numpy as np

from dsenum.converter import DerivativeMultiLatticeHash


def test_converter():
    list_transformations = [
        np.diag((2, 2, 2)),
        np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
        np.array([[-1, 1, 1], [1, -1, 1], [-1, -1, 1]]),  # negative determinat
    ]
    list_num_types = [2, 3]
    list_frac_coords = [np.array([[0, 0, 0]]), np.array([[0, 0, 0], [0.5, 0.5, 0.5]])]

    for transformation in list_transformations:
        for num_types in list_num_types:
            for frac_coords in list_frac_coords:
                converter = DerivativeMultiLatticeHash(transformation, frac_coords)

                index = converter.index
                num_sites = converter.num_sites

                all_factors = converter.get_all_factors()
                assert len(set([tuple((f)) for f in all_factors])) == index

                all_periodic_sites = converter.get_canonical_sites_list()
                assert len(set(all_periodic_sites)) == num_sites

                # TODO: test converter.hash_to_periodic
                # TODO: test converter.augment_permutation

                for csite in all_periodic_sites:
                    ind = converter.ravel_canonical_site(csite)
                    csite2 = converter.unravel_to_canonical_site(ind)
                    assert csite2 == csite
