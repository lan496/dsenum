import numpy as np

from dsenum.converter import DerivativeMultiLatticeHash


def test_converter():
    list_transformations = [
        np.diag((2, 2, 2)),
        np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]),
        np.array([[-1, 1, 1], [1, -1, 1], [-1, -1, 1]]),  # negative determinat
    ]
    list_frac_coords = [np.array([[0, 0, 0]]), np.array([[0, 0, 0], [0.5, 0.5, 0.5]])]

    for transformation in list_transformations:
        for frac_coords in list_frac_coords:
            converter = DerivativeMultiLatticeHash(transformation, frac_coords)

            index = converter.index
            num_sites = converter.num_sites

            all_factors = converter.get_all_factors()
            assert len({tuple(f) for f in all_factors}) == index

            all_lattice_points = converter.get_lattice_points()
            assert len({tuple(f) for f in all_lattice_points}) == index

            all_periodic_sites = converter.get_canonical_sites_list()
            assert len(set(all_periodic_sites)) == num_sites

            all_dsites = converter.get_distinct_derivative_sites_list()
            for dsite in all_dsites:
                csite, derivative_jimage = converter.hash_derivative_site(dsite, return_image=True)
                fcoords1 = converter.get_frac_coords(dsite)
                fcoords2 = (
                    converter.displacement_set[dsite.site_index]
                    + np.dot(converter.left_inv, csite.factor)
                    + np.dot(converter.hnf, derivative_jimage)
                )
                assert np.allclose(fcoords1, fcoords2)

            for csite in all_periodic_sites:
                ind = converter.ravel_canonical_site(csite)
                csite2 = converter.unravel_to_canonical_site(ind)
                assert csite2 == csite
