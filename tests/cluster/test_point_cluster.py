import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie

from dsenum.cluster.point_cluster import EquivalentPointClusterGenerator, PointCluster
from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.site import DerivativeSite
from dsenum.utils import get_symmetry_operations, square2d_lattice_symmetry


def test_point_cluster():
    pcl0 = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))])
    pcl1 = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))])
    pcl2 = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, 1))])
    assert pcl0 == pcl1
    assert pcl0 != pcl2
    assert len(pcl0) == 2

    # string representation
    pcl0_str = "PointCluster[site_index=0, jimage=(0, 0, 0);site_index=0, jimage=(0, 0, 1)]"
    assert str(pcl0) == pcl0_str


@pytest.fixture
def epcg_fcc():
    fcc = Structure(Lattice([[0, 1, 1], [1, 0, 1], [1, 1, 0]]), [DummySpecie("X")], [[0, 0, 0]])
    transformation = np.eye(3, dtype=int)

    epcg = EquivalentPointClusterGenerator.from_structure(fcc, transformation)
    return epcg


@pytest.fixture
def epcg_bcc():
    bcc = Structure(Lattice([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]), [DummySpecie("X")], [[0, 0, 0]])
    transformation = np.eye(3, dtype=int)

    epcg = EquivalentPointClusterGenerator.from_structure(bcc, transformation)
    return epcg


@pytest.fixture
def epcg_hcp():
    # TODO
    raise NotImplementedError


@pytest.fixture
def epcg_perovskite_oxygen_substructure():
    lattice = Lattice(np.eye(3))
    # fmt: off
    perovskite = Structure(lattice,
                           ['Cu', 'La', 'O', 'O', 'O'],
                           [
                               [0, 0, 0],        # B=Cu
                               [0.5, 0.5, 0.5],  # A=La
                               [0.5, 0, 0],      # X=O
                               [0, 0.5, 0],      # X=O
                               [0, 0, 0.5],      # X=O
                           ])
    oxygen_substructure = Structure(lattice,
                                    ['O', 'O', 'O'],
                                    [
                                        [0.5, 0, 0],      # X=O
                                        [0, 0.5, 0],      # X=O
                                        [0, 0, 0.5],      # X=O
                                    ])
    # fmt: on

    transformation = np.eye(3, dtype=int)

    frac_coords = oxygen_substructure.frac_coords
    rotations, translations = get_symmetry_operations(perovskite)
    converter = DerivativeMultiLatticeHash(transformation, frac_coords)
    epcg = EquivalentPointClusterGenerator(
        frac_coords=frac_coords,
        rotations=rotations,
        translations=translations,
        converter=converter,
    )
    return epcg


@pytest.fixture
def epcg_square2d():
    frac_coords = np.array([[0, 0]])
    rotations, translations = square2d_lattice_symmetry()
    transformation = np.eye(2, dtype=int)

    converter = DerivativeMultiLatticeHash(transformation, frac_coords)
    epcg = EquivalentPointClusterGenerator(frac_coords, rotations, translations, converter)
    return epcg


def test_normalize_point_cluster(epcg_fcc):
    pcl = PointCluster([DerivativeSite(1, (0, 0, 0)), DerivativeSite(0, (1, 1, 0))])
    normed = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(1, (-1, -1, 0))])
    assert epcg_fcc.normalize_point_cluster(pcl) == normed


def test_fcc(epcg_fcc):
    distinct_point_clusters = [
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 0, 0))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, -1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 2))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (2, 1, -1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, 1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 2))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (2, 2, -2))]),
    ]
    num_equiv_clusters = [
        6,
        3,
        12,
        6,
        12,
        4,
        24,
        3,
    ]

    for pcl, cnt in zip(distinct_point_clusters, num_equiv_clusters):
        equivs = epcg_fcc.find_equivalent_point_clusters(pcl)
        assert len(equivs) == cnt

    fcc = Structure(Lattice([[0, 1, 1], [1, 0, 1], [1, 1, 0]]), [DummySpecie("X")], [[0, 0, 0]])
    grouped_point_clusters, _ = epcg_fcc.get_all_clusters(fcc, cutoff=4.1, order=2)

    assert len(grouped_point_clusters) == len(num_equiv_clusters)
    for equivs, cnt in zip(grouped_point_clusters, num_equiv_clusters):
        assert len(equivs) == cnt

    # TODO: test order >= 3

    # test for finding distinct point clusters within cutoff
    # upto the 8th NN
    distinct_point_clusters_actual = epcg_fcc.get_distinct_point_clusters(
        fcc, cutoff=4.01, order=2
    )
    grouped_point_clusters_actual = [
        epcg_fcc.find_equivalent_point_clusters(pcl) for pcl in distinct_point_clusters_actual
    ]
    for expect, group in zip(num_equiv_clusters, grouped_point_clusters_actual):
        assert len(group) == expect


def test_bcc(epcg_bcc):
    distinct_point_clusters = [
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 1))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, 2))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 2, 2))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (2, 2, 2))]),
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 2, 2))]),
    ]
    num_equiv_clusters = [
        4,
        3,
        6,
        12,
        4,
        3,
    ]

    for pcl, cnt in zip(distinct_point_clusters, num_equiv_clusters):
        equivs = epcg_bcc.find_equivalent_point_clusters(pcl)
        assert len(equivs) == cnt

    bcc = Structure(Lattice([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]), [DummySpecie("X")], [[0, 0, 0]])
    grouped_point_clusters, _ = epcg_bcc.get_all_clusters(bcc, cutoff=4.1, order=2)

    assert len(grouped_point_clusters) == len(num_equiv_clusters)
    for equivs, cnt in zip(grouped_point_clusters, num_equiv_clusters):
        assert len(equivs) == cnt


@pytest.mark.skip(reason="TODO")
def test_hcp(epcg_hcp):
    pass


def test_perovskite_oxygen_substructure(epcg_perovskite_oxygen_substructure):
    lattice = Lattice(np.eye(3))
    # fmt: off
    oxygen_substructure = Structure(lattice,
                                    ['O', 'O', 'O'],
                                    [
                                        [0.5, 0, 0],      # X=O
                                        [0, 0.5, 0],      # X=O
                                        [0, 0, 0.5],      # X=O
                                    ])
    # fmt: on

    # singlet cluster
    grouped_singlets, _ = epcg_perovskite_oxygen_substructure.get_all_clusters(
        oxygen_substructure, cutoff=1e-8, order=1
    )
    assert len(grouped_singlets) == 1
    assert len(grouped_singlets[0]) == 3  # O1, O2, O3 sites are equivalent

    # pair cluster
    grouped_point_clusters, sizes = epcg_perovskite_oxygen_substructure.get_all_clusters(
        oxygen_substructure, cutoff=2.01, order=2
    )
    sizes_expected = [
        np.sqrt(2) / 2,
        1,
        1,
        np.sqrt(6) / 2,
        np.sqrt(2),
        np.sqrt(2),
        np.sqrt(10) / 2,
        np.sqrt(3),
        np.sqrt(14) / 2,
        2,
        2,
    ]
    assert np.allclose(sizes, sizes_expected)


def test_square2d(epcg_square2d):
    distinct_point_clusters = [
        PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (0, 1))]),
        PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (1, 1))]),
    ]
    num_equiv_clusters = [
        2,
        2,
    ]

    for pcl, cnt in zip(distinct_point_clusters, num_equiv_clusters):
        equivs = epcg_square2d.find_equivalent_point_clusters(pcl)
        assert len(equivs) == cnt
