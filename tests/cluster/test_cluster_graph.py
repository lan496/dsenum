import networkx as nx
import numpy as np

from pymatgen.core.periodic_table import Specie

from dsenum.site import DerivativeSite
from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.cluster.cluster_graph import BinaryPairClusterGraph
from dsenum.cluster.point_cluster import EquivalentPointClusterGenerator, PointCluster
from dsenum.utils import square2d_lattice_symmetry, get_lattice


def test_pair_cluster_graph_square2d():
    frac_coords = np.array([[0, 0]])
    rotations, translations = square2d_lattice_symmetry()
    # 2x2 supercell
    transformation = np.diag([2, 2])

    converter = DerivativeMultiLatticeHash(transformation, frac_coords)
    epcg = EquivalentPointClusterGenerator(frac_coords, rotations, translations, converter)

    # prepare test cases
    graph1NN = nx.Graph()
    graph1NN.add_nodes_from([0, 1, 2, 3])
    graph1NN.add_weighted_edges_from(
        [
            (0, 1, 2),
            (0, 2, 2),
            (1, 3, 2),
            (2, 3, 2),
        ]  # (src, dst, weight)
    )

    graph2NN = nx.Graph()
    graph2NN.add_nodes_from([0, 1, 2, 3])
    graph2NN.add_weighted_edges_from(
        [
            (0, 3, 4),
            (1, 2, 4),
        ]
    )

    testcases = {
        "1st-A3B": {
            "cluster": PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (0, 1))]),
            "graph": graph1NN,
            "composition_ratio": [[3, 1]],
            "loop_offset": 0.0,
            "subtest": [
                {
                    "labeling": [0, 0, 0, 1],
                    "corr": 0.0,
                    "sro": -1.0 / 3.0,
                }
            ],
        },
        "1st-AB": {
            "cluster": PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (0, 1))]),
            "graph": graph1NN,
            "composition_ratio": [[1, 1]],
            "loop_offset": 0.0,
            "subtest": [
                {
                    "labeling": [0, 0, 1, 1],
                    "corr": 0.25,
                    "sro": 0.0,
                },
                {
                    "labeling": [0, 1, 1, 0],
                    "corr": 0.0,
                    "sro": -1.0,
                },
            ],
        },
        "2nd-AB": {
            "cluster": PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (1, 1))]),
            "graph": graph2NN,
            "composition_ratio": [[1, 1]],
            "loop_offset": 0.0,
            "subtest": [
                {
                    "labeling": [0, 0, 1, 1],
                    "corr": 0.0,
                    "sro": -1.0,
                },
                {
                    "labeling": [0, 1, 1, 0],
                    "corr": 0.5,
                    "sro": 1.0,
                },
            ],
        },
    }

    for name, data in testcases.items():
        cluster = data["cluster"]
        composition_ratio = data["composition_ratio"]
        grouped_cluster = epcg.find_equivalent_point_clusters(cluster)
        bpcg = BinaryPairClusterGraph(converter, grouped_cluster, composition_ratio)

        graph_actual = bpcg.graph
        graph_expect = data["graph"]
        assert set(graph_actual.nodes) == set(graph_expect.nodes)
        assert set(graph_actual.edges(data="weight")) == set(graph_expect.edges(data="weight"))

        offset_actual = bpcg.loop_offset
        offset_expect = data["loop_offset"]
        assert np.isclose(offset_actual, offset_expect)

        # test correlations
        for subdata in data["subtest"]:
            labeling = subdata["labeling"]
            corr_expect = subdata["corr"]
            sro_expect = subdata["sro"]
            corr_actual = bpcg.calc_correlation(labeling)
            assert np.isclose(corr_actual, corr_expect)

            sro_actual = bpcg.calc_short_range_order(labeling)
            assert np.isclose(sro_actual, sro_expect)


def test_loop_offset():
    frac_coords = np.array([[0, 0]])
    rotations, translations = square2d_lattice_symmetry()
    # 1x2 supercell
    transformation = np.diag([1, 2])

    converter = DerivativeMultiLatticeHash(transformation, frac_coords)
    epcg = EquivalentPointClusterGenerator(frac_coords, rotations, translations, converter)

    cluster = PointCluster([DerivativeSite(0, (0, 0)), DerivativeSite(0, (0, 1))])
    grouped_cluster = epcg.find_equivalent_point_clusters(cluster)
    composition_ratio = [[1, 1]]
    bpcg = BinaryPairClusterGraph(converter, grouped_cluster, composition_ratio)

    # test loop offset
    loop_offset_expect = 1.0
    assert np.isclose(bpcg.loop_offset, loop_offset_expect)

    # test graph
    graph1NN = nx.Graph()
    graph1NN.add_nodes_from([0, 1])
    graph1NN.add_weighted_edges_from(
        [
            (0, 1, 2),
        ]  # (src, dst, weight)
    )
    graph_actual = bpcg.graph
    assert set(graph_actual.nodes) == set(graph1NN.nodes)
    assert set(graph_actual.edges(data="weight")) == set(graph1NN.edges(data="weight"))


def test_fcc_cluster_graph():
    aristo = get_lattice("fcc")
    composition_ratio = [
        [1, 1],
    ]

    # L1_0 structure
    transformation = np.array(
        [
            [1, 0, 0],
            [1, 2, 0],
            [0, 0, 1],
        ]
    )
    epcg = EquivalentPointClusterGenerator.from_structure(aristo, transformation)
    converter = epcg.get_converter()
    cluster = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))])
    grouped_clusters = epcg.find_equivalent_point_clusters(cluster)
    bpcg = BinaryPairClusterGraph(converter, grouped_clusters, composition_ratio)
    L1_0 = [0, 1]  # AuCu structure
    assert np.isclose(bpcg.loop_offset, 2.0)
    assert np.isclose(bpcg.get_sqs_target_value(), 1.0)
    assert np.isclose(bpcg.calc_correlation(L1_0), 2.0 / 12)
    assert np.isclose(bpcg.calc_short_range_order(L1_0), -1.0 / 3.0)

    # L1_1 structure
    transformation = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 2],
        ]
    )
    epcg = EquivalentPointClusterGenerator.from_structure(aristo, transformation)
    converter = epcg.get_converter()
    cluster = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))])
    grouped_clusters = epcg.find_equivalent_point_clusters(cluster)
    bpcg = BinaryPairClusterGraph(converter, grouped_clusters, composition_ratio)
    L1_1 = [0, 1]  # AuCu structure
    assert np.isclose(bpcg.loop_offset, 3.0)
    assert np.isclose(bpcg.get_sqs_target_value(), 0.0)
    assert np.isclose(bpcg.calc_correlation(L1_1), 3.0 / 12)
    assert np.isclose(bpcg.calc_short_range_order(L1_1), 0.0)
