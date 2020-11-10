import networkx as nx

from pyzdd import Universe
from pyzdd.graph import (
    Graph,
    VertexGraphFrontierManager,
    convert_to_raw_graph
)
from pyzdd.structure import (
    construct_binary_derivative_structures_with_sro,
    enumerate_labelings_with_graph,
)


def test_sro():
    num_sites = 4
    num_types = 2

    dd = Universe()
    composition_constraints = [
        ([0, 1, 2, 3], 2),
    ]

    cluster_graph = nx.Graph()
    cluster_graph.add_nodes_from([0, 1, 2, 3])
    cluster_graph.add_weighted_edges_from([
        (0, 1, 2),
        (1, 2, 2),
        (2, 3, 2),
        (3, 0, 2),
    ])
    raw_graph, _ = convert_to_raw_graph(cluster_graph)
    vgfm = VertexGraphFrontierManager(raw_graph)
    target = 2

    construct_binary_derivative_structures_with_sro(dd, num_sites, num_types, composition_constraints, vgfm, target)
    assert dd.cardinality() == "4"

    actual = set()
    for labeling in enumerate_labelings_with_graph(dd, num_types, raw_graph):
        actual.add(tuple(labeling))

    list_expect = [
        [0, 0, 1, 1],
        [0, 1, 1, 0],
        [1, 0, 0, 1],
        [1, 1, 0, 0],
    ]
    expect = set()
    for labeling in list_expect:
        expect.add(tuple(labeling))

    assert actual == expect
