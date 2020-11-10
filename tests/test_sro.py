from pyzdd import (
    Universe,
    Graph,
    add_undirected_edge,
    VertexGraphFrontierManager,
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
    cluster_graph = Graph(num_sites)
    add_undirected_edge(cluster_graph, 0, 1, 2)
    add_undirected_edge(cluster_graph, 1, 2, 2)
    add_undirected_edge(cluster_graph, 2, 3, 2)
    add_undirected_edge(cluster_graph, 3, 0, 2)
    vgfm = VertexGraphFrontierManager(cluster_graph)
    target = 2

    construct_binary_derivative_structures_with_sro(dd, num_sites, num_types, composition_constraints, vgfm, target)
    assert dd.cardinality() == "4"

    actual = set()
    for labeling in enumerate_labelings_with_graph(dd, num_types, cluster_graph):
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
