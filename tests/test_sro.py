import networkx as nx

from pyzdd import Universe, Permutation, generate_permutation_group
from pyzdd.graph import (
    Graph,
    convert_to_raw_graph,
    get_vertex_order_by_bfs,
)
from pyzdd.structure import (
    construct_derivative_structures_with_sro,
    enumerate_binary_labelings_with_graph,
)


def test_sro():
    num_sites = 4
    num_types = 2

    c4 = Permutation([1, 2, 3, 0])
    m = Permutation([3, 2, 1, 0])
    automorphism = generate_permutation_group([c4, m])
    translations = generate_permutation_group([c4])

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
    vertex_order = get_vertex_order_by_bfs(raw_graph)
    graphs = [
        raw_graph,
    ]
    targets = [
        2,
    ]

    construct_derivative_structures_with_sro(
        dd,
        num_sites,
        num_types,
        vertex_order,
        automorphism,
        translations,
        composition_constraints,
        graphs,
        targets,
        remove_superperiodic=True,
    )
    assert dd.cardinality() == "1"

    actual = set()
    for labeling in enumerate_binary_labelings_with_graph(dd, num_sites, vertex_order):
        actual.add(tuple(labeling))

    list_expect = [
        [0, 0, 1, 1],
    ]
    expect = set()
    for labeling in list_expect:
        expect.add(tuple(labeling))

    assert actual == expect
