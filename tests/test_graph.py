from typing import List, Set

import networkx as nx
from networkx.generators.lattice import grid_graph
from networkx import path_graph, complete_graph

from pyzdd.graph import GraphAuxiliary, convert_to_raw_graph


def test_on_grid_graph():
    N = 3
    grid_2d = grid_graph(dim=[N, N])
    graph, mapping = convert_to_raw_graph(grid_2d)
    graphaux = GraphAuxiliary(graph)

    bags_expect = [
        set([0]),
        set([0, 3]),
        set([1, 3]),
        set([1, 3, 6]),
        set([1, 4, 6]),
        set([1, 4, 6]),
        set([2, 4, 6]),
        set([2, 4, 7]),
        set([2, 4, 7]),
        set([2, 5, 7]),
        set([5, 7]),
        set([5, 8]),
        set([]),
    ]

    for bag, expect in zip(graphaux.bags, bags_expect):
        assert bag == expect
    g = nx.relabel_nodes(grid_2d, mapping)

    width_expect = 3
    assert graphaux.width == width_expect
