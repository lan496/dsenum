from typing import List, Set

import networkx as nx
from networkx.generators.lattice import grid_graph
from networkx import path_graph, complete_graph

from pyzdd.graph import GraphAuxiliary, convert_to_raw_graph


def test_on_grid_graph():
    """
       e1    e5
    v0----v1----v2
    |     |     |
    |e0   |e4   |e9
    |     |     |
    v3----v4----v5
    |     |     |
    |e2   |e7   |e11
    |     |     |
    v6----v7----v8
       e6   e10
    """
    N = 3
    grid_2d = grid_graph(dim=[N, N])
    graph, mapping = convert_to_raw_graph(grid_2d)
    graphaux = GraphAuxiliary(graph)

    g = nx.relabel_nodes(grid_2d, mapping)

    frontiers = [
        [],
        [0, 3],
        [1, 3],
        [1, 3, 6],
        [1, 4, 6],
        [1, 4, 6],
        [2, 4, 6],
        [2, 4, 7],
        [2, 4, 7],
        [2, 5, 7],
        [5, 7],
        [5, 8],
    ]
    introduced = [
        [0, 3],
        [1],
        [6],
        [4],
        [],
        [2],
        [7],
        [],
        [5],
        [],
        [8],
        [],
    ]
    forgotten = [
        [],
        [0],
        [],
        [3],
        [],
        [1],
        [6],
        [],
        [4],
        [2],
        [7],
        [5, 8],
    ]
    frontier_size = 3
    vertex_to_position = [0, 0, 0, 1, 1, 1, 2, 2, 0]

    for eid in range(g.size()):
        assert graphaux.frontier(eid) == frontiers[eid]
        assert graphaux.introduced(eid) == introduced[eid]
        assert graphaux.forgotten(eid) == forgotten[eid]
    for u in range(g.order()):
        assert graphaux.map_vertex_to_position(u) == vertex_to_position[u]

    assert graphaux.max_frontier_size == frontier_size
