from typing import List, Set

import networkx as nx
from networkx.generators.lattice import grid_graph
from networkx import path_graph, complete_graph

from pyzdd.graph import (
    GraphAuxiliary,
    VertexGraphFrontierManager,
    convert_to_raw_graph,
)


def test_on_grid_graph():
    """
       e1    e5
    v0----v1----v2
    |     |     |
    |e0   |e4   |e9
    | e3  |  e8 |
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
        [0, 3],
        [0, 1, 3],
        [1, 3, 6],
        [1, 3, 4, 6],
        [1, 4, 6],
        [1, 2, 4, 6],
        [2, 4, 6, 7],
        [2, 4, 7],
        [2, 4, 5, 7],
        [2, 5, 7],
        [5, 7, 8],
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
    frontier_size = 4
    vertex_to_position = [0, 2, 1, 1, 3, 0, 0, 2, 1]

    for eid in range(g.size()):
        assert graphaux.frontier(eid) == frontiers[eid]
        assert graphaux.introduced(eid) == introduced[eid]
        assert graphaux.forgotten(eid) == forgotten[eid]
    for u in range(g.order()):
        assert graphaux.map_vertex_to_position(u) == vertex_to_position[u]

    assert graphaux.max_frontier_size == frontier_size


def test_vertex_order():
    """
       e0
    v0----v1
    |     |
    |e3   |e1
    |     |
    v2----v3
       e2
    """
    graph = nx.Graph()
    graph.add_nodes_from([0, 1, 2, 3])
    graph.add_edges_from([
        (0, 1),
        (0, 2),
        (2, 3),
        (1, 3),
    ])
    raw_graph, _ = convert_to_raw_graph(graph)
    vertex_order = [0, 1, 3, 2]
    vgfm = VertexGraphFrontierManager(raw_graph, vertex_order)
    assert vgfm.get_vertex_order() == vertex_order
