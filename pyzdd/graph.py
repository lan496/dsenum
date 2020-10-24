from typing import List, Dict, Union, Any

import networkx as nx

from _pyzdd import Edge, Graph, GraphAuxiliary


def convert_to_raw_graph(graph: nx.Graph) -> Union[List[List[Edge]], Dict[Any, int]]:
    """
    convert nx.Graph to List[List[Edge]]

    Parameters
    ----------
    graph: nx.Graph
    mapping: map original nodes to 0-indexed relabeled ones.
    """
    mapping = {}
    for i, u in enumerate(sorted(graph.nodes)):
        mapping[u] = i
    relabeled = nx.relabel_nodes(graph, mapping)

    graph_pb = []
    for src in relabeled.nodes:
        edges_pb = [Edge(src, dst, relabeled[src][dst].get("weight", 1)) for dst in relabeled.neighbors(src)]
        graph_pb.append(edges_pb)

    return graph_pb, mapping
