from typing import Generator, List

from pyzdd import Universe
from pyzdd.graph import Graph, VertexGraphFrontierManager

from _pyzdd import (
    construct_derivative_structures,
    construct_binary_derivative_structures_with_sro,
    convert_to_labeling,
    convert_to_labeling_with_graph,
)


def enumerate_labelings(dd: Universe, num_sites, num_types) -> Generator[List[int], None, None]:
    """
    yield labeling of derivative structure

    Parameters
    ----------
    dd: DD for enumerating derivative structures
    num_sites: the number of sites in supercell
    num_types: the number of kinds of species

    Returns
    -------
    labeling
    """
    if not isinstance(dd, Universe):
        TypeError(f"Given type is not Universe but {type(dd)}")

    itr = dd.begin()
    end = dd.end()
    while itr != end:
        labeling = convert_to_labeling(itr, num_sites, num_types)
        yield labeling
        itr.next()


def enumerate_labelings_with_graph(dd: Universe, num_types: int, cluster_graph: Graph) -> Generator[List[int], None, None]:
    """
    yield labeling of derivative structure with a cluster graph

    Parameters
    ----------
    dd: DD for enumerating derivative structures
    num_types: the number of kinds of species
    cluster_graph

    Returns
    -------
    labeling
    """
    if not isinstance(dd, Universe):
        TypeError(f"Given type is not Universe but {type(dd)}")

    itr = dd.begin()
    end = dd.end()
    vgfm = VertexGraphFrontierManager(cluster_graph)
    while itr != end:
        labeling = convert_to_labeling_with_graph(itr, vgfm, num_types)
        yield labeling
        itr.next()
