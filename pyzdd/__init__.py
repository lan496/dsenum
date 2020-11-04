from typing import Generator, Set, List

from _pyzdd import (
    Universe,
    Combination,
    Choice,
    Permutation,
    generate_permutation_group,
    variable_choice,
    construct_derivative_structures,
    convert_to_labeling,
)
from ._version import get_versions  # type: ignore


__version__ = get_versions()["version"]
del get_versions


def enumerate_sets(universe: Universe, n: int) -> Generator[Set[int], None, None]:
    """
    yield combinations

    Parameters
    ----------
    universe
    n: the number of variables

    Returns
    -------
    items: set of selected levels(1-indexed)
    """
    if not isinstance(universe, Universe):
        TypeError(f"Given type is not Universe but {type(universe)}")

    itr = universe.begin()
    end = universe.end()
    while itr != end:
        choice = variable_choice(itr, n)
        yield choice
        itr.next()


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
