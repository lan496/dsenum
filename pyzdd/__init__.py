from typing import Generator, Set

from _pyzdd import (
    Universe,
    Combination,
    Choice,
    variable_choice,
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
