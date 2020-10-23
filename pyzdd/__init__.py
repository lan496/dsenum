from _pyzdd import Universe, Combination
from ._version import get_versions  # type: ignore


__version__ = get_versions()["version"]
del get_versions


def enumerate_sets(universe: Universe):
    if not isinstance(universe, Universe):
        TypeError(f"Given type is not Universe but {type(universe)}")

    itr = universe.begin()
    end = universe.end()
    while itr != end:
        items = itr.itemset()
        itr.next()
        yield items
