from math import factorial

from pyzdd import Universe, Combination, enumerate_sets


def test_combination():
    n = 4
    k = 2
    spec = Combination(n, k)
    universe = Universe(n)
    universe.zddSubset(spec)
    universe.zddReduce()
    size = universe.size()

    count = 0
    items_expect = [
        set([1, 2]),
        set([1, 3]),
        set([2, 3]),
        set([1, 4]),
        set([2, 4]),
        set([3, 4]),
    ]
    for items, expect in zip(enumerate_sets(universe), items_expect):
        assert items == expect
        count += 1

    count_expect = factorial(n) // factorial(k) // factorial(n - k)
    assert count == count_expect
