from math import factorial

from pyzdd import Universe, Combination, Choice, enumerate_sets


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


def test_choice():
    n = 4
    k = 2
    v = [0, 2, 3]

    spec = Choice(n, k, v)
    universe = Universe(n)
    universe.zddSubset(spec)
    universe.zddReduce()

    count = 0
    for items in enumerate_sets(universe):
        count += 1
        print(items)

    count_expect = factorial(len(v)) // factorial(k) // factorial(len(v) - k) * (2 ** (n - len(v)))
    assert count == count_expect
