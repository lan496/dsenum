from pyzdd import Universe, Combination


def test_combination():
    n = 10
    k = 5
    spec = Combination(n, k)
    universe = Universe(n)
    universe.zddSubset(spec)
