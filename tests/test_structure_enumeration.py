from pyzdd import (
    Universe,
    Permutation,
    generate_permutation_group,
    construct_derivative_structures,
    enumerate_labelings,
)

def test_structure_enumeration():
    num_sites = 4
    num_types = 3

    c4 = Permutation([1, 2, 3, 0])
    m = Permutation([3, 2, 1, 0])
    automorphism = generate_permutation_group([c4, m])
    translations = generate_permutation_group([c4])

    dd = Universe()
    construct_derivative_structures(dd, num_sites, num_types, automorphism)
    assert dd.cardinality() == "21"

    actual = set()
    for labeling in enumerate_labelings(dd, num_sites, num_types):
        actual.add(tuple(labeling))

    list_expect = [
        [0, 0, 0, 0],
        [1, 1, 1, 1],
        [2, 2, 2, 2],
        #
        [1, 1, 1, 0],
        [1, 1, 0, 0],
        [1, 0, 1, 0],
        [1, 0 ,0, 0],
        #
        [2, 2, 2, 0],
        [2, 2, 0, 0],
        [2, 0, 2, 0],
        [2, 0, 0, 0],
        #
        [2, 2, 2, 1],
        [2, 2, 1, 1],
        [2, 1, 2, 1],
        [2, 1, 1, 1],
        #
        [2, 2, 1, 0],
        [2, 1, 2, 0],
        [2, 1, 1, 0],
        [2, 1, 0, 1],
        [2, 0, 1, 0],
        [2, 1, 0, 0],
    ]
    expect = set()
    for labeling in list_expect:
        expect.add(tuple(labeling))

    assert actual == expect

    # remove superperiodic structures
    dd = Universe()
    construct_derivative_structures(
        dd,
        num_sites,
        num_types,
        automorphism,
        translations,
        remove_superperiodic=True
    )
    assert dd.cardinality() == "15"

    actual = set()
    for labeling in enumerate_labelings(dd, num_sites, num_types):
        actual.add(tuple(labeling))

    list_expect = [
        [1, 1, 1, 0],
        [1, 1, 0, 0],
        [1, 0 ,0, 0],
        #
        [2, 2, 2, 0],
        [2, 2, 0, 0],
        [2, 0, 0, 0],
        #
        [2, 2, 2, 1],
        [2, 2, 1, 1],
        [2, 1, 1, 1],
        #
        [2, 2, 1, 0],
        [2, 1, 2, 0],
        [2, 1, 1, 0],
        [2, 1, 0, 1],
        [2, 0, 1, 0],
        [2, 1, 0, 0],
    ]
    expect = set()
    for labeling in list_expect:
        expect.add(tuple(labeling))

    assert actual == expect
