import numpy as np
from pymatgen.core import Structure, Composition
from pymatgen.analysis.structure_matcher import StructureMatcher

from dsenum.site import DerivativeSite
from dsenum.enumerate import StructureEnumerator
from dsenum.utils import get_lattice
from dsenum.cluster.point_cluster import PointCluster
from dsenum.cluster.sro import SROStructureEnumerator


def test_sro_L1_1():
    fcc = get_lattice("fcc")
    composition_ratios = [[1, 1]]
    num_types = 2

    index = 2
    transformation = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [1, 1, 2],
        ]
    )

    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
        # 2nd NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, -1))]),
    ]
    sse = SROStructureEnumerator(fcc, index, num_types, composition_ratios, clusters)
    labelings = sse.generate_with_hnf(transformation)

    # only 1st NN
    clusters = [
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
    ]

    sse = SROStructureEnumerator(fcc, index, num_types, composition_ratios, clusters)
    labelings = sse.generate_with_hnf(transformation)
    assert labelings == [[0, 1]]


def get_fcc_structure(a: float) -> Structure:
    matrix = (
        a
        / 2
        * np.array(
            [
                [0, 1, 1],
                [1, 0, 1],
                [1, 1, 0],
            ]
        )
    )
    frac_coords = np.zeros((1, 3))
    fcc = Structure(matrix, ["Cu"], frac_coords)
    return fcc


def count_bonds(structure: Structure, specie0, specie1, pair_length: float):
    """
    explicitly count bonds between `specie0` and `specie1` in `structure`
    """
    disp = 1e-2  # some small posivive value for finding neighbors

    num_bonds = 0
    for site in structure.sites:
        if str(site.specie) != specie0:
            continue
        neighbors = structure.get_neighbors(site, r=pair_length + disp)
        for nsite in neighbors:
            if str(nsite.specie) != specie1:
                continue
            if np.isclose(nsite.nn_distance, pair_length):
                num_bonds += 1
    return num_bonds / 2


def get_fcc_sqs4_zunger(a: float) -> Structure:
    matrix = a * np.array(
        [
            [-0.5, 0.5, 0],
            [0, 0, 1],
            [1, 1, 0],
        ]
    )
    species = ["Cu"] * 2 + ["Pt"] * 2
    coords = a * np.array(
        [
            [-0.125, -0.125, 0],
            [-0.125, 0.375, 0.5],
            [0.375, 0.375, 0],
            [-0.625, -0.125, 0.5],
        ]
    )
    sqs4 = Structure(matrix, species, coords, coords_are_cartesian=True)
    return sqs4


def get_fcc_sqs8_zunger(a: float) -> Structure:
    matrix = a * np.array(
        [
            [1, 0.5, -0.5],
            [0.5, -0.5, 0],
            [1, 1, 2],
        ]
    )
    species = ["Cu"] * 4 + ["Pt"] * 4
    coords = a * np.array(
        [
            [0, 0, 0],
            [0.5, 0.5, 0],
            [0.5, 0, 1.5],
            [0, 0, 2],
            [0.5, 0, 0.5],
            [0, 0, 1],
            [0.5, 0.5, 1],
            [0.5, 0.5, 2],
        ]
    )
    sqs8 = Structure(matrix, species, coords, coords_are_cartesian=True)
    return sqs8


def test_fcc_sqs4():
    a = 3.0
    fcc = get_fcc_structure(a)

    index = 4
    num_types = 2
    mapping_color_species = ["Cu", "Pt"]
    concentration = [0.5, 0.5]
    target_composition = Composition.from_dict(
        {
            mapping_color_species[0]: index * concentration[0],
            mapping_color_species[1]: index * concentration[1],
        }
    )

    se = StructureEnumerator(
        fcc,
        index,
        num_types,
        mapping_color_species=mapping_color_species,
        color_exchange=False,
        remove_superperiodic=True,
        remove_incomplete=True,
    )
    dstructs = se.generate()

    # 1st NN
    length_1st = a / np.sqrt(2)
    target_num_bonds_1st = 6 * index * concentration[1] * concentration[1]

    # enumerate SQS by brute force
    list_sqs = []
    for struct in dstructs:
        if not struct.composition.almost_equals(target_composition):
            continue
        num_bonds_1st = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_1st
        )
        if np.isclose(num_bonds_1st, target_num_bonds_1st):
            list_sqs.append(struct)

    stm = StructureMatcher()

    # check SQS-4 in the original article is contained
    sqs4 = get_fcc_sqs4_zunger(a)
    assert any([stm.fit(struct, sqs4) for struct in list_sqs])

    # enumerate SQS by DD
    composition_ratios = [[1, 1]]
    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
    ]
    sse = SROStructureEnumerator(
        fcc,
        index,
        num_types,
        composition_ratios,
        clusters,
        mapping_color_species=mapping_color_species,
        remove_superperiodic=True,
    )
    list_sqs_dd = sse.generate_structures()

    # check DD and brute force return the same results
    assert len(list_sqs_dd) == len(list_sqs)
    for s1 in list_sqs_dd:
        assert any([stm.fit(s1, s2) for s2 in list_sqs])


def test_fcc_sqs8():
    a = 3.0
    fcc = get_fcc_structure(a)

    index = 8
    num_types = 2
    mapping_color_species = ["Cu", "Pt"]
    concentration = [0.5, 0.5]
    target_composition = Composition.from_dict(
        {
            mapping_color_species[0]: index * concentration[0],
            mapping_color_species[1]: index * concentration[1],
        }
    )

    se = StructureEnumerator(
        fcc,
        index,
        num_types,
        mapping_color_species=mapping_color_species,
        color_exchange=False,
        remove_superperiodic=True,
        remove_incomplete=True,
    )
    dstructs = se.generate()

    # 1st NN
    length_1st = a / np.sqrt(2)
    target_num_bonds_1st = 6 * index * concentration[1] * concentration[1]
    # 2nd NN
    length_2nd = a
    target_num_bonds_2nd = 3 * index * concentration[1] * concentration[1]

    # enumerate SQS by brute force
    list_sqs = []
    for struct in dstructs:
        if not struct.composition.almost_equals(target_composition):
            continue
        num_bonds_1st = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_1st
        )
        num_bonds_2nd = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_2nd
        )
        if np.isclose(num_bonds_1st, target_num_bonds_1st) and np.isclose(
            num_bonds_2nd, target_num_bonds_2nd
        ):
            list_sqs.append(struct)

    stm = StructureMatcher()

    # check SQS-8 in the original article is contained
    sqs8 = get_fcc_sqs8_zunger(a)
    assert any([stm.fit(struct, sqs8) for struct in list_sqs])

    # enumerate SQS by DD
    composition_ratios = [[1, 1]]
    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
        # 2nd NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, -1))]),
    ]
    sse = SROStructureEnumerator(
        fcc,
        index,
        num_types,
        composition_ratios,
        clusters,
        mapping_color_species=mapping_color_species,
        remove_superperiodic=True,
    )
    list_sqs_dd = sse.generate_structures()

    # check DD and brute force return the same results
    assert len(list_sqs_dd) == len(list_sqs)
    for s1 in list_sqs_dd:
        assert any([stm.fit(s1, s2) for s2 in list_sqs])


def get_bcc_structure(a: float) -> Structure:
    matrix = (
        a
        / 2
        * np.array(
            [
                [-1, 1, 1],
                [1, -1, 1],
                [1, 1, -1],
            ]
        )
    )
    frac_coords = np.zeros((1, 3))
    bcc = Structure(matrix, ["Mo"], frac_coords)
    return bcc


def get_bcc_AB_sqs8(a: float) -> Structure:
    matrix = a * np.array(
        [
            [0.5, 0.5, -1.5],
            [1.5, 0.5, -0.5],
            [0.0, -2.0, 0.0],
        ]
    )
    species = ["Mo"] * 4 + ["Nb"] * 4
    coords = a * np.array(
        [
            [2.0, 0.0, -2.0],
            [0.5, -1.5, -0.5],
            [1.0, -1.0, -1.0],
            [1.5, -0.5, -1.5],
            [2.0, -1.0, -2.0],
            [0.5, -0.5, -0.5],
            [1.0, 0.0, -1.0],
            [1.5, 0.5, -1.5],
        ]
    )
    sqs8 = Structure(matrix, species, coords, coords_are_cartesian=True)
    return sqs8


def get_bcc_AB3_sqs8(a: float) -> Structure:
    matrix = a * np.array(
        [
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, -1.0],
            [0.0, -2.0, -2.0],
        ]
    )
    species = ["Mo"] * 2 + ["Nb"] * 6
    coords = a * np.array(
        [
            [-0.5, -0.5, -1.5],
            [-1.0, -1.0, -2.0],
            [-0.5, -1.5, -2.5],
            [-0.5, 0.5, -1.5],
            [-1.0, -1.0, -3.0],
            [-1.0, 0.0, -1.0],
            [-0.5, -0.5, -2.5],
            [-1.0, 0.0, -2.0],
        ]
    )
    sqs8 = Structure(matrix, species, coords, coords_are_cartesian=True)
    return sqs8


def test_bcc_AB_sqs8():
    a = 3.0
    bcc = get_bcc_structure(a)

    index = 8
    num_types = 2
    mapping_color_species = ["Mo", "Nb"]
    concentration = [0.5, 0.5]
    target_composition = Composition.from_dict(
        {
            mapping_color_species[0]: index * concentration[0],
            mapping_color_species[1]: index * concentration[1],
        }
    )

    se = StructureEnumerator(
        bcc,
        index,
        num_types,
        mapping_color_species=mapping_color_species,
        color_exchange=False,
        remove_superperiodic=True,
        remove_incomplete=True,
    )
    dstructs = se.generate()

    # 1st NN
    length_1st = a / 2 * np.sqrt(3)
    target_num_bonds_1st = 4 * index * concentration[1] * concentration[1]
    # 2nd NN
    length_2nd = a
    target_num_bonds_2nd = 3 * index * concentration[1] * concentration[1]

    # enumerate SQS by brute force
    list_sqs = []
    for struct in dstructs:
        if not struct.composition.almost_equals(target_composition):
            continue
        num_bonds_1st = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_1st
        )
        num_bonds_2nd = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_2nd
        )
        if np.isclose(num_bonds_1st, target_num_bonds_1st) and np.isclose(
            num_bonds_2nd, target_num_bonds_2nd
        ):
            list_sqs.append(struct)

    stm = StructureMatcher()

    # check SQS-8 in the original article is contained
    sqs8 = get_bcc_AB_sqs8(a)
    assert any([stm.fit(struct, sqs8) for struct in list_sqs])

    # enumerate SQS by DD
    composition_ratios = [[1, 1]]
    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
        # 2nd NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 1))]),
    ]
    sse = SROStructureEnumerator(
        bcc,
        index,
        num_types,
        composition_ratios,
        clusters,
        mapping_color_species=mapping_color_species,
        remove_superperiodic=True,
    )
    list_sqs_dd = sse.generate_structures()

    # check DD and brute force return the same results
    assert len(list_sqs_dd) == len(list_sqs)
    for s1 in list_sqs_dd:
        assert any([stm.fit(s1, s2) for s2 in list_sqs])


def test_bcc_AB3_sqs8():
    a = 3.0
    bcc = get_bcc_structure(a)

    index = 8
    num_types = 2
    mapping_color_species = ["Mo", "Nb"]
    concentration = [0.25, 0.75]
    target_composition = Composition.from_dict(
        {
            mapping_color_species[0]: index * concentration[0],
            mapping_color_species[1]: index * concentration[1],
        }
    )

    se = StructureEnumerator(
        bcc,
        index,
        num_types,
        mapping_color_species=mapping_color_species,
        color_exchange=False,
        remove_superperiodic=True,
        remove_incomplete=True,
    )
    dstructs = se.generate()

    # 1st NN
    length_1st = a / 2 * np.sqrt(3)
    target_num_bonds_1st = 4 * index * concentration[1] * concentration[1]

    # enumerate SQS by brute force
    list_sqs = []
    for struct in dstructs:
        if not struct.composition.almost_equals(target_composition):
            continue
        num_bonds_1st = count_bonds(
            struct, mapping_color_species[1], mapping_color_species[1], length_1st
        )
        if np.isclose(num_bonds_1st, target_num_bonds_1st):
            list_sqs.append(struct)

    stm = StructureMatcher()

    # check SQS-8 in the original article is contained
    sqs8 = get_bcc_AB3_sqs8(a)
    assert any([stm.fit(struct, sqs8) for struct in list_sqs])

    # enumerate SQS by DD
    composition_ratios = [[1, 3]]
    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
    ]
    sse = SROStructureEnumerator(
        bcc,
        index,
        num_types,
        composition_ratios,
        clusters,
        mapping_color_species=mapping_color_species,
        remove_superperiodic=True,
    )
    list_sqs_dd = sse.generate_structures()

    # check DD and brute force return the same results
    assert len(list_sqs_dd) == len(list_sqs)
    for s1 in list_sqs_dd:
        assert any([stm.fit(s1, s2) for s2 in list_sqs])
