import numpy as np
from pymatgen.analysis.structure_analyzer import SpacegroupAnalyzer
from pymatgen.analysis.structure_matcher import StructureMatcher
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import DummySpecie
from pymatgen.io.cif import CifWriter
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor


def get_lattice(kind):
    if kind == 'hcp':
        latt = Lattice(np.array([[1, 0, 0],
                                 [0.5, np.sqrt(3) / 2, 0],
                                 [0, 0, 2 * np.sqrt(6) / 3]]))
        coords = [[0, 0, 0], [1 / 3, 1 / 3, 0.5]]
    else:
        coords = [[0, 0, 0]]

    if kind == 'sc':
        latt = Lattice(np.eye(3))
    elif kind == 'fcc':
        latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
    elif kind == 'bcc':
        latt = Lattice(np.array([[-1, 1, 1], [1, -1, 1], [1, 1, -1]]))
    elif kind == 'hex':
        latt = Lattice.hexagonal(1, 2 * np.sqrt(6) / 3)
    elif kind == 'tet':
        latt = Lattice(np.diag([1, 1, 1.2]))

    struct = Structure(latt, [DummySpecie('X')] * len(coords), coords)
    return struct


def get_fcc_with_vacancy():
    latt = Lattice(np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]]))
    displacement_set = [[0, 0, 0],  # lattice point
                        [0.25, 0.25, 0.25],  # tetrahedral site
                        [0.5, 0.5, 0.5],  # octahedral site
                        [0.75, 0.75, 0.75],  # tetrahedral site
                        ]
    struct = Structure(latt, [DummySpecie('X')] * len(displacement_set), displacement_set)
    return struct


def get_symmetry_operations(structure):
    """
    find symmetry operations for a given structure

    Parameters
    ----------
    structure: pymatgen.core.Structure

    Returns
    -------
    rotations: array, (# of symmetry operations, 3, 3)
    translations: array, (# of symmetry operations, 3)
    """
    sym_dataset = SpacegroupAnalyzer(structure).get_symmetry_dataset()
    rotations = sym_dataset['rotations']
    translations = sym_dataset['translations']
    return rotations, translations


def unique_structures(structures):
    """
    unique list of structure by StructureMatcher

    Parameters
    ----------
    structures: list of pymatgen.core.Structure objects

    Returns
    -------
    uniqued: list of pymatgen.core.Structure objects
    """
    uniqued = []
    stm = StructureMatcher()

    for struct in structures:
        rotations, translations = get_symmetry_operations(struct)

        is_unique = True
        for r, t in zip(rotations, translations):
            new_frac_coords = np.dot(r, struct.frac_coords.T) + t[:, np.newaxis]
            struct_tmp = Structure(struct.lattice, struct.species, new_frac_coords.T)
            if any([stm.fit(struct_tmp, unique_struct) for unique_struct in uniqued]):
                is_unique = False
                break

        if is_unique:
            uniqued.append(struct)

    return uniqued


# https://repl.it/@smichr/msp
def msp(items):
    '''
    Yield the permutations of `items` where items is either a list
    of integers representing the actual items or a list of hashable items.
    The output are the unique permutations of the items given as a list
    of integers 0, ..., n-1 that represent the n unique elements in
    `items`.
     Examples
    ========
     >>> for i in msp('xoxox'):
    ...   print(i)
     [1, 1, 1, 0, 0]
    [0, 1, 1, 1, 0]
    [1, 0, 1, 1, 0]
    [1, 1, 0, 1, 0]
    [0, 1, 1, 0, 1]
    [1, 0, 1, 0, 1]
    [0, 1, 0, 1, 1]
    [0, 0, 1, 1, 1]
    [1, 0, 0, 1, 1]
    [1, 1, 0, 0, 1]
     Reference: "An O(1) Time Algorithm for Generating Multiset Permutations", Tadao Takaoka
    https://pdfs.semanticscholar.org/83b2/6f222e8648a7a0599309a40af21837a0264b.pdf
    '''
    def visit(head):
        (rv, j) = ([], head)
        for i in range(N):
            (dat, j) = E[j]
            rv.append(dat)
        return rv
    u = list(set(items))
    E = list(reversed(sorted([u.index(i) for i in items])))
    N = len(E)
    # put E into linked-list format
    (val, nxt) = (0, 1)
    for i in range(N):
        E[i] = [E[i], i + 1]
    E[-1][nxt] = None
    head = 0
    afteri = N - 1
    i = afteri - 1
    yield visit(head)
    while E[afteri][nxt] is not None or E[afteri][val] < E[head][val]:
        j = E[afteri][nxt]  # added to algorithm for clarity
        if j is not None and E[i][val] >= E[j][val]:
            beforek = afteri
        else:
            beforek = i
        k = E[beforek][nxt]
        E[beforek][nxt] = E[k][nxt]
        E[k][nxt] = head
        if E[k][val] < E[head][val]:
            i = k
        afteri = E[i][nxt]
        head = k
        yield visit(head)


def write_cif(filename, struct, refine_cell=False, resize_volume=False):
    struct = refine_and_resize_structure(struct, refine_cell, resize_volume)
    if not struct.is_valid(1e-4):
        return
    cw = CifWriter(struct)
    cw.write_file(filename)


def refine_and_resize_structure(struct, refine_cell=True, resize_volume=True):
    if resize_volume:
        dls = DLSVolumePredictor()
        struct = dls.get_predicted_structure(struct)
        struct.apply_strain(0.5)

    if refine_cell:
        sga = SpacegroupAnalyzer(struct, symprec=1e-6, angle_tolerance=1e-2)
        struct = sga.get_primitive_standard_structure()

    return struct
