from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry
)
from labeling import Labeling
from derivative_structure import DerivativeStructure


def enumerate_derivative_structures(structure, index, num_type):
    A = structure.lattice.matrix.T
    list_HNF = generate_all_superlattices(index)
    sym_dataset = SpacegroupAnalyzer(structure).get_symmetry_dataset()
    rotations = sym_dataset['rotations']
    list_reduced_HNF = \
        reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)

    list_ds = []

    for hnf in list_reduced_HNF:
        labeling = Labeling(hnf, num_type, rotations)
        lbls_tmp = labeling.get_inequivalent_labelings()
        list_ds.extend([DerivativeStructure(hnf, num_type, A, lbl)
                        for lbl in lbls_tmp])

    return list_ds


if __name__ == '__main__':
    from derivative_structure import get_lattice
    structure = get_lattice('fcc')
    index = 12
    num_type = 2

    list_ds = enumerate_derivative_structures(structure, index, num_type)
    print(3875, len(list_ds))
