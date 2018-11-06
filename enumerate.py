from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry
)
from labeling import Labeling
from derivative_structure import DerivativeStructure
from utils import get_symmetry_operations


def enumerate_derivative_structures(structure, index, num_type):
    displacement_set = structure.frac_coords
    num_site_parent = displacement_set.shape[0]
    A = structure.lattice.matrix.T

    list_HNF = generate_all_superlattices(index)
    rotations, translations = get_symmetry_operations(structure)
    list_reduced_HNF = \
        reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)

    list_ds = []

    for hnf in list_reduced_HNF:
        labeling = Labeling(hnf, num_type,
                            num_site_parent, displacement_set,
                            rotations, translations)
        lbls_tmp = labeling.get_inequivalent_labelings()
        list_ds.extend([DerivativeStructure(hnf, num_type, A, lbl,
                                            num_site_parent, displacement_set)
                        for lbl in lbls_tmp])

    return list_ds


if __name__ == '__main__':
    from derivative_structure import get_lattice
    structure = get_lattice('fcc')
    index = 12
    num_type = 2

    list_ds = enumerate_derivative_structures(structure, index, num_type)
    print(3875, len(list_ds))
