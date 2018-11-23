from tqdm import tqdm

from superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry
)
from labeling import Labeling, ConstraintedLabeling
from derivative_structure import DerivativeStructure
from utils import get_symmetry_operations


def enumerate_derivative_structures(structure, index, num_type, constraints=None):
    displacement_set = structure.frac_coords
    num_site_parent = displacement_set.shape[0]
    A = structure.lattice.matrix.T

    list_HNF = generate_all_superlattices(index)
    rotations, translations = get_symmetry_operations(structure)
    list_reduced_HNF = \
        reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)

    list_ds = []

    for hnf in tqdm(list_reduced_HNF):
        if constraints is None:
            labeling = Labeling(hnf, num_type,
                                num_site_parent, displacement_set,
                                rotations, translations)
        else:
            labeling = ConstraintedLabeling(hnf, num_type,
                                            num_site_parent, displacement_set,
                                            rotations, translations,
                                            constraints=constraints)
        lbls_tmp = labeling.get_inequivalent_labelings()
        print(len(lbls_tmp))
        list_ds.extend([DerivativeStructure(hnf, num_type, A, lbl,
                                            num_site_parent, displacement_set)
                        for lbl in lbls_tmp])

    return list_ds


if __name__ == '__main__':
    from utils import get_fcc_with_vacancy
    structure = get_fcc_with_vacancy()
    index = 4
    num_type = 4
    constraints = [
        [0, 1],     # void and anion
        [0, 2, 3],  # void, cation1, and cation2
        [0, 2, 3],
        [0, 2, 3],
    ]

    list_ds = enumerate_derivative_structures(structure, index, num_type, constraints)
    print(list_ds[0].get_structure())
    print(len(list_ds))
