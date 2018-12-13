from tqdm import tqdm

from derivative.superlattice import (
    generate_all_superlattices,
    reduce_HNF_list_by_parent_lattice_symmetry
)
from derivative.labeling import Labeling, LabelGenerator
from derivative.derivative_structure import DerivativeStructure
from derivative.utils import get_symmetry_operations


def enumerate_derivative_structures(structure, index, num_type,
                                    ignore_site_property=False,
                                    constraints=None, oxi_states=None,
                                    n_jobs=-1):
    displacement_set = structure.frac_coords
    num_site_parent = displacement_set.shape[0]
    A = structure.lattice.matrix.T

    list_HNF = generate_all_superlattices(index)
    rotations, translations = get_symmetry_operations(structure)
    list_reduced_HNF = \
        reduce_HNF_list_by_parent_lattice_symmetry(list_HNF, rotations)

    labelgen = LabelGenerator(index, num_type, num_site_parent, constraints, oxi_states,
                              n_jobs=n_jobs)

    list_ds = []

    for hnf in tqdm(list_reduced_HNF):
        labeling = Labeling(hnf, num_type, labelgen,
                            num_site_parent, displacement_set,
                            rotations, translations,
                            ignore_site_property=ignore_site_property)
        lbls_tmp = labeling.get_inequivalent_labelings()
        print("HNF: {}".format(hnf.tolist()))
        print(len(lbls_tmp))
        list_ds.extend([DerivativeStructure(hnf, num_type, A, lbl,
                                            num_site_parent, displacement_set)
                        for lbl in lbls_tmp])
    print('total: {}'.format(len(list_ds)))

    return list_ds


if __name__ == '__main__':
    from utils import get_fcc_with_vacancy
    structure = get_fcc_with_vacancy()
    index = 1
    num_type = 4
    constraints = [
        [0, 1],     # void and anion
        [0, 2, 3],  # void, cation1, and cation2
        [0, 2, 3],
        [0, 2, 3],
    ]

    list_ds = enumerate_derivative_structures(structure, index, num_type,
                                              constraints=constraints)
    print(list_ds[0].get_structure())
    print(len(list_ds))
