from time import time

from tqdm import tqdm
from pymatgen.core.periodic_table import DummySpecie

from dsenum.labeling import Labeling, LabelGenerator, ListBasedLabelGenerator
from dsenum.derivative_structure import DerivativeStructure
from dsenum.utils import get_symmetry_operations

from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.coloring_generator import ColoringGenerator
from dsenum.coloring import SiteColoringEnumerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.derivative_structure import ColoringToStructure


def enumerate_derivative_structures(structure, index, num_type,
                                    ignore_site_property=False, leave_superperiodic=False,
                                    constraints=None, oxi_states=None,
                                    n_jobs=-1):
    displacement_set = structure.frac_coords
    num_site_parent = displacement_set.shape[0]
    A = structure.lattice.matrix.T

    list_reduced_HNF, rotations, translations = \
        generate_symmetry_distinct_superlattices(index, structure, return_symops=True)

    labelgen = LabelGenerator(index, num_type, num_site_parent, constraints, oxi_states,
                              n_jobs=n_jobs)

    list_ds = []

    for hnf in tqdm(list_reduced_HNF):
        labeling = Labeling(hnf, num_type, labelgen,
                            num_site_parent, displacement_set,
                            rotations, translations,
                            ignore_site_property=ignore_site_property,
                            leave_superperiodic=leave_superperiodic)
        lbls_tmp = labeling.get_inequivalent_labelings()
        print("HNF: {}".format(hnf.tolist()))
        list_ds.extend([DerivativeStructure(hnf, num_type, A, lbl,
                                            num_site_parent, displacement_set)
                        for lbl in lbls_tmp])
    print('total: {}'.format(len(list_ds)))

    return list_ds


def enumerate_derivatives(base_structure, index, num_type,
                          mapping_color_species=None,
                          color_exchange=True, leave_superperiodic=False, use_all_colors=True):
    """
    Parameter
    ---------
    base_structure: Structure
    index: int
    num_type: int
    mapping_color_species: if specified, use these species in derivative structures
    color_exchange: identify color-exchanging
    leave_superperiodic: do not discard superperiodic coloring

    Returns
    -------
    list_ds: list of derivative structure
    """
    start = time()

    displacement_set = base_structure.frac_coords
    list_reduced_HNF, rotations, translations = \
        generate_symmetry_distinct_superlattices(index, base_structure, return_symops=True)

    num_sites_base = base_structure.num_sites
    num_sites = num_sites_base * index

    cl_generator = ColoringGenerator(num_sites, num_type)

    if mapping_color_species and len(mapping_color_species) != num_type:
        raise ValueError('mapping_color_species must have num_type species.')
    if mapping_color_species is None:
        mapping_color_species = [DummySpecie(str(i)) for i in range(1, num_type + 1)]

    list_ds = []
    for hnf in tqdm(list_reduced_HNF):
        # print("HNF: {}".format(hnf.tolist()))
        ds_permutaion = DerivativeStructurePermutation(hnf, displacement_set,
                                                       rotations, translations)
        sc_enum = SiteColoringEnumerator(num_type, ds_permutaion, cl_generator,
                                         color_exchange, leave_superperiodic, use_all_colors)
        colorings = sc_enum.unique_colorings()

        # convert to Structure object
        cts = ColoringToStructure(base_structure, ds_permutaion.dhash, mapping_color_species)
        list_ds.extend([cts.convert_to_structure(cl) for cl in colorings])

    end = time()
    print('total: {} (Time: {:.4}sec)'.format(len(list_ds), end - start))

    return list_ds


def remove_symmetry_duplicates(structure, hnf, num_type, list_labelings):
    displacement_set = structure.frac_coords
    num_site_parent = displacement_set.shape[0]
    rotations, translations = get_symmetry_operations(structure)
    labelgen = ListBasedLabelGenerator(list_labelings)

    # discard superperiodic conf, and leave label-exchange duplicates
    labeling = Labeling(hnf, num_type, labelgen,
                        num_site_parent, displacement_set,
                        rotations, translations,
                        ignore_site_property=True, leave_superperiodic=False)

    lbls = labeling.get_inequivalent_labelings()
    return lbls


if __name__ == '__main__':
    from utils import get_lattice
    structure = get_lattice('fcc')
    index = 10
    num_type = 3

    list_ds = enumerate_derivatives(structure, index, num_type,
                                    color_exchange=False,
                                    leave_superperiodic=False,
                                    use_all_colors=True)
