from time import time

from tqdm import tqdm
from pymatgen.core.periodic_table import DummySpecie

from dsenum.utils import get_symmetry_operations
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.coloring_generator import (
    BaseColoringGenerator,
    ColoringGenerator,
    FixedConcentrationColoringGenerator,
    ListBasedColoringGenerator,
)
from dsenum.coloring import SiteColoringEnumerator
from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.derivative_structure import ColoringToStructure


def enumerate_derivative_structures(
    base_structure,
    index,
    num_type,
    mapping_color_species=None,
    composition_constraints=None,
    base_site_constraints=None,
    color_exchange=True,
    leave_superperiodic=False,
    use_all_colors=True,
    method="direct",
    n_jobs=1,
):
    """
    Parameter
    ---------
    base_structure: Structure
    index: int
    num_type: int
    mapping_color_species: (Optional) if specified, use these species in derivative structures
    composition_constraints: (Optional) None or list of int
    base_site_constraints: (Optional) list (num_elements, num_color)
        e.g. site_constraints[2] = [0, 3, 4] means color of site-2 in base_structure must be 0, 3, or 4.
    color_exchange: identify color-exchanging
    leave_superperiodic: do not discard superperiodic coloring
    use_all_colors: bool
    method: "direct" or "lexicographic", so far
    n_jobs: core in lexicographic coset enumeration(only used when method='lexicographic')

    Returns
    -------
    list_ds: list of derivative structure
    """
    start = time()

    list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
        index, base_structure, return_symops=True
    )

    num_sites_base = base_structure.num_sites
    num_sites = num_sites_base * index

    # site constraints
    if base_site_constraints:
        site_constraints = DerivativeMultiLatticeHash.convert_site_constraints(
            base_site_constraints, index
        )
    else:
        site_constraints = None

    # composition constraints
    if composition_constraints is None:
        cl_generator = ColoringGenerator(num_sites, num_type, site_constraints)
    else:
        cl_generator = FixedConcentrationColoringGenerator(
            num_sites, num_type, composition_constraints, site_constraints
        )

    if mapping_color_species and len(mapping_color_species) != num_type:
        raise ValueError("mapping_color_species must have num_type species.")
    if mapping_color_species is None:
        mapping_color_species = [DummySpecie(str(i)) for i in range(1, num_type + 1)]

    list_ds = []
    for hnf in tqdm(list_reduced_HNF):
        list_ds_hnf = enumerate_with_hnf(
            base_structure,
            hnf,
            num_type,
            rotations,
            translations,
            cl_generator,
            mapping_color_species,
            color_exchange,
            leave_superperiodic,
            use_all_colors,
            method,
            n_jobs,
        )
        list_ds.extend(list_ds_hnf)

    end = time()
    print("total: {} (Time: {:.4}sec)".format(len(list_ds), end - start))

    return list_ds


def enumerate_with_hnf(
    base_structure,
    hnf,
    num_type,
    rotations,
    translations,
    cl_generator: BaseColoringGenerator,
    mapping_color_species,
    color_exchange: bool,
    leave_superperiodic: bool,
    use_all_colors: bool,
    method="direct",
    n_jobs=1,
):
    displacement_set = base_structure.frac_coords
    ds_permutation = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
    sc_enum = SiteColoringEnumerator(
        num_type,
        ds_permutation,
        cl_generator,
        color_exchange,
        leave_superperiodic,
        use_all_colors,
        method=method,
        n_jobs=n_jobs,
    )
    colorings = sc_enum.unique_colorings()

    # convert to Structure object
    cts = ColoringToStructure(base_structure, ds_permutation.dhash, mapping_color_species)
    list_ds = [cts.convert_to_structure(cl) for cl in colorings]
    return list_ds


def remove_symmetry_duplicates(
    base_structure,
    hnf,
    num_type,
    list_colorings,
    color_exchange: bool,
    leave_superperiodic: bool,
    use_all_colors: bool,
):
    displacement_set = base_structure.frac_coords
    rotations, translations = get_symmetry_operations(base_structure)
    cl_generator = ListBasedColoringGenerator(num_type, list_colorings)

    ds_permutation = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
    sc_enum = SiteColoringEnumerator(
        num_type, ds_permutation, cl_generator, color_exchange, leave_superperiodic, use_all_colors
    )
    colorings = sc_enum.unique_colorings()
    return colorings


def remove_symmetry_duplicates_from_generator(
    base_structure,
    hnf,
    num_type,
    list_colorings,
    color_exchange: bool,
    leave_superperiodic: bool,
    use_all_colors: bool,
    method="direct",
    n_jobs=1,
):
    displacement_set = base_structure.frac_coords
    rotations, translations = get_symmetry_operations(base_structure)
    cl_generator = ListBasedColoringGenerator(num_type, list_colorings)

    ds_permutation = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
    sc_enum = SiteColoringEnumerator(
        num_type,
        ds_permutation,
        cl_generator,
        color_exchange,
        leave_superperiodic,
        use_all_colors,
        method=method,
        n_jobs=n_jobs,
    )

    colorings = sc_enum.unique_colorings()
    return colorings


if __name__ == "__main__":
    from utils import get_lattice

    structure = get_lattice("fcc")
    index = 15
    num_type = 2

    list_ds = enumerate_derivative_structures(
        structure,
        index,
        num_type,
        color_exchange=True,
        leave_superperiodic=False,
        use_all_colors=True,
        method="direct",
    )
