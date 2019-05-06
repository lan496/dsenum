import os

from pymatgen.io.cif import CifParser
from pymatgen.core import Specie, DummySpecie

from dsenum.enumerate import enumerate_derivative_structures
from dsenum.utils import write_cif


def get_structure(filename):
    cp = CifParser(filename)
    structure = cp.get_structures(primitive=True)[0]
    return structure


if __name__ == '__main__':
    rutile = get_structure('examples/TiO2_mp-2657_conventional_standard.cif')
    index = 1
    # mapping_color_species = [DummySpecie('X', 0), Specie('O', -2), Specie('Sn', +2), Specie('Sn', +4)]
    mapping_color_species = [DummySpecie('X'), Specie('O'), Specie('Sn')]
    num_type = len(mapping_color_species)
    """
    base_site_constraints = [[2, 3],  # 2a
                             [2, 3],  # 2a
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1]]  # 4f
    """
    base_site_constraints = [[2, ],  # 2a
                             [2, ],  # 2a
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1]]  # 4f

    list_dstructs = enumerate_derivative_structures(rutile, index, num_type, mapping_color_species,
                                                    base_site_constraints=base_site_constraints,
                                                    color_exchange=False,
                                                    leave_superperiodic=False,
                                                    use_all_colors=False)

    name = 'SnOx_index={}'.format(index)
    os.makedirs(os.path.join('examples', name), exist_ok=True)
    for i, dstruct in enumerate(list_dstructs):
        # remove void
        dstruct.remove_species([mapping_color_species[0], ])

        filename = os.path.join('examples', name, "SnOx_index={}_{}.cif".format(index, i))
        write_cif(filename, dstruct, refine_cell=True)
