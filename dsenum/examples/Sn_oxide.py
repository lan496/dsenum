from pymatgen.io.cif import CifParser
from pymatgen.core import Specie, DummySpecie

from dsenum.enumerate import enumerate_derivatives


def get_structure(filename):
    cp = CifParser(filename)
    structure = cp.get_structures(primitive=True)[0]
    return structure


if __name__ == '__main__':
    rutile = get_structure('./TiO2_mp-2657_conventional_standard.cif')
    index = 2
    mapping_color_species = [DummySpecie('X', 0), Specie('O', -2), Specie('Sn', +2), Specie('Sn', +4)]
    num_type = len(mapping_color_species)
    base_site_constraints = [[2, 3],  # 2a
                             [2, 3],  # 2a
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1],  # 4f
                             [0, 1]]  # 4f

    list_dstructs = enumerate_derivatives(rutile, index, num_type, mapping_color_species,
                                          base_site_constraints=base_site_constraints,
                                          color_exchange=False,
                                          leave_superperiodic=False,
                                          use_all_colors=True)
    # import pdb; pdb.set_trace()
