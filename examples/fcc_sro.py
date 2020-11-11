import numpy as np
from pymatgen.core import Structure, Lattice

from dsenum.site import DerivativeSite
from dsenum.cluster.point_cluster import PointCluster
from dsenum.cluster.sro import SROStructureEnumerator


def get_fcc_structure(a: float):
    lattice = Lattice.from_parameters(a, a, a, 60, 60, 60)
    species = ["Cu"]
    frac_coords = np.array([[0, 0, 0]])
    structure = Structure(lattice, species, frac_coords)
    return structure


if __name__ == "__main__":
    """
    enumerate FCC derivative structures with the same 1st-NN-SRO of SQS-AB
    """
    aristo = get_fcc_structure(a=2.856)
    num_types = 2
    composition_ratios = [[1, 1]]
    # 1st NN
    cluster = PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))])

    index = 32

    sse = SROStructureEnumerator(aristo, index, num_types, composition_ratios, cluster)
    labelings_with_transformations = sse.count()
