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
    clusters = [
        # 1st NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 1))]),
        # 2nd NN
        PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, -1))]),
        # 3rd NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 1))]),
        # 4th NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 0, 2))]),
        # 5th NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (2, 1, -1))]),
        # # 6th NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (1, 1, 1))]),
        # # 7th NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (0, 1, 2))]),
        # # 8th NN
        # PointCluster([DerivativeSite(0, (0, 0, 0)), DerivativeSite(0, (2, 2, -2))]),
    ]

    index = 8

    sse = SROStructureEnumerator(
        aristo,
        index,
        num_types,
        composition_ratios,
        clusters,
        remove_superperiodic=True,
    )
    labelings_with_transformations = sse.count()
