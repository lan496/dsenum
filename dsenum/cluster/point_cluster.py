from copy import deepcopy
from dataclasses import astuple, dataclass
from itertools import combinations
from typing import List, Tuple

import numpy as np
from pymatgen.core import Structure

from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.site import DerivativeSite
from dsenum.utils import cast_integer_matrix, get_symmetry_operations


@dataclass
class PointCluster:
    """
    cluster class represents the set of DerivativeSite
    Attributes
    ----------
    points: list of Point
    """

    points: List[DerivativeSite]

    def __str__(self):
        points_str = [f"site_index={p.site_index}, jimage={p.jimage}" for p in self.points]
        rep = "PointCluster[" + ";".join(points_str) + "]"
        return rep

    def __eq__(self, other):
        if not isinstance(other, PointCluster):
            return False

        return tuple(self.points) == tuple(other.points)

    def __hash__(self):
        return hash(tuple(self.points))

    def __len__(self):
        return len(self.points)


class EquivalentPointClusterGenerator:
    """
    find equivalent point clusters
    Parameters
    ----------
    frac_coords: np.ndarray, (None, dim)
    rotations: np.ndarray(None, dim, dim)
    translations: np.ndarray(None, dim, dim)
    converter: DerivativeMultiLatticeHash
        converter for periodic boundary condition
    """

    def __init__(
        self,
        frac_coords: np.ndarray,
        rotations: np.ndarray,
        translations: np.ndarray,
        converter: DerivativeMultiLatticeHash,
    ):
        self.frac_coords = frac_coords
        self.rotations = rotations
        self.translations = translations
        self.converter = converter

    @property
    def dim(self):
        return self.frac_coords.shape[1]

    def _get_dsite(self, frac_coord: np.ndarray) -> DerivativeSite:
        """
        return DerivativeSite from fractional coorinates

        Parameters
        ----------
        frac_coords: (dim, )

        Returns
        -------
        dsite: DerivativeSite
        """
        for site_index, fc in enumerate(self.frac_coords):
            jimage = cast_integer_matrix(frac_coord - fc)
            if np.allclose(fc + jimage, frac_coord):
                dsite = DerivativeSite(site_index, tuple(jimage.tolist()))
                return dsite
        raise ValueError(f"invalid fractional coordinates: {frac_coord}")

    def normalize_point_cluster(self, point_cluster: PointCluster) -> PointCluster:
        """
        normalize point-cluster such that,
            - jimage of cluster[0] is (0, 0, 0)
            - points are in ascending order
        """
        # sort sites by (site_index, jimage)
        tuple_points = sorted(astuple(p) for p in point_cluster.points)
        points = [DerivativeSite(site_index, jimage) for site_index, jimage in tuple_points]

        offset = np.array(points[0].jimage)
        shifed_new_points = [
            DerivativeSite(p.site_index, tuple((np.array(p.jimage, dtype=int) - offset).tolist()))
            for p in points
        ]
        return PointCluster(shifed_new_points)

    def operate_point_cluster(
        self, point_cluster: PointCluster, R: np.ndarray, tau: np.ndarray
    ) -> PointCluster:
        """
        operate (R, tau) to a point cluster
        """
        new_points = []
        for point in point_cluster.points:
            site_index, jimage = astuple(point)
            new_frac_coord = np.dot(R, self.frac_coords[site_index] + np.array(jimage)) + tau
            new_point = self._get_dsite(new_frac_coord)
            new_points.append(new_point)
        return PointCluster(new_points)

    def find_equivalent_point_clusters(self, point_cluster: PointCluster) -> List[PointCluster]:
        # (rotations, translations) should contain identity operation
        equiv_clusters = {
            self.normalize_point_cluster(self.operate_point_cluster(point_cluster, R, tau))
            for R, tau in zip(self.rotations, self.translations)
        }

        all_equiv_clusters = []
        # apply translation
        for cluster in equiv_clusters:
            for lp in self.converter.get_lattice_points():
                translated_cluster = self.operate_point_cluster(
                    cluster, np.eye(self.dim), np.array(lp)
                )
                all_equiv_clusters.append(translated_cluster)
        # to account for multiplicity, do not unique all_equiv_clusters
        return all_equiv_clusters

    def get_point_cluster_size(
        self, point_cluster: PointCluster, lattice_matrix: np.ndarray
    ) -> float:
        """
        cluster size is defined as a maximum distance between sites in a point cluster

        Parameters
        ----------
        lattice_matrix: (dim, dim)
            row-wise basis vectors
        """
        size = 0
        for s0, s1 in combinations(point_cluster.points, r=2):
            f0 = self.frac_coords[s0.site_index] + np.array(s0.jimage)
            f1 = self.frac_coords[s1.site_index] + np.array(s1.jimage)
            dist = np.linalg.norm(np.dot(lattice_matrix.T, f1 - f0))
            size = max(size, dist)
        return size

    def get_distinct_point_clusters(
        self, structure: Structure, cutoff: float, order=2, eps=1e-8
    ) -> List[PointCluster]:
        """
        return all distinct `order`-sites point-clusters

        Parameters
        ----------
        structure: Structure
        cutoff: float
            cutoff radius
        order: int

        Returns
        -------
        distinct_point_clusters: List[PointCluster]
            returned list of point clusters are sorted by cluster size in the ascending order.
        """
        # search neighbor points whose distnace from some site in the unit cell is less than cutoff.
        list_points: List[DerivativeSite] = []
        for site in structure:
            neighbors = structure.get_neighbors(site, cutoff)
            for s in neighbors:
                if s.nn_distance - cutoff > eps:
                    continue
                point = self._get_dsite(s.frac_coords)
                list_points.append(point)
        list_points = list(set(list_points))

        # Find distinct singlet clusters
        distinct_point_clusters: List[PointCluster] = []
        found_singlets = set()
        for site_index in range(len(self.frac_coords)):
            cluster = PointCluster([DerivativeSite(site_index, (0, 0, 0))])
            if cluster in found_singlets:
                continue
            distinct_point_clusters.append(cluster)
            for sym_cluster in self.find_equivalent_point_clusters(cluster):
                found_singlets.add(sym_cluster)

        # grow pairs to triplet, quadruple, ...
        for num_nn in range(2, order + 1):
            next_distinct_point_clusters: List[PointCluster] = []
            found = set()

            for added_point in list_points:
                for precluster in distinct_point_clusters:
                    points: List[DerivativeSite] = list(precluster.points) + [added_point]
                    point_cluster: PointCluster = self.normalize_point_cluster(
                        PointCluster(points)
                    )
                    # skip an overlapped cluster
                    if len(set(point_cluster.points)) != num_nn:
                        continue

                    # skip too large point cluster
                    cluster_size = self.get_point_cluster_size(
                        point_cluster, structure.lattice.matrix
                    )
                    if cluster_size - cutoff > eps:
                        continue

                    if point_cluster in found:
                        continue
                    next_distinct_point_clusters.append(point_cluster)

                    equiv_clusters = self.find_equivalent_point_clusters(point_cluster)
                    found.update(equiv_clusters)

            distinct_point_clusters = deepcopy(next_distinct_point_clusters)

        # sort by cluser size
        distinct_point_clusters.sort(
            key=lambda pcl: self.get_point_cluster_size(pcl, structure.lattice.matrix)
        )

        return distinct_point_clusters

    def get_all_clusters(
        self, structure: Structure, cutoff: float, order: int = 2, eps: float = 1e-8
    ) -> Tuple[List[List[PointCluster]], List[float]]:
        """
        return all `order`-sites clusters

        Parameters
        ----------
        structure: Structure
        cutoff: float
            cutoff radius
        order: int

        Returns
        -------
        grouped_point_clusters: List[List[PointCluster]]
            grouped_points_clusters[i] is a list of the i-th equivalent group of clusters
        point_cluster_sizes: List[float]
            point_cluster_sizes[i] is a size of the i-th equivalent group of clusters
        """
        # find distinct point clusters
        distinct_temp_point_clusters = self.get_distinct_point_clusters(
            structure, cutoff, order, eps
        )

        # generate symmetry-equivalent point clusters
        grouped_clusters_sizes = []
        for point_cluster in distinct_temp_point_clusters:
            cluster_size = self.get_point_cluster_size(point_cluster, structure.lattice.matrix)
            equivs = self.find_equivalent_point_clusters(point_cluster)
            grouped_clusters_sizes.append((equivs, cluster_size))

        # sort by cluster size
        grouped_clusters_sizes.sort(key=lambda e: e[1])

        grouped_point_clusters = [gcls for gcls, _ in grouped_clusters_sizes]
        point_cluster_sizes = [size for _, size in grouped_clusters_sizes]

        return grouped_point_clusters, point_cluster_sizes

    @classmethod
    def from_structure(
        cls, structure: Structure, transformation: np.ndarray, symprec: float = 1e-2
    ):
        """
        create EquivalentPointClusterGenerator from pymatgen.core.Structure

        Parameters
        ----------
        structure: pymatgen.core.Structure, base structure
        transformation: transformation matrix of sublattice
        symprec: float
        """
        frac_coords = structure.frac_coords
        rotations, translations = get_symmetry_operations(structure, symprec)
        converter = DerivativeMultiLatticeHash(transformation, frac_coords)
        return cls(
            frac_coords=frac_coords,
            rotations=rotations,
            translations=translations,
            converter=converter,
        )
