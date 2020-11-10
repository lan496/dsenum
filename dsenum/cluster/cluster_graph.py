from typing import List, Tuple
from typing import Counter as CounterType
from collections import Counter

import networkx as nx
from networkx.classes.function import number_of_selfloops
import numpy as np

from dsenum.site import DerivativeSite
from dsenum.converter import DerivativeMultiLatticeHash
from dsenum.cluster.point_cluster import PointCluster


class BinaryPairClusterGraph:
    """
    pair-cluster graph for binary system

    Parameters
    ----------
    converter:
    point_clusters: List[PointCluster]
    composition_ratio: List[List[int]]
        composition_ratio[i] is ratios of the i-th site in the primitive cell.
        composition_ratio[i][0] is for label=0, and composition_ratio[i][1] is for label=1.

    Attributes
    ----------
    graph: nx.Graph
        "weight" attribute remembers the multiplicity of edges
    ratio_sum: List[int]
    loop_offset: float
    """

    def __init__(
        self,
        converter: DerivativeMultiLatticeHash,
        point_clusters: List[PointCluster],
        composition_ratio: List[List[int]],
    ):
        self.converter = converter

        assert len(point_clusters[0]) == 2
        self.point_clusters = point_clusters

        assert len(composition_ratio) == self.converter.num_base_sites
        for ratio in composition_ratio:
            assert len(ratio) == 2
        self.ratio_sum = np.sum(composition_ratio, axis=1)
        self.composition_ratio = composition_ratio

        # vertices
        vertices = [
            int(self.converter.ravel_canonical_site(csite))
            for csite in converter.get_canonical_sites_list()
        ]

        # multiplicity and edges
        edges_counter: CounterType[Tuple[int, ...]] = Counter()
        for pcl in self.point_clusters:
            # points in a cluster are sorted by hashed index
            pair = sorted([self.converter.ravel_derivative_site(point) for point in pcl.points])
            edges_counter[tuple(pair)] += 1

        # number of pair clusters
        self._weight_sum = np.sum([count for _, count in edges_counter.items()])

        # self-loop offset
        # here, edges_counter may have self-loops
        self.loop_offset = 0
        for site_index in range(self.num_base_sites):
            src: int = self.converter.ravel_derivative_site(
                DerivativeSite(site_index, (0,) * self.converter.dim)
            )
            self.loop_offset += (
                self.composition_ratio[site_index][1]
                / self.ratio_sum[site_index]
                * edges_counter[(src, src)]
            )
        self.loop_offset *= self.index

        self._graph = nx.Graph()
        self._graph.add_nodes_from(vertices)
        self._graph.add_edges_from(
            [
                (src, dst, {"weight": weight})
                for (src, dst), weight in edges_counter.items()
                if src < dst
            ]
        )
        # now, graph is simple
        assert number_of_selfloops(self._graph) == 0

    @property
    def index(self):
        return self.converter.index

    @property
    def num_base_sites(self):
        return len(self.composition_ratio)

    @property
    def graph(self) -> nx.Graph:
        return self._graph

    def get_sqs_target_value(self) -> float:
        site_index = self.point_clusters[0].points[0].site_index
        c0 = self.composition_ratio[site_index][1] / self.ratio_sum[site_index]
        return self._weight_sum * c0 - self.loop_offset

    def calc_correlation(self, labeling: List[int]) -> float:
        """
        pair correlation between label=1 and label=1.
        For only binary system.
        """
        corr = np.sum(
            [
                weight * labeling[src] * labeling[dst]
                for src, dst, weight in self.graph.edges(data="weight")
            ]
        )
        corr /= self._weight_sum
        return corr + self.loop_offset

    def calc_short_range_order(self, labeling: List[int]) -> float:
        """
        Warren-Cowley short-range order between label=1 and label=0
        """
        site_index0 = self.point_clusters[0].points[0].site_index
        site_index1 = self.point_clusters[0].points[1].site_index
        c0 = (
            self.composition_ratio[site_index0][1] / self.ratio_sum[site_index0]
        )
        c1 = (
            self.composition_ratio[site_index1][1] / self.ratio_sum[site_index1]
        )
        sro = 1.0 - (c0 - self.calc_correlation(labeling)) / (c0 * (1.0 - c1))
        return sro
