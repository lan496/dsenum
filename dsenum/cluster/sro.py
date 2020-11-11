from typing import List, Tuple, Union

import numpy as np
from pymatgen.core import Structure
from tqdm import tqdm

from pyzdd import Universe
from pyzdd.graph import convert_to_raw_graph, VertexGraphFrontierManager
from pyzdd.structure import (
    construct_binary_derivative_structures_with_sro,
    enumerate_labelings_with_graph,
)

from dsenum.permutation_group import DerivativeMultiLatticeHash
from dsenum.cluster.point_cluster import PointCluster, EquivalentPointClusterGenerator
from dsenum.cluster.cluster_graph import BinaryPairClusterGraph
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.utils import gcd_list


class SROStructureEnumerator:
    """
    Parameters
    ----------
    base_structure: pymatgen.core.Structure
        Aristotype for derivative structures
    index: int
        How many times to expand unit cell
    num_types: int
        The number of species in derivative structures.
        `num_types` may be larger than the number of the kinds of species in `base_structure`: for example, you consider vacancies in derivative structures.
    composition_ratios: List[List[int]]
        composition_ratios[site_index][i] is a ratio of specie-`i` at site `site_index`
    distinct_cluster:
    """

    def __init__(
        self,
        base_structure: Structure,
        index: int,
        num_types: int,
        composition_ratios: List[List[int]],
        distinct_cluster: PointCluster,
    ):
        self._base_structure = base_structure
        self._index = index

        if num_types > 2:
            raise NotImplementedError
        self._num_types = num_types

        # enumerate sublattices
        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, base_structure, return_symops=True
        )
        self._list_reduced_HNF = list_reduced_HNF
        self._rotations = rotations
        self._translations = translations

        # convert composition constraints for supercell
        self._composition_ratios = composition_ratios
        supercell_composition_constraints, possible = convert_binary_composition_ratios(
            composition_ratios, self._index, self._num_types
        )
        self._supercell_composition_constraints = supercell_composition_constraints
        self._possible = possible

        self._distinct_cluster = distinct_cluster

    @property
    def base_structure(self):
        return self._base_structure

    @property
    def num_base_sites(self):
        return self.base_structure.num_sites

    @property
    def num_sites(self):
        return self.num_base_sites * self._index

    def generate(self) -> List[Tuple[np.ndarray, List[List[int]]]]:
        all_labelings_and_transformations = []
        for transformation in tqdm(self._list_reduced_HNF):
            list_labelings: List[List[int]] = []
            if self._possible:
                list_labelings = self.generate_with_hnf(transformation, only_count=False)  # type: ignore

            all_labelings_and_transformations.append((transformation, list_labelings))

        return all_labelings_and_transformations

    def count(self) -> List[Tuple[np.ndarray, int]]:
        all_counts_and_transformations = []
        for transformation in tqdm(self._list_reduced_HNF):
            count: int = 0
            if self._possible:
                count = self.generate_with_hnf(transformation, only_count=True)  # type: ignore

            all_counts_and_transformations.append((transformation, count))

        return all_counts_and_transformations

    def generate_with_hnf(
        self, transformation: np.ndarray, only_count=False
    ) -> Union[List[List[int]], int]:
        """
        return list of labelings with a fixed SRO and concentration
        """
        displacement_set = self.base_structure.frac_coords
        converter = DerivativeMultiLatticeHash(transformation, displacement_set)

        epcg = EquivalentPointClusterGenerator(
            displacement_set, self._rotations, self._translations, converter
        )
        pair_clusters = epcg.find_equivalent_point_clusters(self._distinct_cluster)
        bpcg = BinaryPairClusterGraph(converter, pair_clusters, self._composition_ratios)

        # target value for SQS
        target = bpcg.get_sqs_target_value()
        if not np.isclose(np.around(target), target):
            # if `target` is not integer, impossible to satisfy SRO constraint
            return []
        target = int(np.around(target))

        raw_graph, _ = convert_to_raw_graph(bpcg.graph)
        vgfm = VertexGraphFrontierManager(raw_graph)
        print("Max frontier size:", vgfm.get_max_frontier_size())

        dd = Universe()
        construct_binary_derivative_structures_with_sro(
            dd,
            self.num_sites,
            self._num_types,
            self._supercell_composition_constraints,
            vgfm,
            target,
        )
        print("Cardinality:", dd.cardinality())

        if only_count:
            return int(dd.cardinality())
        else:
            labelings = list(enumerate_labelings_with_graph(dd, self._num_types, raw_graph))
            return labelings


def convert_binary_composition_ratios(
    composition_ratios: List[List[int]], index: int, num_types: int
) -> Tuple[List[Tuple[List[int], int]], bool]:
    """
    Returns
    -------
    composition_ratios_supercell
    possible:
        if false, a given composition is not suited for the index.
    """
    if num_types > 2:
        raise NotImplementedError

    composition_ratios_supercell = []
    possible = True

    num_base_sites = len(composition_ratios)
    for site_index in range(num_base_sites):
        assert len(composition_ratios[site_index]) == 2
        # consistent with the order (site_index, *jimage)
        group = [site_index * index + i for i in range(index)]

        # number of label=1
        ratio_sum = np.sum(composition_ratios[site_index])
        ratio_gcd = gcd_list(composition_ratios[site_index])
        if index % (ratio_sum // ratio_gcd) != 0:
            # unable to satisfy composition ratio
            possible = False
        num_selected_atoms = int(np.around(index * composition_ratios[site_index][1] / ratio_sum))

        composition_ratios_supercell.append((group, num_selected_atoms))

    return (composition_ratios_supercell, possible)
