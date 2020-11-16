from typing import List, Tuple, Union

import numpy as np
from pymatgen.core import Structure
from pymatgen.core.periodic_table import DummySpecie, Specie, Element
from tqdm import tqdm

from pyzdd import Universe, Permutation
from pyzdd.graph import (
    convert_to_raw_graph,
    VertexGraphFrontierManager,
    get_vertex_order_by_bfs,
)
from pyzdd.structure import (
    construct_derivative_structures_with_sro,
    enumerate_binary_labelings_with_graph,
)

from dsenum.permutation_group import DerivativeMultiLatticeHash, DerivativeStructurePermutation
from dsenum.cluster.point_cluster import PointCluster, EquivalentPointClusterGenerator
from dsenum.cluster.cluster_graph import BinaryPairClusterGraph
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.derivative_structure import ColoringToStructure
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
        distinct_clusters: List[PointCluster],
        mapping_color_species: List[Union[str, Element, Specie, DummySpecie]] = None,
        remove_superperiodic=False,
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

        assert len(distinct_clusters) > 0
        self._distinct_clusters = distinct_clusters
        self._remove_superperiodic = remove_superperiodic

        if mapping_color_species and len(mapping_color_species) != self.num_types:
            raise ValueError("mapping_color_species must have num_type species.")
        if mapping_color_species is None:
            mapping_color_species = [DummySpecie(str(i)) for i in range(1, self.num_types + 1)]
        self.mapping_color_species = mapping_color_species

    @property
    def base_structure(self):
        return self._base_structure

    @property
    def num_base_sites(self):
        return self.base_structure.num_sites

    @property
    def num_sites(self):
        return self.num_base_sites * self._index

    @property
    def num_types(self):
        return self._num_types

    def generate_structures(
        self,
        additional_species=None,
        additional_frac_coords=None,
    ) -> List[Structure]:
        displacement_set = self.base_structure.frac_coords
        list_dstructs: List[Structure] = []
        if not self._possible:
            return list_dstructs

        for transformation in tqdm(self._list_reduced_HNF):
            list_labelings = self.generate_with_hnf(transformation, only_count=False)  # type: ignore

            # convert to Structure object
            ds_permutation = DerivativeStructurePermutation(
                transformation,
                displacement_set,
                self._rotations,
                self._translations,
            )
            cts = ColoringToStructure(
                self.base_structure,
                ds_permutation.dhash,
                self.mapping_color_species,
                additional_species=additional_species,
                additional_frac_coords=additional_frac_coords,
            )
            list_ds_hnf = [cts.convert_to_structure(cl) for cl in list_labelings]  # type: ignore
            list_dstructs.extend(list_ds_hnf)

        return list_dstructs

    def generate_labelings(self) -> List[Tuple[np.ndarray, List[List[int]]]]:
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

        # site converter
        ds_permutation = DerivativeStructurePermutation(
            transformation,
            displacement_set,
            self._rotations,
            self._translations,
        )
        converter = ds_permutation.dhash

        # permutation group from symmetry operations
        automorphism = [
            Permutation(sigma) for sigma in ds_permutation.get_symmetry_operation_permutations()
        ]
        translation_group = [Permutation(sigma) for sigma in ds_permutation._prm_t]

        epcg = EquivalentPointClusterGenerator(
            displacement_set, self._rotations, self._translations, converter
        )

        # fix vertex order
        pair_clusters0 = epcg.find_equivalent_point_clusters(self._distinct_clusters[0])
        bpcg0 = BinaryPairClusterGraph(converter, pair_clusters0, self._composition_ratios)
        raw_graph0, _ = convert_to_raw_graph(bpcg0.graph)
        vertex_order = get_vertex_order_by_bfs(raw_graph0)

        graphs = []
        targets = []
        for cluster in self._distinct_clusters:
            pair_clusters = epcg.find_equivalent_point_clusters(cluster)
            bpcg = BinaryPairClusterGraph(converter, pair_clusters, self._composition_ratios)
            raw_graph, _ = convert_to_raw_graph(bpcg.graph)
            graphs.append(raw_graph)

            target = bpcg.get_sqs_target_value()
            # target value for SQS
            if not np.isclose(np.around(target), target):
                # if `target` is not integer, impossible to satisfy SRO constraint
                return []
            target = int(np.around(target))
            targets.append(target)

        dd = Universe()
        construct_derivative_structures_with_sro(
            dd,
            self.num_sites,
            self._num_types,
            vertex_order,
            automorphism,
            translation_group,
            self._supercell_composition_constraints,
            graphs,
            targets,
            self._remove_superperiodic,
        )
        print("Cardinality:", dd.cardinality())

        if only_count:
            return int(dd.cardinality())
        else:
            if self._num_types == 2:
                labelings = list(
                    enumerate_binary_labelings_with_graph(dd, self.num_sites, vertex_order)
                )
            else:
                raise NotImplementedError
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
