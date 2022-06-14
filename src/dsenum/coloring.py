from abc import ABCMeta, abstractmethod
from itertools import permutations
from multiprocessing import Pool, cpu_count
from typing import List, Tuple, cast

from dsenum.coloring_generator import BaseColoringGenerator
from dsenum.core import act_permutation, hash_in_all_configuration  # type: ignore
from dsenum.permutation_group import DerivativeStructurePermutation


class AbstractEnumerator(metaclass=ABCMeta):
    @abstractmethod
    def coset_enumerate(self) -> List[List[int]]:
        raise NotImplementedError


class DirectColoringEnumerator(AbstractEnumerator):
    """
    Parameters
    ----------
    permutation_group: list of permutation
    num_color: int
    cl_generator: BaseColoringGenerator
    color_exchange: bool
    """

    def __init__(
        self,
        permutation_group: List[List[int]],
        num_color: int,
        cl_generator: BaseColoringGenerator,
        color_exchange: bool = True,
    ):
        self.permutation_group = permutation_group
        self.num_color = num_color
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange

        self.list_colorings, self.flags = self.cl_generator.generate_all_colorings()

    def _hash(self, coloring: List[int]) -> int:
        return hash_in_all_configuration(coloring, self.num_color)

    def _walk_orbit(self, coloring: List[int], include_identity: bool = False) -> None:
        offset = 0 if include_identity else 1
        # assume self.permutation_group[0] is identity
        for prm in self.permutation_group[offset:]:
            acted_cl = act_permutation(prm, coloring)
            acted_cl_hash = self._hash(acted_cl)
            self.flags[acted_cl_hash] = False

    def coset_enumerate(self) -> List[List[int]]:
        colorings = []

        for cl in self.list_colorings:
            cl_hash = self._hash(cl)
            # avoid already-visited coloring
            if not self.flags[cl_hash]:
                continue
            colorings.append(cl)

            self._walk_orbit(cl, include_identity=False)

            if self.color_exchange:
                for cl_prm in permutations(range(self.num_color)):
                    exchanged_cl = [cl_prm[c] for c in cl]
                    exchanged_cl_hash = self._hash(exchanged_cl)
                    if exchanged_cl_hash == cl_hash:
                        continue
                    self._walk_orbit(exchanged_cl, include_identity=True)

        return colorings


class LexicographicColoringEnumerator(AbstractEnumerator):
    """
    this algorithm takes the most lexicographically small colorings

    Parameters
    ----------
    permutation_group: list of permutation
    num_color: int
    cl_generator: BaseColoringGenerator
    color_exchange: bool
    n_jobs: int
        when n_jobs > 1, use multiprocessing but it seems slower than signle core......
    """

    def __init__(
        self,
        permutation_group: List[List[int]],
        num_color: int,
        cl_generator: BaseColoringGenerator,
        color_exchange: bool = True,
        n_jobs: int = 1,
    ):
        self.permutation_group = permutation_group
        self.num_color = num_color
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        if n_jobs == -1:
            self.n_jobs = cpu_count()
        else:
            self.n_jobs = n_jobs

    def _hash(self, coloring: List[int]) -> int:
        return hash_in_all_configuration(coloring, self.num_color)

    def _is_champion_coloring(self, coloring: List[int]) -> bool:
        cl_hash = self._hash(coloring)

        if self.color_exchange:
            for cl_prm in permutations(range(self.num_color)):
                exchanged_cl = [cl_prm[c] for c in coloring]
                for prm in self.permutation_group:
                    acted_cl = act_permutation(prm, exchanged_cl)
                    acted_cl_hash = self._hash(acted_cl)
                    if acted_cl_hash < cl_hash:
                        return False
        else:
            for prm in self.permutation_group:
                acted_cl = act_permutation(prm, coloring)
                acted_cl_hash = self._hash(acted_cl)
                if acted_cl_hash < cl_hash:
                    return False
        return True

    def _is_champion_coloring_parallel(self, coloring: List[int]) -> Tuple[List[int], bool]:
        return coloring, self._is_champion_coloring(coloring)

    def coset_enumerate(self) -> List[List[int]]:
        if self.n_jobs != 1:
            with Pool(self.n_jobs) as pool:
                colorings = [
                    cl
                    for cl, keep in pool.map(
                        self._is_champion_coloring_parallel, self.cl_generator.yield_coloring()
                    )
                    if keep
                ]
        else:
            colorings = []
            for cl in self.cl_generator.yield_coloring():
                if self._is_champion_coloring(cl):
                    colorings.append(cl)
        return colorings


class SiteColoringEnumerator(object):
    """
    Parameters
    ----------
    ds_permutation: DerivativeStructurePermutation object
    num_color: int
    color_exchange: bool
    remove_superperiodic: bool
    method: "direct" or "lexicographic", so far
    n_jobs: int, use only when method is "lexicographic"
    """

    def __init__(
        self,
        num_color: int,
        ds_permutation: DerivativeStructurePermutation,
        cl_generator: BaseColoringGenerator,
        color_exchange: bool = True,
        remove_superperiodic: bool = True,
        remove_incomplete: bool = True,
        method: str = "direct",
        n_jobs: int = 1,
    ):
        self.num_color = num_color
        self.ds_permutation = ds_permutation
        self.cl_generator = cl_generator
        self.color_exchange = color_exchange
        self.remove_superperiodic = remove_superperiodic
        self.remove_incomplete = remove_incomplete
        self.method = method
        self.n_jobs = n_jobs

        self.permutation_group = self.ds_permutation.get_symmetry_operation_permutations()

        # typing.cast causes no runtime effect
        if self.method == "direct":
            self.clenum = cast(
                AbstractEnumerator,
                DirectColoringEnumerator(
                    self.permutation_group,
                    self.num_color,
                    self.cl_generator,
                    color_exchange=color_exchange,
                ),
            )
        elif self.method == "lexicographic":
            self.clenum = cast(
                AbstractEnumerator,
                LexicographicColoringEnumerator(
                    self.permutation_group,
                    self.num_color,
                    self.cl_generator,
                    color_exchange=color_exchange,
                    n_jobs=self.n_jobs,
                ),
            )
        else:
            raise ValueError("Unknown method: ", self.method)

    @property
    def translation_permutations(self) -> List[List[int]]:
        return self.ds_permutation.prm_t

    def _hash(self, coloring: List[int]) -> int:
        return hash_in_all_configuration(coloring, self.num_color)

    def unique_colorings(self) -> List[List[int]]:
        symmetry_uniqued_coloring = self.clenum.coset_enumerate()
        colorings = []

        for cl in symmetry_uniqued_coloring:
            if self.remove_incomplete and (not self._has_all_colors(cl)):
                continue
            if self.remove_superperiodic and self._is_superperiodic(cl):
                continue
            colorings.append(cl)

        return colorings

    def _has_all_colors(self, coloring: List[int]) -> bool:
        if len(set(coloring)) == self.num_color:
            return True
        else:
            return False

    def _is_superperiodic(self, coloring: List[int]) -> bool:
        # assume self.translation_permutations[0] is identity
        if any(
            [
                self._hash(act_permutation(prm, coloring)) == self._hash(coloring)
                for prm in self.translation_permutations[1:]
            ]
        ):
            return True
        else:
            return False
