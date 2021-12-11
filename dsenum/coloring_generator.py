from abc import ABCMeta, abstractmethod
from itertools import product
from typing import List

import numpy as np
from sympy.utilities.iterables import multiset_permutations

from dsenum.core import get_composition, hash_in_all_configuration  # type: ignore


class BaseColoringGenerator(metaclass=ABCMeta):
    @abstractmethod
    def generate_all_colorings(self):
        raise NotImplementedError

    @abstractmethod
    def yield_coloring(self):
        raise NotImplementedError


class ColoringGenerator(BaseColoringGenerator):
    def __init__(self, num_elements: int, num_color: int, site_constraints=None):
        """
        Parameters
        ----------
        num_elements:
            the number of elements to color
        num_color:
            the number of color
        site_constraints: list of list of int, optional
            the length of this list should be equal to `num_elements`
        """
        self.num_elements = num_elements
        self.num_color = num_color

        if site_constraints is not None:
            assert len(site_constraints) == self.num_elements
        self.site_constraints = site_constraints

    def generate_all_colorings(self):
        if self.site_constraints:
            list_colorings = []
            flags = dict()
            for cl_compressed in product(*[range(len(sc)) for sc in self.site_constraints]):
                cl = [self.site_constraints[i][idx] for i, idx in enumerate(cl_compressed)]
                list_colorings.append(cl)
                flags[hash_in_all_configuration(cl, self.num_color)] = True
        else:
            list_colorings = list(product(range(self.num_color), repeat=self.num_elements))
            flags = {
                hash_in_all_configuration(list(coloring), self.num_color): True
                for coloring in list_colorings
            }

        return list_colorings, flags

    def yield_coloring(self):
        if self.site_constraints:
            for cl_compressed in product(*[range(len(sc)) for sc in self.site_constraints]):
                cl = [self.site_constraints[i][idx] for i, idx in enumerate(cl_compressed)]
                yield cl
        else:
            for cl in product(range(self.num_color), repeat=self.num_elements):
                yield list(cl)


class ListBasedColoringGenerator(BaseColoringGenerator):
    """
    Parameters
    ----------
    num_color: int
    list_colorings: list of generator
    """

    def __init__(self, num_color, list_colorings):
        self.num_color = num_color
        self.list_colorings = list_colorings

    def generate_all_colorings(self):
        flags = {
            hash_in_all_configuration(coloring, self.num_color): True
            for coloring in self.list_colorings
        }
        return self.list_colorings, flags

    def yield_coloring(self):
        for cl in self.list_colorings:
            yield cl


class FixedConcentrationColoringGenerator(BaseColoringGenerator):
    """
    Parameters
    ----------
    num_elements: int
    num_color: int
    color_ratio: list of int
        length of `color_ratio` should be equal to `num_color`.
        num_elements % sum(color_ratio) == 0
    site_constraints: (Optional), list (num_elements, num_color)
        e.g. site_constraints[2] = [0, 3, 4] means color of site-2 must be 0, 3, or 4.
    """

    def __init__(
        self, num_elements: int, num_color: int, color_ratio: List[float], site_constraints=None
    ):
        self.num_elements = num_elements
        self.num_color = num_color

        self.site_constraints = site_constraints

        self.color_ratio = color_ratio
        ratio_sum = int(np.around(sum(self.color_ratio)))
        if num_elements % ratio_sum != 0:
            raise ValueError("incorrect composition ratio")

        factor = num_elements // ratio_sum
        self.num_elements_each_color = [int(np.around(factor * cr)) for cr in self.color_ratio]

    def generate_all_colorings(self):
        if self.site_constraints:
            list_colorings = []
            flags = dict()
            # Apply site_constraints first
            for cl_compressed in product(*[range(len(sc)) for sc in self.site_constraints]):
                cl = [self.site_constraints[i][idx] for i, idx in enumerate(cl_compressed)]
                if get_composition(cl, self.num_color) == self.num_elements_each_color:
                    list_colorings.append(cl)
                    flags[hash_in_all_configuration(cl, self.num_color)] = True
        else:
            first_coloring = []
            for i in range(self.num_color):
                first_coloring.extend([i for _ in range(self.num_elements_each_color[i])])

            # TODO: inefficient to cast to list
            list_colorings = list(multiset_permutations(first_coloring))
            flags = {
                hash_in_all_configuration(coloring, self.num_color): True
                for coloring in list_colorings
            }
        return list_colorings, flags

    def yield_coloring(self):
        # create one of colorings to use multiset_permutations
        first_coloring = []
        for i in range(self.num_color):
            first_coloring.extend([i for _ in range(self.num_elements_each_color[i])])

        if self.site_constraints:
            for cl in multiset_permutations(first_coloring):
                # TODO: more efficient way to adopt site_constraints
                if satisfy_site_constraints(self.site_constraints, cl):
                    yield cl
        else:
            for cl in multiset_permutations(first_coloring):
                yield cl


def satisfy_site_constraints(site_constraints, coloring):
    if all([(color in site_constraints[i]) for i, color in enumerate(coloring)]):
        return True
    else:
        return False
