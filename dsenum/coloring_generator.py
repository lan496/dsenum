from functools import reduce
from itertools import product

from sympy.utilities.iterables import multiset_permutations


class BaseColoringGenerator:

    def generate_all_colorings(self):
        raise NotImplementedError

    def hash_coloring(self, coloring):
        raise NotImplementedError


def hash_in_all_configuration(coloring, num_color):
    ret = reduce(lambda x, y: x * num_color + y, coloring)
    return ret


class ColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_elements, num_color):
        self.num_elements = num_elements
        self.num_color = num_color

    def generate_all_colorings(self):
        list_colorings = list(product(range(self.num_color), repeat=self.num_elements))
        flags = {self.hash_coloring(coloring): True for coloring in list_colorings}
        return list_colorings, flags

    def hash_coloring(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)


class ListBasedColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_color, list_colorings):
        self.num_color = num_color
        self.list_colorings = list_colorings

    def generate_all_colorings(self):
        flags = {self.hash_coloring(coloring): True for coloring in self.list_colorings}
        return self.list_colorings, flags

    def hash_coloring(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)


class FixedConcentrationColoringGenerator(BaseColoringGenerator):
    """
    Parameters
    ----------
    num_elements: int
    num_color: int
    color_ratio: list of int
        num_elements % sum(color_ratio) == 0
    """

    def __init__(self, num_elements, num_color, color_ratio):
        self.num_elements = num_elements
        self.num_color = num_color
        self.color_ratio = color_ratio

        if num_elements % sum(self.color_ratio) != 0:
            raise ValueError('incorrect composition ratio')

        factor = num_elements // sum(self.color_ratio)
        self.num_elements_each_color = [factor * cr for cr in self.color_ratio]

    def generate_all_colorings(self):
        first_coloring = []
        for i in range(self.num_color):
            first_coloring.extend([i for _ in range(self.num_elements_each_color[i])])

        # TODO: inefficient to cast to list
        list_colorings = list(multiset_permutations(first_coloring))
        flags = {self.hash_coloring(coloring): True for coloring in list_colorings}
        return list_colorings, flags

    def hash_coloring(self, coloring):
        return hash_in_all_configuration(coloring, self.num_color)


class CompositionColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_sites_base, num_color, mapping_color_speices, site_constraints=None):
        # TODO
        pass
