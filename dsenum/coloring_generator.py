from functools import reduce
from itertools import product


class BaseColoringGenerator:

    def generate_all_colorings(self):
        raise NotImplementedError

    def hash_coloring(self, coloring):
        raise NotImplementedError


class ColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_elements, num_color):
        self.num_elements = num_elements
        self.num_color = num_color

    def generate_all_colorings(self):
        list_colorings = list(product(range(self.num_color), repeat=self.num_elements))
        flags = {self.hash_coloring(coloring): True for coloring in list_colorings}
        return list_colorings, flags

    def hash_coloring(self, coloring):
        ret = reduce(lambda x, y: x * self.num_color + y, coloring)
        return ret


class ListBasedColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_color, list_colorings):
        self.num_color = num_color
        self.list_colorings = list_colorings

    def generate_all_colorings(self):
        flags = {self.hash_coloring(coloring): True for coloring in self.list_colorings}
        return self.list_colorings, flags

    def hash_coloring(self, coloring):
        ret = reduce(lambda x, y: x * self.num_color + y, coloring)
        return ret


class FixedConcentrationColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_elements, num_color, color_ratio):
        self.num_elements = num_elements
        self.num_color = num_color
        self.color_ratio = color_ratio

        # TODO


class CompositionColoringGenerator(BaseColoringGenerator):

    def __init__(self, num_sites_base, num_color, mapping_color_speices, site_constraints=None):
        # TODO
        pass
