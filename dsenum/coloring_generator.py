from functools import reduce
from itertools import product


class ColoringGenerator(object):

    def __init__(self, num_elements, num_color):
        self.num_elements = num_elements
        self.num_color = num_color

    def generate_all_colorings(self):
        list_colorings = list(product(range(self.num_color), repeat=self.num_elements))
        flags = {self.hash_coloring(coloring) for coloring in list_colorings}
        return list_colorings, flags

    def hash_coloring(self, coloring):
        ret = reduce(lambda x, y: x * self.num_type + y, coloring)
        return ret
