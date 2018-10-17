import numpy as np

from smith_normal_form import smith_normal_form


class Permutation(object):

    def __init__(self, hnf):
        self.hnf = hnf
        self.shape = tuple(hnf.diagonal())
        self.num = np.prod(self.shape)

        D, L, R = smith_normal_form(self.hnf)
        self.snf = D
        self.left = L
        self.right = R

    def factor(self, point):
        if not isinstance(point, np.ndarray):
            point = np.array(point)

        fct = np.mod(np.dot(self.left, point), self.D.diagonal())
        fct = np.ravel_multi_index(fct, self.snf.diagonal())
        return fct

    def get_factor_list(self, list_):
        pass

    def index_to_point(self, index):
        point = np.unravel_index(index, self.shape)
        return point

    def get_permutation_list(self):
        list_prm = []

        index_list = list(range(self.num))
        point_list = [self.index_to_point(e) for e in index_list]
        factor_list = [self.factor(e) for e in index_list]

        for direction in range(len(self.shape)):
            for step in range(1, self.shape[direction]):
                translated_index_list = np.array([

                ])
