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
        fct = np.mod(np.dot(self.left, point), self.D.diagonal())
        return fct

    def get_factor_list(self, list_):
        pass

    def index_to_point(self, index):
        pass

    def get_permutation_list(self):
        list_prm = []

        e = np.arange(self.num)
        list_prm.append(e)
