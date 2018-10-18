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

    def point_to_factor(self, point):
        if not isinstance(point, np.ndarray):
            point = np.array(point)

        fct = np.mod(np.dot(self.left, point), self.snf.diagonal())
        fct = np.ravel_multi_index(fct, self.snf.diagonal())
        return fct

    def index_to_point(self, index):
        point = np.unravel_index(index, self.shape)
        return np.array(point, dtype=np.int)

    def point_to_index(self, point):
        index = np.ravel_multi_index(point, self.shape)
        return index

    def inverse_permutation(self, prm):
        inv = [0 for _ in range(self.num)]
        for i, e in enumerate(prm):
            inv[e] = i
        return inv

    def product_permutations(self, prm1, prm2):
        ret = [0 for _ in prm1]
        for i in range(self.num):
            ret[i] = prm2[prm1[i]]
        return ret

    def get_permutation_list(self):
        list_prm = []

        index_list = list(range(self.num))
        point_list = [self.index_to_point(index) for index in index_list]
        factor_list = [self.point_to_factor(point) for point in point_list]
        inv_factor_list = self.inverse_permutation(factor_list)

        for origin in point_list:
            displacement = origin
            translated_point_list = [point + displacement for point in point_list]
            translated_factor_list = [self.point_to_factor(point)
                                      for point in translated_point_list]

            prm = self.product_permutations(translated_factor_list, inv_factor_list)
            list_prm.append(prm)

        return list_prm


if __name__ == '__main__':
    hnf = np.array([
        [2, 0],
        [1, 4]
    ])

    permutation = Permutation(hnf)
    print(permutation.shape)
    print(permutation.snf)
    print(permutation.left)
    print(permutation.right)
    list_prm = permutation.get_permutation_list()
    print(len(list_prm))
    for prm in list_prm:
        print(prm)
