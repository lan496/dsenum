from functools import lru_cache
from itertools import product
from typing import List

from scipy.special import binom


def polya_counting(permutation_group: List[List[int]], num_color: int) -> int:
    """
    count the number of coloring with num_color kinds of colors and permutations

    Parameters
    ----------
    permutation_group:
        the i-th permutation permutation_group[i] permutes j to permutation_group[i][j]
        for j = 0,...,len(permutation_group[0])-1

    Returns
    -------
    cnt: int
    """
    cnt = 0
    for perm in permutation_group:
        type_of_perm = get_type_of_permutation(perm)
        cnt += num_color ** sum(type_of_perm)

    assert cnt % len(permutation_group) == 0
    cnt //= len(permutation_group)
    return cnt


@lru_cache(maxsize=None)
def get_inventory_coefficient(type_of_perm: tuple, didx, num_elements_of_each_color: tuple):
    # termination
    if didx == 0:
        if all([e == 0 for e in num_elements_of_each_color]):
            return 1
        else:
            return 0

    ret = 0

    for nec in product(*[range(0, e + 1, didx) for e in num_elements_of_each_color]):
        list_k = [e // didx for e in nec]
        if sum(list_k) != type_of_perm[didx - 1]:
            continue
        complement = tuple(e1 - e2 for e1, e2 in zip(num_elements_of_each_color, nec))
        ret += get_multinomial_coefficient(list_k) * get_inventory_coefficient(
            type_of_perm, didx - 1, complement
        )

    return ret


# ref: https://stackoverflow.com/questions/46374185/does-python-have-a-function-which-computes-multinomial-coefficients
def get_multinomial_coefficient(params):
    if len(params) == 1:
        return 1
    if params[-1] == 0:
        return get_multinomial_coefficient(params[:-1])
    return int(binom(sum(params), params[-1])) * get_multinomial_coefficient(params[:-1])


def polya_fixed_degrees_counting(
    permutation_group: List[List[int]], num_color: int, num_elements_of_each_color: List[int]
):
    """
    count the number of coloring with num_color kinds of colors and permutations.
    In addition, the number of each colors is fixed with num_elements_of_each_color.

    Parameters
    ----------
    permutation_group:
        the i-th permutation permutation_group[i] permutes j to permutation_group[i][j]
        for j = 0,...,len(permutation_group[0])-1
    num_color: int
    num_elements_of_each_color: List of int, (num_color)

    Returns
    -------
    coeffs: int
    """
    num_elements = len(permutation_group[0])
    coeffs = 0

    for perm in permutation_group:
        type_of_perm = tuple(get_type_of_permutation(perm))
        coeffs_perm = get_inventory_coefficient(
            type_of_perm, num_elements, tuple(num_elements_of_each_color)
        )
        coeffs += coeffs_perm

    assert coeffs % len(permutation_group) == 0
    coeffs //= len(permutation_group)
    return coeffs


def get_type_of_permutation(permutation: List[int]):
    num_elements = len(permutation)
    type_of_perm = [0 for _ in range(num_elements)]

    flags = [False for _ in range(num_elements)]
    for i in range(num_elements):
        if flags[i]:
            continue
        flags[i] = True
        pos = permutation[i]
        cnt = 1
        while pos != i:
            flags[pos] = True
            pos = permutation[pos]
            cnt += 1

        type_of_perm[cnt - 1] += 1

    return type_of_perm
