from typing import List
from queue import Queue

import joblib
import matplotlib as mpl

mpl.use("Agg")  # noqa
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sympy.combinatorics import Permutation
from sympy.combinatorics.perm_groups import PermutationGroup
from tqdm import tqdm

from dsenum.permutation_group import DerivativeStructurePermutation
from dsenum.superlattice import generate_all_superlattices
from dsenum.utils import get_lattice, get_symmetry_operations, get_fcc_with_vacancy


def get_generators(perms: List[Permutation]):
    identity = Permutation([], size=perms[0].size)

    gens = []
    perms_tmp = [
        identity,
    ]

    for g in perms:
        if g in perms_tmp:
            continue

        gens.append(g)

        que = Queue()
        for h in perms_tmp:
            que.put(h)

        while not que.empty():
            h = que.get()

            gh = g * h
            if gh not in perms_tmp:
                que.put(gh)
                perms_tmp.append(gh)

            hg = h * g
            if hg not in perms_tmp:
                que.put(hg)
                perms_tmp.append(hg)

    assert len(perms_tmp) == len(perms)
    return gens


def plot_permutation_groups(base_structure, output_name):
    rotations, translations = get_symmetry_operations(base_structure)
    displacement_set = base_structure.frac_coords

    list_data = []
    for index in range(2, 10 + 1):
        for hnf in tqdm(generate_all_superlattices(index)):
            dsperm = DerivativeStructurePermutation(hnf, displacement_set, rotations, translations)
            perms = [Permutation(g) for g in dsperm.get_symmetry_operation_permutations()]
            gens = get_generators(perms)
            G = PermutationGroup(gens)
            data = {
                "index": index,
                "hnf": hnf,
                "generators": gens,
                "group": G,
                "order": G.order(),
                "orbits": G.orbits(),
            }
            list_data.append(data)

    df = pd.DataFrame(list_data)
    sns.swarmplot(x="index", y="order", data=df)
    plt.title(output_name)
    plt.savefig(output_name + ".png")
    return df


if __name__ == "__main__":
    Permutation.print_cyclic = False

    print("fcc")
    base_structure = get_lattice("fcc")
    df_fcc = plot_permutation_groups(base_structure, "fcc")
    joblib.dump(df_fcc, "fcc.df.joblib")

    print("fcc with tetragonal and octhedral sites")
    base_structure = get_fcc_with_vacancy()
    df_fcc_ot = plot_permutation_groups(base_structure, "fcc_oct_tet")
    joblib.dump(df_fcc_ot, "fcc_oct_tet.df.joblib")
