import os
from time import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pymatgen.core import Lattice, Structure

from dsenum import ZddStructureEnumerator


# General ploter utility
def get_custom_rcparams():
    # based on seaborn-paper context
    customrc = {
        "axes.labelsize": 8.8,
        "axes.titlesize": 9.6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        # default grid.linewidth=0.8 is too bold
        "grid.linewidth": 0.2,
        "lines.linewidth": 1.4,
        "patch.linewidth": 0.24,
        "lines.markersize": 5.6,
        "lines.markeredgewidth": 0,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "xtick.major.pad": 5.6,
        "ytick.major.pad": 5.6,
        # font
        "font.family": "serif",
        "font.serif": "Times New Roman",
        "mathtext.fontset": "cm",  # Computer Modern
        # spines
        "axes.spines.top": False,
        "axes.spines.right": False,
        # grid
        "axes.axisbelow": True,
        # ticks
        "xtick.direction": "in",
        "ytick.direction": "in",
    }
    return customrc


def get_fcc_structure():
    a = 2.856

    lattice = Lattice.from_parameters(a, a, a, 60, 60, 60)
    species = ["Cu"]
    frac_coords = np.array([[0, 0, 0]])
    structure = Structure(lattice, species, frac_coords)
    return structure


def measure_listing_time(aristo, index, num_types):
    zse = ZddStructureEnumerator(
        aristo,
        index,
        num_types,
        remove_superperiodic=True,
        remove_incomplete=True,
    )

    start = time()
    list_dstructs = zse.generate(output="poscar")
    end = time()
    return end - start, len(list_dstructs)


def measure_counting_time(aristo, index, num_types):
    zse = ZddStructureEnumerator(
        aristo,
        index,
        num_types,
        remove_superperiodic=True,
        remove_incomplete=True,
    )

    start = time()
    count = zse.count()
    end = time()
    return end - start, count


def plotter(fig, ax, df, list_num_types):
    for i, num_types in enumerate(list_num_types):
        df_tmp = df[(df["count"] > 0) & (df["num_types"] == num_types)].sort_values(by=["count"])

        list_count = df_tmp["count"]
        list_listing_sec = df_tmp["listing_sec"]
        list_counting_sec = df_tmp["counting_sec"]

        ax[i].scatter(list_count, list_listing_sec, label="Listing", color="C0", marker="o")
        ax[i].plot(list_count, list_listing_sec, color="C0", linestyle="--")

        ax[i].scatter(list_count, list_counting_sec, label="Counting", color="C1", marker="s")
        ax[i].plot(list_count, list_counting_sec, color="C1", linestyle="--")

        ax[i].set_xscale("log")
        ax[i].set_xlabel("Number of derivative structures")
        ax[i].set_title(f"Number of Species = {num_types}")

    ax[0].set_ylabel("Run time (sec)")
    ax[0].set_yscale("log")
    ax[0].legend()

    yticks = [10**i for i in range(-3, 2 + 1)]
    ax[0].set_yticks(yticks)
    ax[0].set_ylim(min(yticks), max(yticks))

    xticks = [10**i for i in range(0, 6 + 1)]
    ax[0].set_xticks(xticks)
    ax[0].set_xlim(min(xticks), max(xticks))
    return fig, ax


if __name__ == "__main__":
    aristo = get_fcc_structure()

    all_data = []

    num_types_and_max_index = [(2, 15), (3, 10), (4, 10)]
    for num_types, max_index in num_types_and_max_index:
        for index in range(1, max_index + 1):
            counting_elapsed, count = measure_counting_time(aristo, index, num_types)
            listing_elapsed, count2 = measure_listing_time(aristo, index, num_types=num_types)
            assert count == count2
            data = {
                "index": index,
                "num_types": num_types,
                "count": count,
                "counting_sec": counting_elapsed,
                "listing_sec": listing_elapsed,
            }
            all_data.append(data)

    df = pd.DataFrame(all_data)

    with plt.rc_context(get_custom_rcparams()):
        list_num_types = [2, 3, 4]
        fig, ax = plt.subplots(
            1,
            len(list_num_types),
            figsize=(3.0 * len(list_num_types), 3.375),
            dpi=100,
            sharey=True,
            sharex=True,
        )
        fig, ax = plotter(fig, ax, df, list_num_types)
        plt.tight_layout()
        fname = "fcc_enumeration_runtime.pdf"
        fig.savefig(fname, bbox_inches="tight", pad_inches=0.0, dpi=500)
