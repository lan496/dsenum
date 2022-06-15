import pytest

from dsenum.permutation_group import (
    DerivativeStructurePermutation,
    is_permutation_group,
)
from dsenum.superlattice import generate_symmetry_distinct_superlattices
from dsenum.utils import get_lattice


@pytest.mark.parametrize("kind", ["fcc", "bcc", "sc", "hex", "tet", "hcp"])
def test_permutations(kind):
    structure = get_lattice(kind)
    frac_coords = structure.frac_coords
    for index in range(1, 10 + 1):
        list_reduced_HNF, rotations, translations = generate_symmetry_distinct_superlattices(
            index, structure, return_symops=True
        )
        for hnf in list_reduced_HNF:
            dsperm = DerivativeStructurePermutation(hnf, frac_coords, rotations, translations)
            assert is_permutation_group(dsperm.prm_t)
            prm_all = dsperm.get_symmetry_operation_permutations()
            assert is_permutation_group(prm_all)
