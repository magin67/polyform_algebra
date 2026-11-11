"""PMG Spectra of Isomers, Spectral Transfer to a New Boundary.

Source: https://github.com/magin67/pmg-notes/tree/a55218dafc004ce1d8bd30dca60ea417d6eb5a74/pmg-library/en/Graph%20Transformations%20and%20Gluing/Spectra%20of%20Isomers
"""

import pytest
from sympy import Matrix, Symbol, expand

from polyform_algebra.boundary_data import spectral_pair_passport, spectral_pair_transfer
from polyform_algebra.core import network


@pytest.mark.parametrize("exit", [("b", "a"), ("a", "d"), ("d", "e")])
@pytest.mark.parametrize("shared", [0, 2])
def test_spectral_transfer_matches_direct_minors(exit, shared):
    x = Symbol("x")
    entry = ("a", "b")
    seam = [("a", "b", shared)] if shared else []
    left_edges = [("a", "c", 1), ("c", "b", 1)]
    right_edges = [("a", "d", 1), ("d", "e", 1), ("e", "b", 1)]
    old = network(tuple("cab"), left_edges + seam)
    kernel = network(tuple("daeb"), right_edges + seam)
    glued = network(tuple("dabce"), left_edges + right_edges + seam)
    state = Matrix(spectral_pair_passport(old, entry, x))
    expected = Matrix(spectral_pair_passport(glued, exit, x))
    transfer = spectral_pair_transfer(kernel, entry, exit, x, shared_conductance=shared)
    assert (transfer * state - expected).applyfunc(expand) == Matrix.zeros(5, 1)
    assert state[4].subs(x, 2) == 0
    assert transfer.subs(x, 2) * state.subs(x, 2) == expected.subs(x, 2)
