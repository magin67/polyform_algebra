"""Boundary Polynomials and Chains of Cores, PMG notes v5 at 4899950.

https://github.com/magin67/pmg-notes/tree/4899950bb3620eb2294ad3d1a5872152d071015a/pmg-library/en/Graph%20Transformations%20and%20Gluing
"""

import pytest
from sympy import Symbol, expand

from polyform_algebra.boundary_data import (
    forest_passport, join_forest_states, pair_forest_passports,
    pair_spectral_passports, spectral_pair_passport,
)
from polyform_algebra.core import InvalidInput, network
from polyform_algebra.core_chains import cycle_chain_tree_count


def test_pair_passports_and_boundary_order():
    first = network(tuple("axb"), [("a", "x", 1), ("x", "b", 1)])
    second = network(tuple("aby"), [("a", "y", 1), ("y", "b", 1)])
    combined = network(tuple("axby"), first.edges + second.edges)
    boundary = tuple("ba")
    left, right = (forest_passport(net, boundary) for net in (first, second))
    assert left == {(("b",), ("a",)): {0: 1, 1: 2}, (("b", "a"),): {2: 1}}
    assert pair_forest_passports(boundary, left, right) == forest_passport(combined, boundary)
    assert join_forest_states(boundary, (("a", "b"),), (("b", "a"),)) is None
    with pytest.raises(InvalidInput):
        forest_passport(first, ("x", "x"))


@pytest.mark.parametrize("shared", [0, 2])
def test_spectral_pair_gluing(shared):
    x = Symbol("x")
    seam = [("a", "b", shared)] if shared else []
    first = network(tuple("acb"), [("a", "c", 1), ("c", "b", 1)] + seam)
    second = network(tuple("dba"), [("a", "d", 3), ("d", "b", 2)] + seam)
    combined = network(tuple("acbd"), [("a", "c", 1), ("c", "b", 1), ("a", "d", 3), ("d", "b", 2)] + seam)
    boundary = tuple("ba")
    left, right = (spectral_pair_passport(net, boundary, x) for net in (first, second))
    assert pair_spectral_passports(left, right, x, shared_conductance=shared) == spectral_pair_passport(combined, boundary, x)
    determinant, minor_a, minor_b, cross, interior = left
    assert expand(determinant * interior - minor_a * minor_b + cross**2) == 0


def test_cycle_chain_continuant():
    assert cycle_chain_tree_count([4, 4, 4], [1, 1]) == 56
    assert cycle_chain_tree_count([6, 6, 6], [1, 1]) == 204
    assert cycle_chain_tree_count([3, 4], [0]) == 12
    with pytest.raises(ValueError):
        cycle_chain_tree_count([3, 3, 3], [1, 1])
