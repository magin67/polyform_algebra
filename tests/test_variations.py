"""Exact finite examples from PMG Variations in pmg-notes v5 at 63d4a45.

https://github.com/magin67/pmg-notes/tree/63d4a45da9f3a3cbe46f57f4e91c838f9a828c0c/pmg-library/en/Graph%20Transformations%20and%20Gluing
"""

from fractions import Fraction as Q

import pytest
from sympy import Matrix

from polyform_algebra.core import DomainError, InvalidInput, green, network, spanning_tree_weight
from polyform_algebra.variations import (
    conductance_greenian_update, resistance_constraint_polynomials,
    vary_conductances, vary_greenian_coefficients, vary_resistances,
)


def test_edge_update_and_disconnect():
    path = network(tuple("abc"), [("a", "b", 1), ("b", "c", 1)])
    changes = {("b", "a"): -1, ("c", "a"): 2}
    updated = vary_conductances(path, changes)
    result, ratio = conductance_greenian_update(path, changes)
    assert result == Matrix(green(updated).value.values)
    assert ratio == spanning_tree_weight(updated).value == 2
    assert conductance_greenian_update(path, {}) == (Matrix(green(path).value.values), 1)
    with pytest.raises(DomainError):
        conductance_greenian_update(path, {("a", "b"): -1})
    with pytest.raises(InvalidInput):
        vary_conductances(path, {("a", "b"): 1, ("b", "a"): 1})


def test_inverse_variations():
    path = network(tuple("abcd"), [("a", "b", 1), ("b", "c", 1), ("c", "d", 1)])
    updated = vary_greenian_coefficients(path, {("a", "b"): Q(1, 5)})
    assert updated.edges == (("a", "b", Q(2, 5)), ("a", "c", Q(1, 5)),
                             ("b", "c", Q(7, 10)), ("c", "d", 1))
    triangle = network(tuple("abc"), [("a", "b", 1), ("a", "c", 1), ("b", "c", 1)])
    updated = vary_resistances(triangle, {("a", "b"): Q(1, 3), ("a", "c"): Q(1, 3)})
    assert updated.edges == (("a", "b", Q(3, 5)), ("a", "c", Q(3, 5)), ("b", "c", Q(6, 5)))
    with pytest.raises(InvalidInput):
        vary_resistances(triangle, {("a", "b"): 2})


def test_exact_constraints_clear_tree_denominator():
    path = network(tuple("abc"), [("a", "b", 1), ("b", "c", 1)])
    (c,), tree, (residual,) = resistance_constraint_polynomials(path, [("a", "b")], {("a", "c"): Q(3, 2)})
    assert tree == c
    assert residual == 1 - c/2
    assert residual.subs(c, 2) == 0 and tree.subs(c, 2) == 2
