"""Examples 2–6 of the PMG Matroid Parameters series in pmg-notes v5.

Source: https://github.com/magin67/pmg-notes/tree/8f30a5b9440c20ed841b50c9e0ce40241dd65e02/pmg-library/en/Invariants%2C%20Decompositions%2C%20and%20Spectra/Matroid%20parameters%20and%20forest%20characteristics%20of%20graphs
"""

from fractions import Fraction as Q
from itertools import combinations

import pytest

from polyform_algebra.core import DomainError, circumsphere, network
from polyform_algebra.matroid import (
    edge_forest_probabilities, forest_coefficients, forest_power_sums,
    forest_resistance_radius, forest_resistance_scales, normalized_forest_coefficients,
)
from polyform_algebra.toponomes import laplacian_polynomial


def test_ordinary_and_rooted_forest_scales():
    complete = network(tuple("abcd"), [(a, b, 1) for a, b in combinations("abcd", 2)])
    ordinary = forest_coefficients(complete)
    assert ordinary == (1, 6, 15, 16)
    assert forest_resistance_scales(ordinary) == (Q(15, 16), Q(2, 5), Q(1, 6))
    assert forest_resistance_scales(laplacian_polynomial(complete)) == (Q(3, 4), Q(1, 4), Q(1, 12))


def test_normalization_and_newton_sums():
    tree = network(tuple("abcd"), [("a", "b", 1), ("b", "c", 2), ("c", "d", 3)])
    triangle = network(tuple("abc"), [("a", "b", 1), ("a", "c", 1), ("b", "c", 1)])
    normalized = normalized_forest_coefficients(tree)
    assert normalized == (1, Q(11, 6), 1, Q(1, 6))
    assert forest_power_sums(normalized, 5) == tuple(
        sum(Q(1, weight)**k for weight in (1, 2, 3)) for k in range(1, 6)
    )
    assert forest_power_sums(normalized_forest_coefficients(triangle), 4) == (1, Q(1, 3), 0, Q(-1, 9))
    assert forest_power_sums((1,), 2) == (0, 0)
    with pytest.raises(ValueError):
        forest_power_sums(normalized, -1)


def test_edge_probabilities_and_radius():
    triangle = network(tuple("abc"), [("a", "b", 1), ("a", "c", 2), ("b", "c", 3)])
    assert edge_forest_probabilities(triangle) == {("a", "b"): Q(5, 11), ("a", "c"): Q(8, 11), ("b", "c"): Q(9, 11)}
    assert forest_resistance_radius(triangle) == circumsphere(triangle).value.resistance_radius
    assert forest_resistance_radius(network(tuple("a"), [])) == 0
    disconnected = network(tuple("abcd"), [("a", "b", 2), ("c", "d", 3)])
    assert forest_coefficients(disconnected) == (1, 5, 6)
    assert edge_forest_probabilities(disconnected) == {("a", "b"): 1, ("c", "d"): 1}
    with pytest.raises(DomainError):
        forest_resistance_radius(disconnected)
