"""PMG Metric Distributions and Edge Influence, pmg-notes v5 at 7cebaa5.

https://github.com/magin67/pmg-notes/tree/7cebaa55451d801c44be39ca1149c69bd00dd6fd/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20Distributions%20and%20Edge%20Influence
"""

from fractions import Fraction as Q

import pytest
from sympy import zeros

from polyform_algebra.core import DomainError, InvalidInput, network
from polyform_algebra.metric_distributions import metric_distribution, pair_blade


def test_weighted_tree_mixed_field_and_orientation():
    path = network(tuple("abcde"), [("a", "b", 3), ("b", "c", 2), ("c", "d", Q(1, 2)), ("d", "e", 4)])
    data = metric_distribution(path)
    first, second = (pair_blade(path, (pair,)) for pair in (("a", "d"), ("b", "e")))
    weights = data.coefficients(first, second)
    assert tuple(weights.values()) == (0, Q(1, 2), 2, 0)
    assert sum(weights.values()) == data.evaluate(first, second) == Q(5, 2)
    assert data.coefficients(-first, second) == {key: -value for key, value in weights.items()}


def test_second_grade_sum_and_scalar_grade():
    graph = network(tuple("abcd"), [("a", "b", 1), ("b", "c", 1), ("c", "d", 1), ("d", "a", 1)])
    data = metric_distribution(graph, 2)
    assert sum(data.functionals.values(), zeros(3)) == data.metric
    blade = pair_blade(graph, (("a", "b"), ("a", "c")))
    assert data.evaluate(blade) == sum(data.coefficients(blade).values()) == Q(1, 2)
    assert metric_distribution(graph, 0).coefficients([2], [3]) == {(): 6}
    with pytest.raises(InvalidInput):
        metric_distribution(graph, 4)
    with pytest.raises(DomainError):
        metric_distribution(network(tuple("ab"), []))
