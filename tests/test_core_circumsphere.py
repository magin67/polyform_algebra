"""Exact center and radius examples from the PMG metric parameters series.

https://github.com/magin67/pmg-notes/blob/de21abc72638400e642cb5559af19d833a16f2e6/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20parameters/Sphere%20Centers%20and%20Graph%20Heterogeneity.md
"""

from fractions import Fraction as Q
from itertools import combinations

import pytest
from sympy import Matrix, ones, zeros

from polyform_algebra.core import (
    DomainError, UnsupportedOperation, circumsphere, kron_reduce, network, resistance, to_numeric,
)


@pytest.mark.parametrize(('vertices', 'edges', 'center', 'radius'), [
    (('a', 'b', 'c', 'd'), (('a', 'b', 1), ('b', 'c', 1), ('c', 'd', 1)),
     (Q(1, 2), 0, 0, Q(1, 2)), Q(3, 4)),
    (('a', 'b', 'c', 'd'), (('a', 'b', 1), ('a', 'c', 1), ('a', 'd', 1)),
     (Q(-1, 2), Q(1, 2), Q(1, 2), Q(1, 2)), Q(3, 4)),
    (('a', 'b', 'c'), (('a', 'b', 1), ('a', 'c', 2), ('b', 'c', 3)),
     (Q(9, 22), Q(8, 22), Q(5, 22)), Q(15, 121)),
])
def test_exact_circumsphere(vertices, edges, center, radius):
    """The exact note values also solve the independent equal-radius system."""
    net = network(vertices, edges)
    result = circumsphere(net).value
    assert result.center_coordinates.context == net.context
    assert result.center_coordinates.values == center
    assert result.resistance_radius == radius
    assert isinstance(result.resistance_radius, Q)

    R = Matrix(resistance(net).value.values)
    column = ones(len(vertices), 1)
    system = R.row_join(-2 * column).col_join(column.T.row_join(zeros(1)))
    independent = system.LUsolve(Matrix([0] * len(vertices) + [1]))
    assert tuple(independent) == (*center, radius)

    scaled = network(vertices, tuple((a, b, 7 * c) for a, b, c in edges))
    scaled_result = circumsphere(scaled).value
    assert scaled_result.center_coordinates.values == center
    assert scaled_result.resistance_radius == radius / 7


def test_circumsphere_after_kron_reduction():
    """A three-vertex face of the unit K4 has side resistance 1/2 and radius 1/6."""
    labels = ('a', 'b', 'c', 'd')
    tetrahedron = network(labels, tuple((a, b, 1) for a, b in combinations(labels, 2)))
    reduced = kron_reduce(tetrahedron, ('c', 'a', 'b')).value.network
    sphere = circumsphere(reduced).value
    assert sphere.center_coordinates.context.labels == ('c', 'a', 'b')
    assert sphere.center_coordinates.values == (Q(1, 3),) * 3
    assert sphere.resistance_radius == Q(1, 6)


def test_circumsphere_domain():
    """A point has radius zero; a finite global simplex needs connected exact data."""
    point = network(('a',), ())
    sphere = circumsphere(point).value
    assert sphere.center_coordinates.values == (1,)
    assert sphere.resistance_radius == 0
    for net in (network((), ()), network(('a', 'b'), ())):
        with pytest.raises(DomainError, match='connected'):
            circumsphere(net)
    with pytest.raises(UnsupportedOperation, match='exact'):
        circumsphere(to_numeric(point).value)
