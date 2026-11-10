"""Check the two exact calculations from the PMG metric parameters series.

https://github.com/magin67/pmg-notes/blob/de21abc72638400e642cb5559af19d833a16f2e6/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20parameters/Degree%20Potential%20and%20Centralities.md
"""

from fractions import Fraction as Q

import pytest
from sympy import Matrix, ones, zeros

from polyform_algebra.core import (
    DomainError, UnsupportedOperation, degree_potential, laplacian, laplacian_circumsphere,
    network, to_numeric,
)


@pytest.mark.parametrize(('edges', 'potential', 'center', 'radius'), [
    ((('a', 'b', 1), ('b', 'c', 1)),
     (Q(-1, 9), Q(2, 9), Q(-1, 9)), (Q(5, 18), Q(4, 9), Q(5, 18)), Q(25, 18)),
    ((('a', 'b', 7), ('a', 'c', 1), ('b', 'c', 1)),
     (Q(2, 3), Q(2, 3), Q(-4, 3)), (Q(2, 3), Q(2, 3), Q(-1, 3)), Q(8)),
    ((('a', 'b', 1), ('a', 'c', 1), ('b', 'c', 1)),
     (0, 0, 0), (Q(1, 3),) * 3, Q(2)),
])
def test_laplacian_circumsphere(edges, potential, center, radius):
    """The note values also solve the Poisson and equal-radius equations independently."""
    net = network(('a', 'b', 'c'), edges)
    result = laplacian_circumsphere(net).value
    assert result.degree_potential == degree_potential(net).value
    assert result.center_coordinates.context == net.context
    assert result.degree_potential.values == potential
    assert result.center_coordinates.values == center
    assert result.squared_radius == radius
    assert isinstance(result.squared_radius, Q)

    L = Matrix(laplacian(net).value.values)
    diagonal = Matrix(L.diagonal()).T
    unit = ones(3, 1)
    assert L * Matrix(potential) == diagonal - unit * sum(diagonal) / 3
    assert sum(potential) == 0
    squared_distances = diagonal * unit.T + unit * diagonal.T - 2 * L
    system = squared_distances.row_join(-2 * unit).col_join(unit.T.row_join(zeros(1)))
    assert tuple(system.LUsolve(Matrix([0, 0, 0, 1]))) == (*center, radius)

    scale = Q(7, 3)
    scaled = network(net.context.labels, tuple((a, b, scale * c) for a, b, c in edges))
    scaled_result = laplacian_circumsphere(scaled).value
    assert scaled_result.degree_potential.values == potential
    assert scaled_result.center_coordinates.values == center
    assert scaled_result.squared_radius == scale * radius


def test_degree_geometry_domain():
    """A single point is well defined; a global result needs connected exact data."""
    point = network(('a',), ())
    result = laplacian_circumsphere(point).value
    assert result.squared_radius == 0
    assert result.center_coordinates.values == (1,)
    assert degree_potential(point).value.values == (0,)
    for operation in (degree_potential, laplacian_circumsphere):
        for net in (network((), ()), network(('a', 'b'), ())):
            with pytest.raises(DomainError, match='connected'):
                operation(net)
        with pytest.raises(UnsupportedOperation, match='exact'):
            operation(to_numeric(point).value)
