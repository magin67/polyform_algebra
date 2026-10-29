"""Geometry checks for PMG notes 5 and 9; source links are in the example."""
from itertools import permutations
from math import sqrt

import numpy as np
import pytest

from examples.introductory.note05_triangle_geometry import calculate
from polyform_algebra import Point, quForm


def test_triangle_geometry_values_and_normalization():
    """K3 has side squared 2/3, Gram determinant 1/3, and triangle area squared 1/12."""
    result = calculate()
    expected = {
        "tau": 3, "R12": 2/3, "R13": 2/3, "R23": 2/3,
        "uv": 1/3, "uw": -1/3, "vw": 1/3, "reversed_uv": -1/3,
        "gram_determinant": 1/3, "boundary_norm": 1/3,
        "boundary_linear_norm": 1/sqrt(3), "area_squared": 1/12,
        "area": sqrt(3)/6, "product_form_norm": 1/3, "product_of_norms": 4/9,
    }
    for name, value in expected.items():
        assert result[name] == pytest.approx(value), name
    np.testing.assert_allclose(result["gram"], [[2/3, 1/3], [1/3, 2/3]])
    # Evaluation is not multiplicative: the missing cross term is (u dot v)^2.
    assert result["product_of_norms"]-result["product_form_norm"] == pytest.approx(1/9)


def test_triangle_geometry_against_independent_matrix_and_coordinates():
    """A grounded Laplacian and an explicit equilateral embedding give independent references."""
    result = calculate()
    # Ground vertex 1 in the explicit K3 Laplacian, without using polyform conversion.
    reduced = np.array([[2, -1], [-1, 2]], dtype=float)
    gram = np.linalg.solve(reduced, np.eye(2))
    np.testing.assert_allclose(result["gram"], gram)
    assert result["boundary_norm"] == pytest.approx(np.linalg.det(gram))
    # Resistance distances are sqrt(2/3), not the unit edge resistance.
    side = sqrt(2/3)
    coordinates = np.array([[0, 0], [side, 0], [side/2, side*sqrt(3)/2]])
    edges = (coordinates[1:] - coordinates[0]).T
    np.testing.assert_allclose(result["gram"], edges.T @ edges)
    area = abs(np.linalg.det(edges))/2
    assert result["area"] == pytest.approx(area)
    assert result["area_squared"] == pytest.approx(3*side**4/16)


@pytest.mark.parametrize("order", list(permutations(range(3))))
def test_boundary_norm_is_independent_of_vertex_order(order):
    """Reversing orientation or choosing another base vertex preserves the squared area."""
    points = Point.create_list(["1", "2", "3"])
    a, b, c = points
    metric = (quForm([(a, b)]) + quForm([(a, c)]) + quForm([(b, c)])).exp()
    x, y, z = (points[i] for i in order)
    boundary_norm = metric.norm(quForm([(x, y, z)]))
    u, v = y-x, z-x
    gram_determinant = metric.norm(u @ u)*metric.norm(v @ v)-metric.norm(u @ v)**2
    assert boundary_norm == pytest.approx(1/3)
    assert boundary_norm == pytest.approx(gram_determinant)
