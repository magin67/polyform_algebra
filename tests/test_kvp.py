"""Exact examples from KVP Decomposition and Oriented Extensions, PMG v5."""

import pytest
from sympy import Matrix, Rational, ones, zeros

from polyform_algebra.kvp import decompose, distance_table, pair_laplacian, signed_forest_toponome


def test_path_projections():
    """The P4 degree deviations give K, V, and a zero-row-sum residual."""
    table = Matrix([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]])
    parts = decompose(table)
    assert parts.K[0, 1] == parts.V[1, 2] == parts.P[0, 1] == Rational(1, 2)
    assert parts.P * ones(4, 1) == zeros(4, 1)
    assert parts.K + parts.V + parts.P == table
    assert parts.O == parts.U == zeros(4)
    assert distance_table(pair_laplacian(table)) == 8 * parts.K + 4 * parts.V + 2 * parts.P


def test_one_arc_and_two_vertex_boundary():
    """An isolated third vertex admits circulation; on two vertices only K and U survive."""
    table = Matrix([[0, 1, 0], [0, 0, 0], [0, 0, 0]])
    parts = decompose(table)
    assert parts.O[0, 1] == Rational(1, 6)
    assert parts.U[0, 1] == Rational(1, 3)
    assert parts.O * ones(3, 1) == zeros(3, 1)
    assert parts.K + parts.V + parts.P + parts.O + parts.U == table
    two = decompose(table[:2, :2])
    assert two.V == two.P == two.O == zeros(2)


def test_signed_path_forest_cancellation():
    """The grade-three type (2,1) cancels for the P component of P5."""
    table = zeros(5)
    for i in range(4):
        table[i, i + 1] = table[i + 1, i] = 1
    assert signed_forest_toponome(decompose(table).P) == {
        (): 1, (2,): -2, (1, 1): 1, (3,): Rational(-2, 27), (4,): Rational(1, 9),
    }


def test_pair_table_validation():
    """The projectors require zero-diagonal exact pair data."""
    with pytest.raises(ValueError):
        decompose([[1, 0], [0, 1]])
    with pytest.raises(ValueError):
        decompose([[0, 0.5], [0, 0]])
    with pytest.raises(ValueError):
        signed_forest_toponome([[0, 1], [0, 0]])
