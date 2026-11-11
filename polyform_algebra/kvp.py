"""Exact KVP and KVPOU projections for zero-diagonal pair tables."""

from dataclasses import dataclass
from itertools import combinations

from sympy import ImmutableMatrix, Matrix, diag, ones

from .core.data import _number
from .toponomes import _forest_types


def _square(values):
    matrix = Matrix(values)
    if matrix.rows != matrix.cols or matrix.rows < 2:
        raise ValueError("expected a square matrix of order at least two")
    return matrix.applyfunc(lambda value: _number(value, "exact"))


def _pair_table(values):
    matrix = _square(values)
    if any(matrix.diagonal()):
        raise ValueError("pair table must have zero diagonal")
    return matrix


@dataclass(frozen=True)
class PairDecomposition:
    """Five components in the same vertex order as the input table."""

    K: ImmutableMatrix
    V: ImmutableMatrix
    P: ImmutableMatrix
    O: ImmutableMatrix
    U: ImmutableMatrix


def decompose(values):
    """Split an exact ordered pair table into K, V, P, O, and U."""
    matrix = _pair_table(values)
    n = matrix.rows
    unit = ones(n, 1)
    symmetric, skew = (matrix + matrix.T) / 2, (matrix - matrix.T) / 2
    mean = sum(matrix) / (n * (n - 1))
    deviations = symmetric * unit - (n - 1) * mean * unit
    imbalance = skew * unit
    K = Matrix(n, n, lambda i, j: mean if i != j else 0)
    V = Matrix(n, n, lambda i, j:
               (deviations[i] + deviations[j]) / (n - 2) if i != j and n > 2 else 0)
    U = (imbalance * unit.T - unit * imbalance.T) / n
    return PairDecomposition(*(ImmutableMatrix(component)
                               for component in (K, V, symmetric - K - V, skew - U, U)))


def pair_laplacian(values):
    """Return the outgoing Laplacian of an exact pair table, including signed weights."""
    matrix = _pair_table(values)
    return ImmutableMatrix(diag(*map(sum, matrix.tolist())) - matrix)


def distance_table(values):
    """Return the table A_ii + A_jj - 2 A_ij of an exact square matrix."""
    matrix = _square(values)
    n = matrix.rows
    return ImmutableMatrix(n, n, lambda i, j: matrix[i, i] + matrix[j, j] - 2 * matrix[i, j])


def signed_forest_toponome(values):
    """Sum forest weight products by positive type for a signed symmetric pair table."""
    matrix = _pair_table(values)
    if matrix != matrix.T:
        raise ValueError("signed undirected forests require a symmetric table")
    edges = [(i, j, matrix[i, j]) for i, j in combinations(range(matrix.rows), 2) if matrix[i, j]]
    return _forest_types(matrix.rows, edges)
