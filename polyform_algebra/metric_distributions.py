"""Exact metric distributions on exterior coordinates of a connected network."""

from dataclasses import dataclass
from itertools import combinations
from math import prod

from sympy import ImmutableMatrix, Matrix, eye, zeros

from .core import DomainError, InvalidInput, laplacian
from .core.data import _number
from .toponomes import _require_exact_network


def _dimension(net):
    _require_exact_network(net)
    if len(net.components) != 1:
        raise DomainError("metric distributions require a nonempty connected network")
    return len(net.context.labels) - 1


def _direction(labels, pair):
    first, second = pair
    if first == second or first not in labels or second not in labels:
        raise InvalidInput("expected a pair of distinct vertex labels")
    # Coordinates in the basis e_i - e_ref, with the last vertex as reference.
    return Matrix(len(labels) - 1, 1, [int(v == second) - int(v == first) for v in labels[:-1]])


def _blade(columns, indices):
    return Matrix(len(indices), 1, [columns.extract(rows, range(columns.cols)).det() for rows in indices])


def pair_blade(net, pairs):
    """Return coordinates of the wedge of the supplied ordered vertex differences."""
    dimension = _dimension(net)
    directions = tuple(_direction(net.context.labels, pair) for pair in pairs)
    if len(directions) > dimension:
        raise InvalidInput("exterior grade exceeds the zero-sum dimension")
    columns = Matrix.hstack(*directions) if directions else zeros(dimension, 0)
    return ImmutableMatrix(_blade(columns, tuple(combinations(range(dimension), len(directions)))))


@dataclass(frozen=True)
class MetricDistribution:
    """Exterior metric and its coefficient matrices indexed by edge subsets."""

    metric: ImmutableMatrix
    blades: dict
    functionals: dict

    def _inputs(self, left, right):
        vectors = []
        for values in (left, left if right is None else right):
            column = Matrix(values)
            if column.shape != (self.metric.rows, 1):
                raise InvalidInput("expected a column of exterior coordinates")
            vectors.append(column.applyfunc(lambda value: _number(value, "exact")))
        return vectors

    def evaluate(self, left, right=None):
        """Evaluate the exterior metric on two coordinate columns."""
        left, right = self._inputs(left, right)
        return (left.T * self.metric * right)[0]

    def coefficients(self, left, right=None):
        """Evaluate all coefficients of a quadratic or mixed exterior input."""
        left, right = self._inputs(left, right)
        return {edges: (left.T * row * right)[0] for edges, row in self.functionals.items()}


def metric_distribution(net, grade=1):
    """Return the exact exterior metric and all grade-sized edge contributions.

    Coordinates use the zero-sum basis relative to the last network vertex.
    The coefficients sum to the metric evaluation. Edge-subset enumeration
    is intended for small graphs.
    """
    dimension = _dimension(net)
    if type(grade) is not int or not 0 <= grade <= dimension:
        raise InvalidInput("grade must lie between zero and the zero-sum dimension")
    indices = tuple(combinations(range(dimension), grade))
    if grade:
        inverse = Matrix(laplacian(net).value.values)[:-1, :-1].inv()
        metric = ImmutableMatrix([[inverse.extract(rows, cols).det() for cols in indices] for rows in indices])
    else:
        metric = ImmutableMatrix(eye(1))

    labels = net.context.labels
    directions = {(a, b): _direction(labels, (a, b)) for a, b, _ in net.edges}
    weights = {(a, b): weight for a, b, weight in net.edges}
    blades, functionals = {}, {}
    for edges in combinations(directions, grade):
        columns = Matrix.hstack(*(directions[edge] for edge in edges)) if edges else zeros(dimension, 0)
        blade = ImmutableMatrix(_blade(columns, indices))
        image = metric * blade
        blades[edges] = blade
        functionals[edges] = ImmutableMatrix(prod(weights[edge] for edge in edges) * image * image.T)
    return MetricDistribution(metric, blades, functionals)
