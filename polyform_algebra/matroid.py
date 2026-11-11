"""Exact forest invariants of undirected conductance networks."""

from fractions import Fraction

from .core import DomainError, resistance
from .core.data import _number
from .toponomes import _require_exact_network, forest_toponome


def forest_coefficients(net):
    """Return ordinary spanning-forest coefficients in ascending edge count.

    The highest grade is the rank of the cycle matroid. Like forest_toponome,
    this calculation enumerates forests and is intended for small networks.
    """
    by_type = forest_toponome(net)
    coefficients = [Fraction(0)] * (1 + max(map(sum, by_type)))
    for parts, weight in by_type.items():
        coefficients[sum(parts)] += weight
    return tuple(coefficients)


def normalized_forest_coefficients(net):
    """Reverse the ordinary forest coefficients and divide by the top coefficient."""
    coefficients = forest_coefficients(net)
    return tuple(value / coefficients[-1] for value in reversed(coefficients))


def forest_resistance_scales(coefficients):
    """Return adjacent ratios, starting at the highest degree of a forest polynomial."""
    values = tuple(_number(value, "exact") for value in coefficients)
    if not values or values[0] != 1 or any(value <= 0 for value in values):
        raise ValueError("expected positive forest coefficients starting with one")
    return tuple(values[i - 1] / values[i] for i in range(len(values) - 1, 0, -1))


def forest_power_sums(normalized_coefficients, count):
    """Compute the first count Newton sums from the normalized forest polynomial."""
    coefficients = tuple(_number(value, "exact") for value in normalized_coefficients)
    if not coefficients or coefficients[0] != 1:
        raise ValueError("expected normalized coefficients starting with one")
    if type(count) is not int or count < 0:
        raise ValueError("count must be a nonnegative integer")
    degree = len(coefficients) - 1
    sums = []
    for order in range(1, count + 1):
        value = sum(
            (-1) ** (j - 1) * coefficients[j] * sums[order - j - 1]
            for j in range(1, min(order - 1, degree) + 1)
        )
        if order <= degree:
            value += (-1) ** (order - 1) * order * coefficients[order]
        sums.append(Fraction(value))
    return tuple(sums)


def _edge_resistances(net):
    _require_exact_network(net)
    pairs = tuple((a, b) for a, b, _ in net.edges)
    return resistance(net, pairs=pairs).value.values


def edge_forest_probabilities(net):
    """Map each network edge to its inclusion probability in a weighted spanning forest.

    A spanning tree is sampled in each connected component. Parallel
    conductances already combined by network() are treated as one edge.
    """
    return {(a, b): conductance * value for (a, b, conductance), value
            in zip(net.edges, _edge_resistances(net), strict=True)}


def forest_resistance_radius(net):
    """Return the squared resistance circumradius via the two-forest formula."""
    _require_exact_network(net)
    if len(net.components) != 1:
        raise DomainError("forest resistance radius requires a nonempty connected network")
    if not net.edges:
        return Fraction(0)
    coefficients = forest_coefficients(net)
    xi = coefficients[-2] / coefficients[-1]
    h2 = sum(conductance * value**2 for (_, _, conductance), value
             in zip(net.edges, _edge_resistances(net), strict=True))
    return xi / 3 - h2 / 12
