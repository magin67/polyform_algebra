"""Exact partition toponomes and the second-kind product of types."""

from collections import Counter
from itertools import combinations, permutations
from math import prod

from sympy import Matrix

from .combinatorics import _lengths
from .core.data import InvalidInput, Network, _labels, _number
from .core.network_math import laplacian


def _positive_type(parts):
    return tuple(sorted(_lengths(parts), reverse=True))


def _collect(terms):
    coefficients = Counter()
    for kind, coefficient in terms:
        coefficients[kind] += coefficient
    return {kind: coefficient for kind, coefficient in coefficients.items() if coefficient}


def partition_toponome(terms):
    """Reduce chosen terms given as (block sizes, rational coefficient) pairs.

    Each listed block size is positive. Singleton blocks are omitted in the
    result, and a block of size b contributes a part of grade b - 1.
    """
    return _collect(
        (tuple(size - 1 for size in _positive_type(sizes) if size > 1),
         _number(coefficient, "exact"))
        for sizes, coefficient in terms
    )


def _require_exact_network(net):
    if not isinstance(net, Network) or net.context.arithmetic != "exact":
        raise TypeError("expected a network with exact conductances")


def forest_toponome(net):
    """Sum exact conductance products over spanning forests by positive type.

    The empty forest has coefficient one. The enumeration is exponential in
    general and is intended for small undirected graphs.
    """
    _require_exact_network(net)
    labels = net.context.labels
    positions = {label: position for position, label in enumerate(labels)}
    edges = [(positions[a], positions[b], weight) for a, b, weight in net.edges]

    def forests(start, components, weight):
        sizes = Counter(components).values()
        yield tuple(sorted((size - 1 for size in sizes if size > 1), reverse=True)), weight
        for offset in range(start, len(edges)):
            a, b, conductance = edges[offset]
            first, second = components[a], components[b]
            if first != second:
                merged = tuple(first if component == second else component for component in components)
                yield from forests(offset + 1, merged, weight * conductance)

    return _collect(forests(0, tuple(range(len(labels))), 1))


def _matched_types(left, right):
    for size in range(min(len(left), len(right)) + 1):
        for first in combinations(range(len(left)), size):
            for second in permutations(range(len(right)), size):
                parts = [left[i] + right[j] for i, j in zip(first, second)]
                parts.extend(value for i, value in enumerate(left) if i not in first)
                parts.extend(value for j, value in enumerate(right) if j not in second)
                yield tuple(sorted(parts, reverse=True))


def partition_product(left, right):
    """Multiply two rational type-to-coefficient mappings by partial matching.

    Parts are positive integers. Occurrences of equal parts are distinguished
    during matching. The unit is {(): 1}.
    """
    first = _collect((_positive_type(kind), _number(value, "exact")) for kind, value in left.items())
    second = _collect((_positive_type(kind), _number(value, "exact")) for kind, value in right.items())
    return _collect(
        (kind, a * b)
        for left_kind, a in first.items()
        for right_kind, b in second.items()
        for kind in _matched_types(left_kind, right_kind)
    )


def color_class_toponome(net, vertex_weights=None):
    """Count unordered proper color-class partitions by total class weights.

    A weight mapping, if supplied, must give every vertex a positive integer
    weight. Unit vertex weights are used by default. Edge conductances are
    irrelevant to this count.
    """
    if not isinstance(net, Network):
        raise TypeError("expected a network")
    labels = net.context.labels
    if vertex_weights is None:
        vertex_weights = dict.fromkeys(labels, 1)
    if set(vertex_weights) != set(labels):
        raise ValueError("vertex weights must cover exactly the network vertices")
    weights = dict(zip(labels, _lengths(vertex_weights[label] for label in labels)))
    neighbours = {label: set() for label in labels}
    for a, b, _ in net.edges:
        neighbours[a].add(b)
        neighbours[b].add(a)

    def colorings(position, classes):
        if position == len(labels):
            totals = (sum(weights[label] for label in block) for block in classes)
            yield tuple(sorted(totals, reverse=True)), 1
            return
        label = labels[position]
        singleton = {label}
        yield from colorings(position + 1, classes + (singleton,))
        for i, block in enumerate(classes):
            if block.isdisjoint(neighbours[label]):
                extended = classes[:i] + (block | singleton,) + classes[i + 1:]
                yield from colorings(position + 1, extended)

    return _collect(colorings(0, ()))


def coefficient_table(toponome):
    """Collect a rational toponome by (total grade, number of positive parts)."""
    terms = ((_lengths(kind), _number(coefficient, "exact"))
             for kind, coefficient in toponome.items())
    return _collect(
        ((sum(parts), len(parts)), coefficient) for parts, coefficient in terms
    )


def vertex_toponome(toponome):
    """Apply the root-choice weight to each type before any grouping.

    A part of grade k contributes k + 1 possible roots. The empty type has
    weight one. Coefficients may be signed rational numbers.
    """
    terms = ((_positive_type(kind), _number(coefficient, "exact"))
             for kind, coefficient in toponome.items())
    return _collect((parts, coefficient * prod(k + 1 for k in parts))
                    for parts, coefficient in terms)


def laplacian_polynomial(net):
    """Return the exact coefficients of det(I + t L), constant term first.

    The matrix route uses the characteristic polynomial of L, without
    enumerating forests. Trailing zero coefficients are omitted.
    """
    _require_exact_network(net)
    coefficients = Matrix(laplacian(net).value.values).charpoly().all_coeffs()
    while len(coefficients) > 1 and coefficients[-1] == 0:
        coefficients.pop()
    return tuple(_number((-1)**k * value, "exact") for k, value in enumerate(coefficients))


def rooted_forest_weight(net, roots):
    """Sum weights of spanning forests with exactly one prescribed root per component.

    Roots are unique vertex labels. The principal minor uses the rows and
    columns of nonroots in the original Laplacian. If every vertex is a root,
    the empty minor has determinant one.
    """
    _require_exact_network(net)
    roots = set(_labels(roots))
    labels = net.context.labels
    if not roots.issubset(labels):
        raise InvalidInput("roots must be network vertex labels")
    keep = [i for i, label in enumerate(labels) if label not in roots]
    matrix = Matrix(laplacian(net).value.values)
    return _number(matrix.extract(keep, keep).det(), "exact")
