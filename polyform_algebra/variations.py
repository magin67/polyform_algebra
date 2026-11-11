"""Exact finite variations of connected conductance networks."""

from sympy import Matrix, Symbol, diag, expand, eye, ones, zeros

from .core import (
    DomainError, InvalidInput, green, laplacian, network_from_laplacian,
    network_from_resistance, resistance,
)
from .core.data import _number
from .toponomes import _require_exact_network


def _connected(net):
    _require_exact_network(net)
    if len(net.components) != 1:
        raise DomainError("variations require a nonempty connected network")
    return net.context.labels


def _pair_indices(labels, pairs):
    positions = {label: i for i, label in enumerate(labels)}
    seen = set()
    for first, second in pairs:
        if first == second or first not in positions or second not in positions:
            raise InvalidInput("a variation requires two distinct known labels")
        indices = tuple(sorted((positions[first], positions[second])))
        if indices in seen:
            raise InvalidInput("duplicate unordered pair")
        seen.add(indices)
        yield indices


def _direction(size, pair):
    first, second = pair
    result = zeros(size, 1)
    result[first], result[second] = 1, -1
    return result


def _rank_update(labels, changes):
    pairs = tuple(_pair_indices(labels, changes))
    columns = Matrix.hstack(*(_direction(len(labels), pair) for pair in pairs)) if pairs else zeros(len(labels), 0)
    return columns, diag(*(_number(value, "exact") for value in changes.values()))


def vary_conductances(net, changes):
    """Apply simultaneous conductance increments and validate the new network.

    An absent edge starts at zero. The returned network may be disconnected.
    """
    labels = _connected(net)
    columns, increments = _rank_update(labels, changes)
    updated = Matrix(laplacian(net).value.values) + columns * increments * columns.T
    return network_from_laplacian(updated.tolist(), labels, policy=net.policy).value


def conductance_greenian_update(net, changes):
    """Return the exact updated Greenian and tree-weight ratio via Woodbury.

    The updated conductances must define a connected network. The Greenian
    is returned as an exact SymPy matrix in the original vertex order.
    """
    labels = _connected(net)
    updated = vary_conductances(net, changes)
    if len(updated.components) != 1:
        raise DomainError("finite Greenian update requires connected output")
    columns, increments = _rank_update(labels, changes)
    initial = Matrix(green(net).value.values)
    if not changes:
        return initial, 1
    projection = initial * columns
    system = eye(columns.cols) + increments * columns.T * projection
    return initial - projection * system.inv() * increments * projection.T, system.det()


def vary_greenian_coefficients(net, changes):
    """Update signed Greenian coefficients gamma_ab = -G_ab and invert.

    A singular or indefinite Greenian, or an inverse with negative
    conductances, does not realize an exact connected network.
    """
    labels = _connected(net)
    columns, increments = _rank_update(labels, changes)
    updated = Matrix(green(net).value.values) + columns * increments * columns.T
    if len(labels) > 1 and updated[:-1, :-1].is_positive_definite is not True:
        raise InvalidInput("updated Greenian has no positive definite zero-sum metric")
    gauge = ones(len(labels)) / len(labels)
    inverse = (updated + gauge).inv() - gauge
    return network_from_laplacian(inverse.tolist(), labels, policy=net.policy).value


def vary_resistances(net, changes):
    """Update specified pair resistances and reconstruct the exact network."""
    labels = _connected(net)
    updated = Matrix(resistance(net).value.values)
    for (i, j), delta in zip(_pair_indices(labels, changes), changes.values(), strict=True):
        updated[i, j] += _number(delta, "exact")
        updated[j, i] = updated[i, j]
    return network_from_resistance(updated.tolist(), labels, policy=net.policy).value


def resistance_constraint_polynomials(net, variable_pairs, targets):
    """Return (variables, tree polynomial, residuals) for conductance constraints.

    Each variable is an absolute conductance; other pairs retain their
    original conductances. Targets map pairs to desired resistances. A
    solution needs nonnegative variables, positive tree weight, and zero
    residuals, in the order of the supplied targets.
    """
    labels = _connected(net)
    size = len(labels)
    initial = Matrix(laplacian(net).value.values)
    pairs = tuple(_pair_indices(labels, variable_pairs))
    variables = tuple(Symbol(f"c_{i}", real=True) for i in range(len(pairs)))
    grounded = initial[:-1, :-1]
    for (i, j), variable in zip(pairs, variables, strict=True):
        direction = _direction(size, (i, j))[:-1, :]
        grounded += (variable + initial[i, j]) * direction * direction.T
    tree = expand(grounded.det(method="domain-ge"))
    residuals = []
    for pair, target in zip(_pair_indices(labels, targets), targets.values(), strict=True):
        direction = _direction(size, pair)[:-1, :]
        separating = (grounded + direction * direction.T).det(method="domain-ge") - tree
        residuals.append(expand(separating - _number(target, "exact") * tree))
    return variables, tree, tuple(residuals)
