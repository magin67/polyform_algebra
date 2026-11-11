"""Exact forest and spectral data at marked vertices of a network."""

from collections import Counter

from sympy import Dummy, Matrix, Poly, Symbol, expand

from .core import InvalidInput, laplacian
from .core.data import _labels, _number
from .toponomes import _require_exact_network, _spanning_forest_components


def _boundary_positions(net, boundary):
    _require_exact_network(net)
    labels = _labels(boundary)
    positions = {label: i for i, label in enumerate(net.context.labels)}
    if not labels or not set(labels) <= positions.keys():
        raise InvalidInput("boundary must contain known vertex labels")
    return labels, tuple(positions[label] for label in labels)


def forest_passport(net, boundary):
    """Map each boundary partition to its nonzero ordinary coefficients by grade.

    All spanning forests contribute, including those with internal components.
    The exponential enumeration is intended for small networks.
    """
    labels, positions = _boundary_positions(net, boundary)
    vertex_index = {label: i for i, label in enumerate(net.context.labels)}
    edges = [(vertex_index[a], vertex_index[b], weight) for a, b, weight in net.edges]
    result = {}
    for components, grade, weight in _spanning_forest_components(len(vertex_index), edges):
        blocks = {}
        for label, position in zip(labels, positions, strict=True):
            blocks.setdefault(components[position], []).append(label)
        partition = tuple(tuple(block) for block in blocks.values())
        result.setdefault(partition, Counter())[grade] += weight
    return {partition: dict(sorted(poly.items())) for partition, poly in result.items()}


def join_forest_states(boundary, first, second):
    """Return the joined partition, or None if the two states create a cycle.

    Both partitions must cover the same ordered boundary. The parts represent
    edge-disjoint components that intersect only in the boundary vertices.
    """
    labels = _labels(boundary)
    components = {label: label for label in labels}
    for partition in (first, second):
        if set(_labels(tuple(label for block in partition for label in block))) != set(labels):
            raise InvalidInput("each forest state must partition the boundary")
        if any(not block for block in partition):
            raise InvalidInput("forest state blocks must be nonempty")

    # A tree on each block represents that block's independent connections.
    for partition in (first, second):
        for anchor, *rest in partition:
            for label in rest:
                origin, target = components[anchor], components[label]
                if origin == target:
                    return None
                for vertex in labels:
                    if components[vertex] == target:
                        components[vertex] = origin
    blocks = {}
    for label in labels:
        blocks.setdefault(components[label], []).append(label)
    return tuple(tuple(block) for block in blocks.values())


def pair_forest_passports(boundary, first, second):
    """Compose passports for edge-disjoint parts sharing exactly the boundary."""
    result = {}
    for first_state, first_poly in first.items():
        for second_state, second_poly in second.items():
            state = join_forest_states(boundary, first_state, second_state)
            if state is not None:
                poly = result.setdefault(state, Counter())
                for i, a in first_poly.items():
                    for j, b in second_poly.items():
                        poly[i + j] += a * b
    return {partition: dict(sorted(poly.items())) for partition, poly in result.items()}


def spectral_pair_passport(net, boundary, variable=Symbol("x")):
    """Return determinant, two diagonal cofactors, cross cofactor and interior minor.

    The pencil is xI-L. The order of the two diagonal cofactors follows
    the order of the supplied boundary vertices.
    """
    _, positions = _boundary_positions(net, boundary)
    if len(positions) != 2:
        raise InvalidInput("spectral pair passport requires two boundary vertices")
    return _pencil_passport(_laplacian_pencil(net, variable), positions)


def _laplacian_pencil(net, variable):
    matrix = Matrix(laplacian(net).value.values)
    return variable * Matrix.eye(matrix.rows) - matrix


def _pencil_passport(pencil, positions):
    a, b = positions
    indices = list(range(pencil.rows))
    without_a = [i for i in indices if i != a]
    without_b = [i for i in indices if i != b]
    interior = [i for i in indices if i not in positions]

    def minor(rows, columns):
        return expand(pencil.extract(rows, columns).det(method="domain-ge"))

    return (
        minor(indices, indices), minor(without_a, without_a),
        minor(without_b, without_b), (-1)**(a+b) * minor(without_b, without_a),
        minor(interior, interior),
    )


def spectral_pair_transfer(kernel, entry, exit, variable=Symbol("x"), *, shared_conductance=0):
    """Return the five-coordinate transfer for attaching a two-terminal kernel.

    Entry labels identify the previous graph's ordered boundary. Exit labels
    may overlap the entry or be disjoint. Supply the conductance of a shared
    terminal edge to count that edge once; otherwise the conductances add.
    Polynomial minors define the map even at zeros of the old interior minor.
    """
    _, entry_positions = _boundary_positions(kernel, entry)
    _, exit_positions = _boundary_positions(kernel, exit)
    if len(entry_positions) != 2 or len(exit_positions) != 2:
        raise InvalidInput("transfer requires two entry and two exit vertices")
    shared = _number(shared_conductance, "exact")
    if shared < 0:
        raise InvalidInput("shared conductance must be nonnegative")

    first, second = entry_positions
    pencil = _laplacian_pencil(kernel, variable)
    u, v, w = Dummy("u"), Dummy("v"), Dummy("w")
    pencil[first, first] += shared - variable + u
    pencil[second, second] += shared - variable + w
    pencil[first, second] += v - shared
    pencil[second, first] += v - shared
    rows = []
    for minor in _pencil_passport(pencil, exit_positions):
        coefficients = Poly(minor, u, v, w)
        coefficient = coefficients.coeff_monomial
        # Replace the old Schur entries by C/E, -D/E, B/E and its
        # determinant by A/E using Jacobi's identity.
        rows.append([coefficient(u*w), coefficient(w), coefficient(u),
                     -coefficient(v), coefficient(1)])
    return Matrix(rows)


def pair_spectral_passports(first, second, variable=Symbol("x"), *, shared_conductance=0):
    """Glue spectral pair passports with the same ordered terminals.

    A physical terminal edge present in both inputs can be counted once by
    supplying its conductance. With the default, edge contributions add.
    """
    a, b, c, d, e = first
    aa, bb, cc, dd, ee = second
    common = _number(shared_conductance, "exact")
    if common < 0:
        raise InvalidInput("shared conductance must be nonnegative")
    interior = e * ee
    left, right, cross = b * ee + e * bb, c * ee + e * cc, d * ee + e * dd
    full = a * ee + e * aa + b * cc + c * bb - 2 * d * dd
    shift = common - variable
    return tuple(map(expand, (
        full + shift * (left + right) - 2 * common * cross
        + (variable**2 - 2 * common * variable) * interior,
        left + shift * interior, right + shift * interior,
        cross + common * interior, interior,
    )))
