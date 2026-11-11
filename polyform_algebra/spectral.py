"""Exact spectral decompositions and projections of conductance Laplacians."""

from sympy import Matrix, eye, simplify, zeros

from .core import InvalidInput, laplacian
from .toponomes import _require_exact_network


def laplacian_projectors(net):
    """Return one exact projector per distinct Laplacian eigenvalue.

    A projective block includes all directions with the given eigenvalue.
    This symbolic calculation is intended for small exact networks.
    """
    _require_exact_network(net)
    matrix = Matrix(laplacian(net).value.values)
    basis_by_value = matrix.eigenvects()
    result = {}
    for value, multiplicity, vectors in basis_by_value:
        if len(vectors) != multiplicity:
            raise InvalidInput("could not resolve the exact eigenspace")
        basis = Matrix.hstack(*vectors)
        projector = basis * (basis.T * basis).inv() * basis.T
        result[value] = projector.applyfunc(simplify)
    return result


def chebyshev_laplacian(net, step):
    """Evaluate 2(I - T_step(I - L/2)) on an exact Laplacian.

    The result is a symmetric zero-row-sum matrix. Its off-diagonal signs
    can be checked with network_from_laplacian when a graph is expected.
    """
    _require_exact_network(net)
    if type(step) is not int or step < 0:
        raise InvalidInput("Chebyshev step must be a nonnegative integer")
    matrix = Matrix(laplacian(net).value.values)
    identity = eye(matrix.rows)
    if step == 0:
        return zeros(matrix.rows)
    angle_operator = identity - matrix / 2
    previous, current = identity, angle_operator
    for _ in range(1, step):
        previous, current = current, 2 * angle_operator * current - previous
    return 2 * (identity - current)


def fixed_basis_projection(net, modes):
    """Return coordinates and matrix of the diagonal projection in a fixed basis.

    Modes form an orthogonal basis of the zero-sum space in vertex order.
    They need not be normalized. The resulting matrix may have signed
    conductances even when the input network has nonnegative conductances.
    """
    _require_exact_network(net)
    size = len(net.context.labels)
    matrix = Matrix(laplacian(net).value.values)
    vectors = tuple(Matrix(mode) for mode in modes)
    if len(vectors) != max(0, size - 1):
        raise InvalidInput("expected a basis of the zero-sum space")
    projection = zeros(size)
    coordinates = []
    for index, vector in enumerate(vectors):
        if (vector.shape != (size, 1)
                or any(value.is_real is not True or value.is_algebraic is not True or value.is_Float
                       for value in vector)
                or simplify(sum(vector)) != 0
                or any(simplify(vector.dot(other)) != 0 for other in vectors[:index])):
            raise InvalidInput("expected exact real orthogonal zero-sum modes")
        length = simplify(vector.dot(vector))
        if length == 0:
            raise InvalidInput("basis modes must be nonzero")
        coordinate = simplify(vector.dot(matrix * vector) / length)
        coordinates.append(coordinate)
        projection += coordinate * vector * vector.T / length
    return tuple(coordinates), projection.applyfunc(simplify)
