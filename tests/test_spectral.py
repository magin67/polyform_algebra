"""Spectral decomposition examples from PMG notes v5 at 0f478b3.

https://github.com/magin67/pmg-notes/tree/0f478b36dcc35c20a06c8d6c3041493996742910/pmg-library/ru/Инварианты%2C%20разложения%20и%20спектры
"""

import pytest
from sympy import Matrix, eye, zeros

from polyform_algebra.core import InvalidInput, green, laplacian, network
from polyform_algebra.spectral import chebyshev_laplacian, fixed_basis_projection, laplacian_projectors


def test_repeated_eigenvalue_is_one_projector():
    graph = network(tuple("abcd"), [("a", "b", 1), ("b", "c", 1), ("c", "d", 1), ("d", "a", 1)])
    blocks = laplacian_projectors(graph)
    assert {value: block.trace() for value, block in blocks.items()} == {0: 1, 2: 2, 4: 1}
    for value, block in blocks.items():
        assert block == block.T == block**2
        assert all(block * other == zeros(4) for key, other in blocks.items() if key != value)
    assert sum(blocks.values(), zeros(4)) == eye(4)
    assert 2*blocks[2] + 4*blocks[4] == Matrix(laplacian(graph).value.values)
    assert blocks[2]/2 + blocks[4]/4 == Matrix(green(graph).value.values)


def test_chebyshev_doubling():
    labels = tuple("abcdef")
    path = network(labels, [(a, b, 1) for a, b in zip(labels, labels[1:])])
    order = [0, 1, 3, 5, 4, 2]
    expected = Matrix(6, 6, lambda i, j: 2 if i == j else -int((i-j) % 6 in (1, 5)))
    assert chebyshev_laplacian(path, 2).extract(order, order) == expected
    assert chebyshev_laplacian(path, 0) == zeros(6)
    with pytest.raises(InvalidInput):
        chebyshev_laplacian(path, -1)


def test_fixed_basis_projection_validates_modes():
    path = network(tuple("abc"), [("a", "b", 1), ("b", "c", 1)])
    coordinates, result = fixed_basis_projection(path, [(1, 0, -1), (1, -2, 1)])
    assert coordinates == (1, 3)
    assert result == Matrix(laplacian(path).value.values)
    with pytest.raises(InvalidInput):
        fixed_basis_projection(path, [(1, 0, -1), (1, -1, 0)])
