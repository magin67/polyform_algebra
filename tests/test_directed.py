"""Examples from the Balanced Graphs series in PMG v5."""

import pytest
from sympy import Matrix, Rational, eye, ones, zeros

from polyform_algebra.core import DomainError, InvalidInput
from polyform_algebra.directed import DirectedNetwork


def test_cycle_with_parallel_arcs():
    """Parallel arcs add, and each principal cofactor counts trees toward its root."""
    graph = DirectedNetwork(tuple("abc"), (("a", "b", 1), ("a", "b", 1),
                                          ("b", "a", 1), ("b", "c", 1), ("c", "a", 1)))
    assert graph.laplacian() == Matrix([[2, -2, 0], [-1, 2, -1], [-1, 0, 1]])
    assert graph.tree_potentials() == (2, 2, 2)
    assert graph.laplacian() * graph.centered_greenian() == eye(3) - ones(3) / 3
    assert graph.one_way_resistances() == Matrix([[0, 1, 4], [2, 0, 3], [2, 3, 0]]) / 6
    assert graph.vertex_gram() == Matrix([[0, -4, -5], [-2, 0, -7], [-7, -5, 0]]) / 12


def test_exact_balancing_and_domain():
    """Mean-one scaling of the three ratings rows preserves their proportions."""
    graph = DirectedNetwork(tuple("abc"), (("a", "b", 1), ("a", "c", 2), ("b", "a", 2),
                                          ("b", "c", 1), ("c", "a", 1), ("c", "b", 1)))
    balanced, factors = graph.balance()
    assert factors == tuple(Rational(value, 16) for value in (15, 12, 21))
    assert ones(1, 3) * balanced.laplacian() == zeros(1, 3)
    assert balanced.tree_potentials() == (Rational(315, 64),) * 3
    with pytest.raises(DomainError):
        graph.centered_greenian()
    with pytest.raises(DomainError):
        DirectedNetwork(tuple("ab"), (("a", "b", 1),)).balance()


def test_directed_forests_require_one_successor_and_no_cycle():
    """The three-cycle contributes three trees and excludes the full cycle."""
    graph = DirectedNetwork(tuple("abc"), (("a", "b", 1), ("b", "c", 1), ("c", "a", 1)))
    assert graph.forest_toponome() == {(): 1, (1,): 3, (2,): 3}
    fork = DirectedNetwork(tuple("abc"), (("a", "b", 1), ("a", "c", 2)))
    assert fork.forest_toponome() == {(): 1, (1,): 3}
    assert DirectedNetwork(tuple("ab"), (("a", "b", 2), ("b", "a", 3))).forest_toponome() == {
        (): 1, (1,): 5,
    }


def test_isolated_vertex_and_invalid_arcs():
    """The singleton has one empty tree, and invalid exact weights are rejected."""
    graph = DirectedNetwork(("a",), ())
    assert graph.tree_potentials() == (1,)
    assert graph.centered_greenian() == zeros(1)
    for arc in (("a", "a", 1), ("a", "c", 1), ("a", "b", -1), ("a", "b", 0.5)):
        with pytest.raises(InvalidInput):
            DirectedNetwork(tuple("ab"), (arc,))
