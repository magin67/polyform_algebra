"""Partition Toponomes in PMG v5.

https://magin67.github.io/pmg-notes/en/polyform-algebra-and-combinatorics/partition-toponomes/
"""

from fractions import Fraction
from itertools import combinations

import pytest

from polyform_algebra.core import network
from polyform_algebra.toponomes import (
    coefficient_table, color_class_toponome, forest_toponome,
    laplacian_polynomial, partition_product, partition_toponome,
    rooted_forest_weight, vertex_toponome,
)


def test_chosen_representation_collects_signed_terms():
    """Permuting block sizes and adding singletons preserves the positive type."""
    assert partition_toponome([
        ((1, 3, 2), Fraction(3, 2)),
        ((2, 3), -1),
        ((1, 1), 1),
        ((4,), 2),
        ((1, 4), -2),
    ]) == {(2, 1): Fraction(1, 2), (): 1}


def test_forest_toponome_skips_cycles_and_retains_isolates():
    """A weighted triangle and a separate edge give independent forest choices."""
    graph = network(tuple("abcdef"), (
        ("a", "b", Fraction(1, 2)), ("a", "c", 2), ("b", "c", 3), ("d", "e", 4),
    ))
    assert forest_toponome(graph) == {
        (): 1, (1,): Fraction(19, 2), (2,): Fraction(17, 2),
        (1, 1): 22, (2, 1): 34,
    }
    assert forest_toponome(network(("a", "b"), ())) == {(): 1}


def test_partial_matchings_and_bilinearity():
    """Two pairs yield one empty, four single and two double matchings."""
    pair = {(1, 1): 1}
    expected = {(1, 1, 1, 1): 1, (2, 1, 1): 4, (2, 2): 2}
    assert partition_product(pair, pair) == expected
    assert partition_product({(): 2, (1, 1): -1}, pair) == {
        (1, 1): 2, **{kind: -value for kind, value in expected.items()},
    }
    assert partition_product({(): 1}, {(1, 2): 3, (2, 1): -3}) == {}
    assert coefficient_table(expected) == {(4, 4): 1, (4, 3): 4, (4, 2): 2}


def test_color_classes_are_unordered_and_weighted():
    """The two disconnected edges have seven class partitions without color names."""
    graph = network(tuple("abcd"), (("a", "b", 5), ("c", "d", 7)))
    assert color_class_toponome(graph) == partition_product({(1, 1): 1}, {(1, 1): 1})
    triangle = network(tuple("abc"), tuple((a, b, 1) for a, b in combinations("abc", 2)))
    assert color_class_toponome(triangle, {"a": 2, "b": 3, "c": 2}) == {(3, 2, 2): 1}
    assert color_class_toponome(network(tuple("ab"), ())) == {(1, 1): 1, (2,): 1}


def test_exact_input_validation():
    """Float forest conductances and malformed sizes cannot enter exact counts."""
    with pytest.raises(TypeError):
        forest_toponome(network(tuple("ab"), (("a", "b", 1),), arithmetic="numeric"))
    with pytest.raises(ValueError):
        partition_toponome([((2, 0), 1)])
    with pytest.raises(ValueError):
        color_class_toponome(network(tuple("ab"), ()), {"a": 1})


def test_root_weights_precede_grouping():
    """P6 cell (4, 2) combines size products 4 * 2 and 3 * 3."""
    path = network(tuple("abcdef"), tuple((a, b, 1) for a, b in zip("abcde", "bcdef")))
    forests = forest_toponome(path)
    assert coefficient_table(forests)[4, 2] == 3
    assert coefficient_table(vertex_toponome(forests))[4, 2] == 25
    assert vertex_toponome({(): 1, (1,): Fraction(-1, 3)}) == {(): 1, (1,): Fraction(-2, 3)}


def test_laplacian_polynomial_preserves_exact_weights_and_degree():
    """Two weighted edges and an isolate give (1 + t)(1 + 6t), of degree two."""
    graph = network(tuple("abcde"), (("a", "b", Fraction(1, 2)), ("c", "d", 3)))
    assert laplacian_polynomial(graph) == (1, 7, 6)
    assert laplacian_polynomial(network(tuple("ab"), ())) == (1,)
    with pytest.raises(TypeError):
        laplacian_polynomial(network(tuple("ab"), (), arithmetic="numeric"))


def test_prescribed_roots_use_the_original_laplacian():
    """Opposite roots in C4 leave diagonal entries 2 in the principal minor."""
    cycle = network(tuple("abcd"), (("a", "b", 1), ("b", "c", 1), ("c", "d", 1), ("a", "d", 1)))
    assert rooted_forest_weight(cycle, ("a", "c")) == 4
    assert rooted_forest_weight(cycle, ("a", "b")) == 3
    assert rooted_forest_weight(cycle, tuple("abcd")) == 1
    assert rooted_forest_weight(cycle, ()) == 0
    with pytest.raises(ValueError):
        rooted_forest_weight(cycle, ("z",))
    with pytest.raises(ValueError):
        rooted_forest_weight(cycle, ("a", "a"))
