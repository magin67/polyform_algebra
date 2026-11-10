"""Partition Toponomes in PMG v5.

https://magin67.github.io/pmg-notes/en/polyform-algebra-and-combinatorics/partition-toponomes/
"""

from fractions import Fraction
from itertools import combinations

import pytest

from polyform_algebra.core import network
from polyform_algebra.toponomes import (
    coefficient_table, color_class_toponome, forest_toponome,
    partition_product, partition_toponome,
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
