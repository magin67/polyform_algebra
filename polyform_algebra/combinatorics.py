"""Exact combinatorial potentials of partitions and directed states."""

from collections import Counter
from math import factorial, prod
from operator import index


def _lengths(values):
    lengths = []
    for value in values:
        if isinstance(value, bool):
            raise TypeError("component lengths must be integers")
        length = index(value)
        if length <= 0:
            raise ValueError("component lengths must be positive")
        lengths.append(length)
    return lengths


def second_kind_potential(block_sizes):
    """Count set partitions of distinguishable labels with given block sizes.

    Block sizes are positive integers in any order. The empty type counts
    the unique empty partition.
    """
    sizes = _lengths(block_sizes)
    denominator = prod(map(factorial, sizes))
    denominator *= prod(map(factorial, Counter(sizes).values()))
    return factorial(sum(sizes)) // denominator


def path_cycle_potential(cycle_lengths=(), path_lengths=()):
    """Count states with distinct arc labels and indistinguishable vertices.

    The arguments contain positive arc counts. Length one is valid for a
    loop or a one-arc path. Isolated vertices are omitted. Equal lengths in
    different kinds of components have separate multiplicities.
    """
    cycles = _lengths(cycle_lengths)
    paths = _lengths(path_lengths)
    denominator = prod(cycles) * prod(map(factorial, Counter(cycles).values()))
    denominator *= prod(map(factorial, Counter(paths).values()))
    return factorial(sum(cycles) + sum(paths)) // denominator


def first_kind_potential(cycle_lengths):
    """Count permutations of distinguishable labels with given cycle lengths."""
    return path_cycle_potential(cycle_lengths)
