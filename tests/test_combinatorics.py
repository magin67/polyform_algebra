"""PMG Combinatorial Potentials (v5): sizes, cycle types, and mixed states."""

from math import factorial

import pytest

from polyform_algebra.combinatorics import (
    first_kind_potential,
    path_cycle_potential,
    second_kind_potential,
)


def test_second_kind_normalizes_equal_blocks():
    """Swapping equal blocks leaves the same partition (Partitions, one type)."""
    assert second_kind_potential((2, 2)) == 3
    assert second_kind_potential(iter((3, 2, 3))) == 280
    assert second_kind_potential((30, 30)) == factorial(60) // (2 * factorial(30)**2)


def test_first_kind_adds_cyclic_orders():
    """Each triple has two cyclic orders (Cyclic Partitions, cycle type)."""
    assert first_kind_potential((3, 3, 2)) == 4 * second_kind_potential((3, 3, 2)) == 1120
    assert first_kind_potential((2, 1, 2)) == path_cycle_potential((2, 2, 1)) == 15


def test_mixed_state_keeps_component_kinds_distinct():
    """A loop and an open arc have separate multiplicities (Mixed States, examples)."""
    assert path_cycle_potential((1,), (1,)) == 2
    assert path_cycle_potential((1, 1), ()) == 1
    assert path_cycle_potential((), (1, 1)) == 1
    assert path_cycle_potential((2,), (2,)) == 12
    assert path_cycle_potential((), (2,)) == 2


def test_empty_types():
    """The empty type contributes one object in all three models."""
    assert second_kind_potential(()) == first_kind_potential(()) == path_cycle_potential() == 1


@pytest.mark.parametrize("potential", [second_kind_potential, first_kind_potential,
                                      lambda lengths: path_cycle_potential(path_lengths=lengths)])
def test_component_lengths_require_positive_integers(potential):
    """Zero-sized blocks and zero-arc components are outside the input convention."""
    for value in (0, -1):
        with pytest.raises(ValueError):
            potential((value,))
    for value in (True, 1.5, "2"):
        with pytest.raises(TypeError):
            potential((value,))
