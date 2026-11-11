"""Exact Greenian forms and complementary affine coefficients."""

from fractions import Fraction
from itertools import combinations

from .core.algebra import (
    Boundary, Exterior, Form, _accumulate, _answer, _closed_pairs, _forms, _limits,
    form_product, partition_form, polar_form,
)
from .core.data import DomainError, InvalidInput, _labels, _ordered
from .core.network_math import green


def greenian_form(net, *, limits=None):
    """Represent the centered Greenian by signed pair-boundary squares."""
    values = green(net).value.values
    limits = _limits(limits)
    terms = (
        (Boundary(net.context, (i, j), limits=limits), -values[i][j])
        for i, j in combinations(range(len(values)), 2) if values[i][j]
    )
    return _answer(Form(net.context, tuple(terms), limits), "greenian_form", limits)


def cut_partition_form(context, blocks, *, limits=None):
    """Square the wedge of centered block sums, with the last block omitted.

    Blocks must cover the whole vertex context. Omitting a different block
    leaves the same form because the sum of all centered block sums is zero.
    """
    unit = partition_form(context, (), limits=limits).value
    blocks = tuple(_labels(block) for block in _ordered(blocks))
    labels = context.labels
    flat = tuple(label for block in blocks for label in block)
    if not labels or not all(blocks) or len(flat) != len(labels) or set(flat) != set(labels):
        raise InvalidInput("cut blocks must be a full partition of a nonempty context")

    result = unit
    for block in blocks[:-1]:
        mean = Fraction(len(block), len(labels))
        coefficients = tuple(((i,), int(label in block) - mean) for i, label in enumerate(labels))
        centered = Exterior(context, coefficients, limits=unit.limits)
        factor = polar_form(centered, centered, limits=limits).value * Fraction(1, 2)
        result = form_product(result, factor, limits=limits).value
    return _answer(result, "cut_partition_form", result.limits)


def affine_coefficients(form, *, limits=None):
    """Return coefficients in the difference basis based at the last vertex.

    Index tuples refer to the remaining vertices, in context order. Each
    tensor slot must lie in the boundary subalgebra.
    """
    _forms(form)
    limits = _limits(limits, form)
    pairs = _closed_pairs(form, limits, vertex_space=True)
    last = len(form.context.labels) - 1
    return {pair: weight for pair, weight in pairs.items()
            if last not in pair[0] and last not in pair[1]}


def dual_coefficients(form, top, *, limits=None):
    """Delete both blades of each term from a pure top form.

    Return the resulting coefficient mapping in the same affine basis.
    Deletion uses the natural vertex-coordinate pairing with centered dual
    elements. The two exterior signs and the top coefficient are retained.
    """
    _forms(form, top)
    limits = _limits(limits, form, top)
    blade = tuple(range(len(form.context.labels) - 1))
    top_pairs = affine_coefficients(top, limits=limits)
    if set(top_pairs) != {(blade, blade)}:
        raise DomainError("expected a nonzero pure top form")
    reference = _closed_pairs(top, limits, vertex_space=True)
    terms = _closed_pairs(form, limits, vertex_space=True)
    last = len(blade)

    def delete(full, removed):
        positions = {i: position for position, i in enumerate(full)}
        if any(i not in positions for i in removed):
            return None
        remaining = tuple(i for i in full if i not in removed)
        if last in remaining:
            return None
        inversions = sum(positions[i] - position for position, i in enumerate(removed))
        return remaining, (-1) ** inversions

    result = {}
    for (left, right), weight in reference.items():
        for (a, b), value in terms.items():
            dl, dr = delete(left, a), delete(right, b)
            if dl is not None and dr is not None:
                _accumulate(result, (dl[0], dr[0]), weight * value * dl[1] * dr[1], limits)
    return result
