"""Checks for partition forms from the PMG tetrahedron example.

https://github.com/magin67/pmg-notes/blob/c1c017893e96e221061b1e33bfc034891bcd0518/pmg-library/en/Polyform%20Algebra%20and%20Combinatorics/Structural%20identities/Partition%20Forms%20and%20the%20Tetrahedron%20Identity.md
"""

from fractions import Fraction

import pytest

from polyform_algebra.core import (
    AlgebraLimits, Boundary, Context, InvalidInput, ResourceLimit, UnsupportedOperation,
    equivalent, form_coefficients, form_product, partition_form, square,
)


def test_partition_form_ignores_order_and_singletons():
    """Exterior signs disappear on squaring and singleton boundaries equal one."""
    context = Context(('a', 'b', 'c', 'd', 'e'))
    q = partition_form(context, (('a', 'b'), ('c', 'd'))).value
    reordered = partition_form(context, (('d', 'c'), ('e',), ('b', 'a'))).value
    assert q == reordered
    edge_ab = square(Boundary(context, (0, 1))).value
    edge_cd = square(Boundary(context, (2, 3))).value
    assert equivalent(q, form_product(edge_ab, edge_cd).value).value
    unit = partition_form(context, (('a',),)).value
    assert form_coefficients(unit).value == ((((), ()), Fraction(1)),)


def test_form_coefficients_are_sparse_and_ordered():
    """The square of b-a has two diagonal coefficients and both cross terms."""
    context = Context(('a', 'b'))
    edge = partition_form(context, (('a', 'b'),)).value
    assert form_coefficients(edge).value == (
        (((0,), (0,)), Fraction(1)), (((0,), (1,)), Fraction(-1)),
        (((1,), (0,)), Fraction(-1)), (((1,), (1,)), Fraction(1)),
    )
    assert form_coefficients(edge - edge).value == ()
    with pytest.raises(ResourceLimit):
        form_coefficients(edge, limits=AlgebraLimits(max_terms=3))


@pytest.mark.parametrize('blocks', [((),), (('b',), ('a', 'b')), (('c',),), (('a', 'a'),)])
def test_partition_blocks_must_be_valid(blocks):
    """Only known vertices in disjoint nonempty blocks define a partition."""
    with pytest.raises(InvalidInput):
        partition_form(Context(('a', 'b')), blocks)


def test_partition_form_requires_unit_affine_weights():
    """The omitted-singleton convention assumes each singleton boundary is one."""
    context = Context(('a', 'b'), epsilon=(1, 0))
    with pytest.raises(UnsupportedOperation, match='unit affine weights'):
        partition_form(context, (('a',),))
