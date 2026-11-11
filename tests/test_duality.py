"""Finite exact checks from Greenian and Polyform Duality (PMG v5)."""

from fractions import Fraction

import pytest

from polyform_algebra.core import (
    Context, equivalent, form_exp, laplacian_form, network, partition_form, spanning_coefficient,
)
from polyform_algebra.duality import affine_coefficients, cut_partition_form, dual_coefficients, greenian_form


def test_weighted_path_mutual_duality():
    """Complementary grades of the weighted P4 exponentials agree in both directions."""
    graph = network(tuple("abcd"), (("a", "b", 2), ("b", "c", 3), ("c", "d", 5)))
    direct = form_exp(laplacian_form(graph).value).value
    inverse = form_exp(greenian_form(graph).value).value
    top = partition_form(graph.context, (tuple("abcd"),)).value
    assert spanning_coefficient(inverse).value == Fraction(1, 480)
    assert dual_coefficients(direct, top * Fraction(1, 480)) == affine_coefficients(inverse)
    assert dual_coefficients(inverse, 30 * top) == affine_coefficients(direct)


def test_cut_duality_and_omitted_block():
    """A partition into three blocks gives the same cut form with any omitted block."""
    context = Context(tuple("abcd"))
    blocks = (("a", "b"), ("c",), ("d",))
    top = partition_form(context, (context.labels,)).value
    direct = partition_form(context, blocks).value
    cut = cut_partition_form(context, blocks).value
    assert equivalent(cut, cut_partition_form(context, blocks[::-1]).value).value
    assert dual_coefficients(direct, top * Fraction(1, 16)) == affine_coefficients(cut)
    assert dual_coefficients(cut, top) == affine_coefficients(direct)
    with pytest.raises(ValueError):
        cut_partition_form(context, (("a",), ("c",)))
    with pytest.raises(ValueError):
        dual_coefficients(direct, direct)
