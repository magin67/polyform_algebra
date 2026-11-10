"""Checks of the isotropic extension against the PMG Spaces example.

https://github.com/magin67/pmg-notes/blob/6a80e80026413c7842e2d37de508b67c7ad7b687/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Nonlocal%20Elements%2C%20Spheres%2C%20and%20Their%20Centers.md#example-nonlocal-vertices-and-an-unchanged-boundary-norm
"""

from fractions import Fraction as Q

import pytest

from polyform_algebra.core import (
    ArithmeticMismatch, ContextMismatch, DomainError, InvalidInput, IsotropicExtension,
    UnsupportedOperation, coordinates, network, to_numeric,
)


@pytest.fixture
def unit_path():
    return network(('a', 'b', 'c'), (('a', 'b', 1), ('b', 'c', 1)))


def test_nonlocal_vertices_and_unchanged_boundary_norm(unit_path):
    """Shifts by z change element norms by twice the shift and preserve the edge Gram matrix."""
    extension = IsotropicExtension(unit_path)
    a, b, c = (extension.vertex(label) for label in unit_path.context.labels)
    z = extension.z
    a_shifted, b_shifted, c_shifted = a + z, b - 2*z, c + 3*z

    assert extension.context.epsilon == (1, 1, 1, 0)
    assert extension.norm(z) == extension.affine_weight(z) == 0
    assert extension.inner_product(a, z) == 1
    assert extension.inner_product(a, b) == Q(-1, 2)
    assert extension.inner_product(a_shifted, b_shifted) == Q(-3, 2)
    assert tuple(extension.norm(s) for s in (a_shifted, b_shifted, c_shifted)) == (2, -4, 6)
    assert tuple(extension.center(s) for s in (a_shifted, b_shifted, c_shifted)) == (a, b, c)

    u, v = b_shifted - a_shifted, c_shifted - a_shifted
    assert (u, v) != (b - a, c - a)
    assert extension.norm(u) == extension.norm(b - a) == 1
    assert extension.norm(v) == extension.norm(c - a) == 2
    assert extension.inner_product(u, v) == extension.inner_product(b - a, c - a) == 1
    assert extension.norm(u) * extension.norm(v) - extension.inner_product(u, v)**2 == 1


def test_center_of_an_affine_combination(unit_path):
    """The midpoint has norm -R_ac/4 and needs a correction of R_ac*z/8."""
    extension = IsotropicExtension(unit_path, null_label='n')
    a, c = extension.vertex('a'), extension.vertex('c')
    midpoint = (a + c) / 2
    assert extension.context.labels == ('a', 'b', 'c', 'n')
    assert extension.norm(midpoint) == Q(-1, 2)
    center = extension.center(midpoint)
    assert center == midpoint + extension.z / 4
    assert extension.norm(center) == 0 and extension.affine_weight(center) == 1
    assert extension.center(center + Q(5, 3) * extension.z) == center


def test_isotropic_extension_domain_and_context(unit_path):
    """The resistance construction needs a connected exact network and matching coordinates."""
    extension = IsotropicExtension(unit_path)
    with pytest.raises(DomainError, match='connected'):
        IsotropicExtension(network(('a', 'b'), ()))
    with pytest.raises(UnsupportedOperation, match='exact arithmetic'):
        IsotropicExtension(to_numeric(unit_path).value)
    with pytest.raises(InvalidInput, match='unique'):
        IsotropicExtension(unit_path, null_label='a')
    with pytest.raises(InvalidInput, match='vertex label'):
        extension.vertex('z')
    with pytest.raises(ContextMismatch):
        extension.norm(IsotropicExtension(unit_path).vertex('a'))
    with pytest.raises(ArithmeticMismatch):
        extension.norm(to_numeric(extension.vertex('a')).value)
    with pytest.raises(DomainError, match='affine weight one'):
        extension.center(coordinates(extension.context, (2, 0, 0, 0)))
