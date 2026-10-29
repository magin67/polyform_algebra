"""Affine means and terminology, without a chosen metric.

Sources (pmg-notes, branch v5):
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Foundations/PMG%20terminology.md
Sections 'Centroid', 'Affine weight', and 'Earlier names and ambiguous usages'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Projections,%20Bi-coordinates,%20and%20Di-coordinates.md
Section 'Bi-coordinates and frames': epsilon(y)=sum b_i epsilon(h_i).
"""

import pytest

from polyform_algebra import Basis, Element, Point, Vector, Polysimplex, Simplex, Polyform


def test_point_centroid_returns_element_and_frame():
    """The mean of three weight-one points has coefficients 1/3 and weight one."""
    a, b, c = Point.create_list("abc")
    basis = Basis([a, b, c])
    centroid, frame = basis.centroid(name="h")
    assert isinstance(centroid, Element)
    assert centroid.affine_weight == pytest.approx(1)
    assert frame == (a + b + c) / 3
    assert centroid[0] is frame
    assert centroid.basis_for_frame == [basis]
    assert centroid.boundary() == Polysimplex.one()
    assert basis.elements == [a, b, c]


@pytest.mark.parametrize("weights", [(0, 0), (1, 0), (2, -1, 0.5)])
def test_mixed_centroid_weight_is_mean_of_element_weights(weights):
    """Linearity gives epsilon(mean)=mean(epsilon), not an unconditional one."""
    elements = [Element(str(i), affine_weight=w) for i, w in enumerate(weights)]
    centroid, frame = Basis(elements).centroid()
    assert centroid.affine_weight == pytest.approx(sum(weights) / len(weights))
    assert frame.affine_weight() == centroid.affine_weight
    for elem in elements:
        assert frame.terms[Simplex([elem])] == pytest.approx(1 / len(elements))


def test_singleton_centroid_frame_and_unnamed_elements():
    """A singleton mean represents its sole element with coefficient one."""
    a = Point("a")
    centroid, frame = Basis([a]).centroid()
    another, _ = Basis([a]).centroid()
    assert frame == a.as_polysimplex()
    assert centroid.affine_weight == a.affine_weight
    # An unnamed derived element cannot collide with a named basis element.
    assert centroid != a
    assert another != centroid


@pytest.mark.parametrize("elements, error", [
    ([], ValueError), (["a", "b"], TypeError),
    ([Point("a"), Point("a")], ValueError),
])
def test_centroid_rejects_undefined_or_unsupported_input(elements, error):
    """Empty means are undefined; label-only bases need an explicit element adapter."""
    with pytest.raises(error):
        Basis(elements).centroid()


def test_affine_weight_of_boundary_and_signed_frame():
    """Boundary differences have epsilon=0; oriented atomic terms retain their sign."""
    a, b = Point.create_list("ab")
    assert Polysimplex({Simplex([(a, b)]): 1}).affine_weight() == 0
    assert Polysimplex({Simplex([a], sign=-1): 2}).affine_weight() == -2
    assert (2 * a - 3 * b).affine_weight() == -1
    assert Polysimplex.zero().affine_weight() == 0


@pytest.mark.parametrize("frame", [
    Polysimplex.one(), Polysimplex({Simplex(Point.create_list("ab")): 1}),
    Polyform.zero(),
])
def test_affine_weight_does_not_rename_higher_grade_metadata(frame):
    """Epsilon is a functional on linear elements, not on scalars or formal forms."""
    with pytest.raises((TypeError, ValueError)):
        frame.affine_weight()


def test_affine_weight_constructor_and_legacy_spelling_share_storage():
    """Legacy multiplicity remains compatible with equality, hashing and frames."""
    old = Element("p", multiplicity=0.5)
    new = Element("p", affine_weight=0.5)
    assert old == new
    assert hash(old) == hash(new)
    assert old.affine_weight == new.multiplicity == 0.5
    assert Element("p", multiplicity=0.5, affine_weight=0.5) == new
    assert Point("a").affine_weight == 1
    assert Vector("u").affine_weight == 0
    assert Element("default").affine_weight == 1
    with pytest.raises(ValueError, match="must agree"):
        Element("conflict", multiplicity=1, affine_weight=0)


def test_frame_validation_uses_affine_weight():
    """A half-point frame has epsilon=1/2, so declaring weight one must fail."""
    a = Point("a")
    frame = a / 2
    assert Element(frame=frame, affine_weight=0.5)[0] == frame
    assert frame.multiplicity() == frame.affine_weight()  # Legacy atomic case.
    with pytest.raises(ValueError, match="affine weight"):
        Element(frame=frame, affine_weight=1)
