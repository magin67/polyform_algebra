import pytest

from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.simplex import Simplex

# ---------- Points ----------
a, b, c, d, e, f, g, h, r, x, y, z = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'r', 'x', 'y', 'z'])

# ---------- Point ----------

def test_point_mul():
    """Two distinct points multiply to their oriented atomic simplex."""
    smp = a * b
    assert smp == Simplex([a, b])

@pytest.mark.parametrize("element_type", [Point, Vector])
@pytest.mark.parametrize("same_instance", [True, False])
def test_equal_elements_add_coefficients(element_type, same_instance):
    """Equal elements add coefficients: u+u=2u, for identical and distinct instances.

    A dictionary literal previously discarded one of the equal summands."""
    left = element_type("u")
    right = left if same_instance else element_type("u")

    result = left + right

    assert result.terms == {Simplex([left]): 2}
    assert result == 2 * left


@pytest.mark.parametrize("element_type", [Point, Vector])
def test_distinct_elements_add_preserves_terms(element_type):
    """Fixing u+u must preserve both distinct summands in u+v."""
    left, right = element_type.create_list(["u", "v"])

    assert (left + right).terms == {
        Simplex([left]): 1,
        Simplex([right]): 1,
    }


# ---------- Vector ----------

def test_vector_comps():
    """The frame b-a+2(d-c) retains all four signed coordinates."""
    v = Vector('v', frame=b - a - 2*c + 2*d)
    # b - a + 2(d - c) = -a + b -2c + 2d
    expected_terms = {Simplex([a]): -1, Simplex([b]): 1, Simplex([c]): -2, Simplex([d]): 2}
    assert v[0].terms == expected_terms

def test_vector_from_vertices():
    """A difference of points has coordinates +1 and -1 and affine weight zero."""
    v = Vector('v', frame=(a - b))
    assert v[0].terms == {Simplex([a]): 1, Simplex([b]): -1}

def test_vector_lists():
    """Coordinate lists expose the atoms and coefficients of a stored vector frame."""
    v = Vector('v', frame=(a + b - 2*c))
    els, coeffs = v[0].terms_lists()
    assert set(els) == set([a, b, c])
    assert set(coeffs) == set([1, 1, -2])

def test_vector_bilinear_form():
    """Polarization of b-a and c-b involves the three pairwise edge squares."""
    v = Vector('v', frame=(b - a))
    w = Vector('w', frame=(c - b))
    qf = v[0] @ w[0]
    # Expected edge coefficients: ab=-0.5, ac=0.5, bc=-0.5.
    assert len(qf.terms) == 3
    forms = {frozenset(k.components[0]) for k in qf.terms}
    expected = {frozenset({a, b}), frozenset({a, c}), frozenset({b, c})}
    assert forms == expected
    # Their coefficient sum is -0.5.

def test_vector_to_polyform():
    """The scalar polyform of b-a with itself is exactly the edge square [(ab)]^2."""
    v = Vector('v', frame=(b - a))
    qf = v[0] @ v[0]
    assert len(qf.terms) == 1
    form = list(qf.terms.keys())[0]
    # Orientation disappears in a quadratic form.
    assert set(form.components[0]) == {a, b}
    assert qf.terms[form] == 1
