import pytest
import sys
import os

from polyform_algebra.objects.element import Point
from polyform_algebra.objects.simplex import Simplex
from polyform_algebra.combinations.polysimplex import Polysimplex

# ---------- Points ----------
a, b, c, d, e, f, g, h, r, x, y, z = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'r', 'x', 'y', 'z'])

# ---------- Simplex ----------

def test_simplex_grade():
    """Two atomic points have grade two; their boundary has grade one."""
    s = Simplex([a, b])
    assert s.grade() == 2   # number of atomic factors
    f = Simplex([(a, b)])
    assert f.grade() == 1   # len(comp)-1

def test_simplex_creation():
    """An atomic pair retains two components and is a nonzero grade-two simplex."""
    s = Simplex([a, b])
    assert s.size() == 2
    assert s.grade() == 2
    assert not s.is_zero()
    assert s.components == (a, b)

def test_simplex_mul():
    """Disjoint atomic products concatenate; a repeated point makes the wedge zero."""
    s1 = Simplex([a, b])
    s2 = Simplex([c, d])
    s3 = s1 * s2
    assert s3.components == (a, b, c, d)
    # Repeating b makes the exterior product zero.
    s4 = Simplex([a, b]) * Simplex([b, c])
    assert s4.is_zero()

def test_simplex_boundary():
    """Boundaries are alternating sums: d(ab)=b-a and d(abc)=bc-ac+ab."""
    s = Simplex([a, b])
    bnd = s.boundary()
    assert isinstance(bnd, Polysimplex)
    # Boundary [a,b] = [b] - [a].
    assert bnd.terms == {Simplex([b]): 1, Simplex([a]): -1}
    s2 = Simplex([a, b, c])
    bnd2 = s2.boundary()
    # [b,c] - [a,c] + [a,b]
    expected = {Simplex([b,c]): 1, Simplex([a,c]): -1, Simplex([a,b]): 1}
    assert bnd2.terms == expected

# ---------- PolySimplex ----------

def test_polysimplex_eq():
    """Reversing boundary orientation and its coefficient preserves the object."""
    s1 = Polysimplex({(a, b, c): 1, (d, e): 1})
    s2 = Polysimplex({(b, a, c): -1, (d, e): 1})
    assert s1 == s2

def test_polysimplex_add():
    """Oppositely oriented copies cancel in addition."""
    s1 = Polysimplex({(a, b, c): 1, (d, e): 1})
    s2 = Polysimplex({(b, a, c): 1, (e, d): 1})
    assert (s1 + s2).is_zero()

def test_polysimplex_sub():
    """Subtracting an orientation-adjusted copy gives zero."""
    s1 = Polysimplex({(a, b, c): 1, (d, e): 1})
    s2 = Polysimplex({(b, a, c): -1, (e, d): -1})
    assert (s1 - s2).is_zero()

def test_polysimplex_mul():
    """The product of two disjoint atomic pairs is one four-point simplex."""
    s1 = Polysimplex({(a, b): 1})
    s2 = Polysimplex({(c, d): 1})
    prod = s1 * s2
    assert len(prod.terms) == 1
    assert list(prod.terms.keys())[0].components == (a, b, c, d)

def test_polysimplex_boundary():
    """Boundary is linear: the coefficient two multiplies every oriented face."""
    s = Polysimplex({Simplex([a, b ,c]): 2})
    bnd = s.boundary()
    assert bnd.terms[Simplex([b, c])] == 2
    assert bnd.terms[Simplex([a, c])] == -2
    assert bnd.terms[Simplex([a, b])] == 2

# ---------- Face ----------

def test_face_creation():
    """Two disjoint edge boundaries use four vertices and have total grade two."""
    f = Simplex([(a, b), (c, d)])
    assert f.size() == 4
    assert f.grade() == 2
    assert not f.is_zero()

def test_face_mul_single_overlap():
    """Boundaries with one common vertex merge: (ab)(bc)=(abc)."""
    f1 = Simplex([(a, b)])
    f2 = Simplex([(b, c)])
    prod = f1 * f2
    # (a,b)*(b,c) -> (a,b,c)
    assert len(prod.components) == 1
    assert tuple(prod.components[0]) == (a,b,c)

def test_face_mul_no_overlap():
    """Disjoint edge boundaries remain two factors of total grade two."""
    f1 = Simplex([(a, b)])
    f2 = Simplex([(c, d)])
    prod = f1 * f2
    assert prod.grade() == 2
    assert prod.components[0] == (a,b)
    assert prod.components[1] == (c,d)

def test_face_mul_simp():
    """Multiplication by an endpoint gives a(b-a)=ab because aa=0."""
    f1 = Simplex([a])
    f2 = Simplex([(a, b)])
    prod = f1 * f2
    assert prod == Simplex([a, b])

def test_face_mul_two_common():
    """Two common vertices force the boundary product to vanish."""
    f1 = Simplex([(a, b, c)])
    f2 = Simplex([(a, b)])
    prod = f1 * f2
    # Two shared vertices make the boundary product zero.
    assert prod.is_zero()
