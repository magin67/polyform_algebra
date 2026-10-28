"""Regression tests for equivalence, exponentials and basis link invalidation.

PMG definitions: https://github.com/magin67/pmg-notes/blob/v5/pmg-library/ru/Основания/Спецификация%20ядра%20PMG.md
Sections: formal forms, multiplication, nilpotency and the metric polyform.
Only supported quadratic forms are exposed by these tests.
"""
from math import exp

import numpy as np
import pytest

from polyform_algebra import Basis, Point, Vector, Simplex, Polysimplex, Polyform, quForm


@pytest.mark.parametrize("kind", ["unit", "point", "vector", "edge", "top"])
def test_nonzero_forms_are_not_equivalent_to_zero(kind):
    """Annihilation by every available edge does not imply a form is zero."""
    a, b, c = Point.create_list("abc")
    forms = {
        "unit": quForm.one(), "point": quForm([a]),
        "vector": quForm([Vector("v")]), "edge": quForm([(a, b)]),
        "top": quForm([(a, b, c)]),
    }
    form = forms[kind].as_polyform()
    assert not form.is_equivalent(Polyform.zero())
    assert not Polyform.zero().is_equivalent(form)
    assert form.is_equivalent(Polyform(form))
    assert not form.is_equivalent(2 * form)


def test_equivalence_respects_grades_and_scalar_component():
    """Neither a scalar nor a top-grade difference may disappear in a comparison."""
    a, b, c = Point.create_list("abc")
    edge = quForm([(a, b)]).as_polyform()
    top = quForm([(a, b, c)]).as_polyform()
    assert not (edge + top).is_equivalent(edge)
    assert not (edge + Polyform.one()).is_equivalent(edge)
    assert not top.is_equivalent(edge)
    assert not edge.is_equivalent(1)


def test_structural_identity_survives_multiplication():
    """The four-vertex relation and its product remain valid in quadratic coordinates."""
    from itertools import combinations
    a, b, c, d, e = Point.create_list("abcde")
    triangles = Polyform.from_pairs([([t], 1) for t in combinations([a,b,c,d], 3)])
    pairs = Polyform.from_pairs([([(a,b),(c,d)],1), ([(a,c),(b,d)],1), ([(a,d),(b,c)],1)])
    assert triangles != pairs  # Distinct stored sums can represent the same form.
    assert triangles.is_equivalent(pairs)
    assert (triangles * quForm([(d,e)])).is_equivalent(pairs * quForm([(d,e)]))
    assert not triangles.is_equivalent(pairs + quForm([(a,b,c)]))


def test_equivalence_uses_absolute_coefficient_tolerance():
    """The requested tolerance applies after expansion, including scalar coefficients."""
    x = 2 * Polyform.one()
    y = (2 + 1e-6) * Polyform.one()
    assert x.is_equivalent(y, tolerance=1e-5)
    assert not x.is_equivalent(y, tolerance=1e-7)
    assert x.is_equivalent(x, tolerance=0)
    with pytest.raises(ValueError, match="Tolerance"):
        x.is_equivalent(y, tolerance=-1)


@pytest.mark.parametrize("cls", [Polysimplex, Polyform])
@pytest.mark.parametrize("scalar", [0, 1, -2, 0.5])
def test_scalar_exponential(cls, scalar):
    """The zero exponential is the unit; scalar units use the ordinary real exponential."""
    result = (scalar * cls.one()).exp()
    assert len(result.terms) == 1
    assert next(iter(result.terms)).is_one()
    assert next(iter(result.terms.values())) == pytest.approx(exp(scalar))


def test_polysimplex_exponential_keeps_factorials_and_constant():
    """For N=ab+cd, N^2=2abcd and exp(N)=1+ab+cd+abcd."""
    a,b,c,d = Point.create_list("abcd")
    n = a*b + c*d
    expected = Polysimplex.one() + n + a*b*c*d
    assert n.exp() == expected
    assert (2*Polysimplex.one()+n).exp() == exp(2)*expected
    assert n.exp(max_terms=0).is_one()
    assert n.exp(max_terms=1) == Polysimplex.one()+n
    assert n.exp(max_terms=2) == expected


@pytest.mark.parametrize("limit", [-1, 1.5, True])
def test_exponential_rejects_invalid_truncation(limit):
    """A truncation limit counts non-negative integer powers of the nilpotent part."""
    with pytest.raises(ValueError, match="max_terms"):
        Polysimplex.zero().exp(max_terms=limit)


def test_polyform_exponential_mixed_grades_and_context():
    """Disjoint simple squares commute; their cross term has coefficient 2*3=6."""
    a,b,c,d,e,f = Point.create_list("abcdef")
    x = quForm([(a,b)])
    y = quForm([(c,d,e)])
    n = Polyform(2*x+3*y, vertices=[a,b,c,d,e,f])
    expected = Polyform.one()+2*x+3*y+6*(x*y)
    assert n.exp() == expected
    assert (n+Polyform.one()).exp() == exp(1)*expected
    assert (n+Polyform.one()).exp().vertices == n.vertices
    assert (n.exp()*(-n).exp()).is_equivalent(Polyform.one())


def linked_bases():
    """Create two genuinely linked point bases with a nontrivial transition matrix."""
    a,b = Point.create_list("ab")
    p = Point("p", frame=a.as_polysimplex())
    q = Point("q", frame=(a+b)/2)
    source, target = Basis([p,q]), Basis([a,b])
    source.link_to(target)
    return source, target


@pytest.mark.parametrize("side", [0,1])
@pytest.mark.parametrize("operation", ["add", "remove"])
def test_basis_mutation_invalidates_both_directions(side, operation):
    """Changing either endpoint must reject both forward and reverse stale transforms."""
    source, target = linked_bases()
    changed = (source,target)[side]
    if operation == "add":
        changed.add_element(Point("extra"))
    else:
        changed.remove_element(changed[0])
    assert target not in source.transitions
    assert source not in target.transitions
    for left,right in [(source,target),(target,source)]:
        with pytest.raises(ValueError):
            left[0].as_polysimplex().transform(left,right)
    assert all(changed.index(e) == i for i,e in enumerate(changed.elements))


def test_noop_mutations_keep_links_and_relink_restores_transforms():
    """No-op edits retain valid links; restoring a basis permits explicit relinking."""
    source,target = linked_bases()
    forward = source.transitions[target][0].copy()
    source.add_element(source[0])
    source.remove_element(Point("absent"))
    np.testing.assert_allclose(source.transitions[target][0],forward)
    extra = Point("extra")
    target.add_element(extra)
    target.remove_element(extra)
    source.link_to(target)
    p,q = source.elements
    a,b = target.elements
    assert q.as_polysimplex().transform(source,target) == (a+b)/2
    assert b.as_polysimplex().transform(target,source) == 2*q-p


def test_basis_invalidation_preserves_unrelated_links():
    """Editing A removes A-B and A-C but leaves a valid independent B-C link."""
    source,target = linked_bases()
    third = Basis(list(target.elements))
    source.link_to(third)
    target.link_to(third)
    source.add_element(Point("extra"))
    assert source.transitions == {}
    assert source not in target.transitions and source not in third.transitions
    assert third in target.transitions and target in third.transitions


def test_exterior_exponential_respects_odd_grade_signs():
    """Odd atoms anticommute, so exp(a+b)=1+a+b without an ab term."""
    a,b,c = Point.create_list("abc")
    assert (a+b).exp() == Polysimplex.one()+a+b
    # An even bc commutes with a, giving N^2/2=abc for N=a+bc.
    assert (a.as_polysimplex()+b*c).exp() == Polysimplex.one()+a+b*c+a*b*c
