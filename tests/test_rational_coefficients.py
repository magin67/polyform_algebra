"""T29: exact scalar operations reuse the existing exterior/quadratic algebra.

Source: PMG v5 Core Specification rev. 27, exterior product, boundary and
polyform exponential. See the source URL in test_core_data.py.
"""
from fractions import Fraction as F
from numbers import Rational

import pytest
import sympy as sp

from polyform_algebra import Point, Simplex, Polysimplex, Polyform, quForm


@pytest.mark.parametrize('q', [F(1, 10**12), sp.Rational(1, 10**12)])
def test_tiny_rationals_survive_all_scalar_paths(q):
    """A coefficient below the legacy float tolerance stays nonzero and exact."""
    a, b = Point.create_list('ab')
    for term in (a, Simplex([a, b]), quForm([(a, b)])):
        p = q * term
        assert p == term * q
        assert not p.is_zero()
        assert list(p.terms.values()) == [F(1, 10**12)]
        assert list((p / 3).terms.values()) == [F(1, 3 * 10**12)]
        assert list((-p).terms.values()) == [-F(1, 10**12)]
        assert (p-p).is_zero()
        assert list((p+p).terms.values()) == [F(2, 10**12)]
        assert p.extract_grade(p.grade()) == p


def test_exact_boundary_and_factorial_series():
    """For X=ab+cd, exp(X)=1+ab+cd+abcd; the 1/2! cancels exactly."""
    a, b, c, d = Point.create_list('abcd')
    x = F(1, 3)*(a*b) + F(1, 7)*(c*d)
    expected = Polysimplex.one() + x + F(1, 21)*Simplex([a, b, c, d])
    result = x.exp()
    assert result == expected
    assert all(isinstance(v, Rational) for v in result.terms.values())
    assert (F(1, 10**12)*Simplex([a,b,c])).boundary().boundary().is_zero()
    assert (a/3).terms == (F(1,3)*a).terms


def test_rational_form_exp_potential_and_polarization():
    """The weighted two-edge tree has tau=1/6 and endpoint resistance 5 exactly."""
    a, b, c = Point.create_list('abc')
    lap = Polyform({quForm([(a,b)]): F(1,2), quForm([(b,c)]): F(1,3)}, vertices=[a,b,c])
    metric = lap.exp()
    assert metric.potential(Polyform.one()) == F(1,6)
    assert metric.norm(quForm([(a,c)])) == F(5)
    assert all(isinstance(v, Rational) for v in metric.terms.values())
    polar = Polyform.polar_form(F(1,10**12)*a, b)
    assert not polar.is_zero()
    assert all(isinstance(v, Rational) for v in polar.terms.values())


def test_legacy_float_rounding_still_available():
    """Existing approximate callers keep their float pruning; rationals opt into exactness."""
    a = Point('a')
    assert (1e-12*a).is_zero()
    assert not (F(1,10**12)*a).is_zero()


def test_exact_equivalence_does_not_treat_a_small_form_as_zero():
    """The zero test must agree with preserved rational coefficients (T29)."""
    a, b = Point.create_list('ab')
    form = F(1, 10**12) * quForm([(a, b)])
    assert not form.is_equivalent(Polyform.zero())
    assert form.is_equivalent(Polyform(form))


def test_legacy_metric_outputs_remain_numpy_compatible():
    """Integer graph callers retain float norms; explicit fractions retain exact norms."""
    import numpy as np
    a, b = Point.create_list('ab')
    pair = quForm([(a, b)])
    numeric_metric = Polyform({pair: 3}, vertices=[a, b]).exp()
    value = numeric_metric.norm(pair)
    assert isinstance(value, float)
    assert np.linalg.det(np.array([[value]])) == pytest.approx(1/3)
    exact_metric = Polyform({pair: F(3)}, vertices=[a, b]).exp()
    assert isinstance(exact_metric.norm(pair), F)
