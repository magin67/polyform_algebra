"""Formal contracts A01-A12, independent coordinate and forest references.

PMG v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0:
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md
Sections: exterior objects, forms, fixed spanning-tree coefficient, metric identity.
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Polyform%20Algebra%20and%20Combinatorics/Structural%20identities/Partition%20Forms%20and%20the%20Tetrahedron%20Identity.md
Equation (1) and its products; the form identity is independent of any metric.
"""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
from itertools import combinations, permutations
from math import factorial

import pytest
import sympy as sp

from polyform_algebra.core import (
    Context, Coordinates, AlgebraLimits, Exterior, SimpleTerm, Form,
    Atom, Boundary, Unit, exterior, simple_term, form,
    wedge, boundary, grade_part, square, form_product, equivalent, polar_form,
    form_exp, laplacian_form, spanning_coefficient, potential, evaluate_form,
    network, spanning_tree_weight, resistance, inner_product,
    InvalidInput, DomainError, UnsupportedOperation, ResourceLimit, ContextMismatch,
)


def sq(x):
    return square(x).value


def mul(x, y):
    return form_product(x, y).value


def eq(x, y):
    return equivalent(x, y).value


def test_exterior_normalization_and_context():
    """T19/T29: permutation signs accumulate; only exact zeros disappear."""
    c = Context(('a', 'b', 'c'))
    x = exterior(c, [((1, 0), F(1, 3)), ((0, 1), F(1, 3)), ((0, 0), 9), ((2,), F(1, 10**12))])
    assert x.coefficients == (((2,), F(1, 10**12)),)
    assert (x / 3 * 3 - x).coefficients == ()
    assert (x-x).context is c
    with pytest.raises(FrozenInstanceError):
        x.context = Context(('a', 'b', 'c'))
    assert grade_part(x, 0).value == Exterior(c, ())
    assert grade_part(x, 1).value == x


def test_wedge_algebra_laws():
    """T19: graded signs, associativity, bilinearity, identity and repeated atoms."""
    c = Context(tuple('abcd'))
    x = exterior(c, {(0,): 2, (1,): F(1, 3)})
    y = exterior(c, {(1,): 3, (2,): -2})
    z = Atom(c, 3).as_exterior()
    xy = wedge(x, y).value
    assert xy == -wedge(y, x).value
    assert wedge(xy, z).value == wedge(x, wedge(y, z).value).value
    assert wedge(x+y, z).value == wedge(x, z).value + wedge(y, z).value
    assert wedge(x, x).value.coefficients == ()
    assert wedge(x, Unit(c)).value == x
    assert wedge(x, exterior(c, {})).value == exterior(c, {})
    assert wedge(xy, z).value == wedge(z, xy).value  # (-1)^(2*1)=+1.


@pytest.mark.parametrize('epsilon', [(1, 1, 1, 1), (1, 2, F(1, 3), -2), (0, 0, 0, 0)])
def test_boundary_squared_and_leibniz(epsilon):
    """T20: contraction squares to zero and obeys the graded Leibniz identity."""
    c = Context(tuple('abcd'), epsilon=epsilon)
    x = exterior(c, {(0, 1): 2, (1, 2): F(1, 3)})
    y = exterior(c, {(2,): 3, (3,): -2})
    assert boundary(boundary(x).value).value == exterior(c, {})
    lhs = boundary(wedge(x, y).value).value
    rhs = wedge(boundary(x).value, y).value + wedge(x, boundary(y).value).value
    assert lhs == rhs  # x has even grade two.
    lhs = boundary(wedge(y, x).value).value
    rhs = wedge(boundary(y).value, x).value - wedge(y, boundary(x).value).value
    assert lhs == rhs  # y has odd grade one.
    assert boundary(Unit(c)).value == exterior(c, {})


def test_weighted_boundary_regression():
    """T30: epsilon(b)=2 gives partial(ab)*partial(bc)=2*partial(abc)."""
    c = Context(tuple('abc'), epsilon=(1, 2, 1))
    ab = boundary(exterior(c, {(0, 1): 1})).value
    bc = boundary(exterior(c, {(1, 2): 1})).value
    abc = boundary(exterior(c, {(0, 1, 2): 1})).value
    assert wedge(ab, bc).value == 2*abc
    with pytest.raises(UnsupportedOperation):
        Boundary(c, (0, 1))


def determinant_coordinates(columns, n):
    """Independent expansion of a decomposable object by minors of its column matrix."""
    if len(columns) > n:
        return {}
    if not columns:
        return {(): F(1)}
    matrix = sp.Matrix.hstack(*columns)
    result = {}
    for rows in combinations(range(n), len(columns)):
        value = F(matrix.extract(rows, range(len(columns))).det())
        if value:
            result[rows] = value
    return result


def test_1600_oriented_simple_products_against_minors():
    """T21/T41: 40 oriented atoms/boundaries give 1600 products checked by determinants."""
    c = Context(tuple('abcd'))
    basis = [sp.eye(4)[:, i] for i in range(4)]
    cases = [(Atom(c, i), [basis[i]]) for i in range(4)]
    for size in (2, 3):
        for indices in permutations(range(4), size):
            cases.append((Boundary(c, indices), [basis[i]-basis[indices[0]] for i in indices[1:]]))
    assert len(cases) == 40
    for x, a in cases:
        for y, b in cases:
            expected = determinant_coordinates(a+b, 4)
            actual = wedge(x, y).value
            compact = simple_term(c, (x, y)).as_exterior()
            assert dict(actual.coefficients) == expected
            assert compact == actual


def test_simple_sign_square_and_commuting_forms():
    """T22: exterior orientation survives products but cancels in quadratic squares."""
    c = Context(tuple('abcd'))
    a, b = Boundary(c, (1, 0)), Boundary(c, (2, 3))
    assert a.as_exterior() == -Boundary(c, (0, 1)).as_exterior()
    assert simple_term(c, (a, b)).as_exterior() == -simple_term(c, (b, a)).as_exterior()
    assert eq(sq(a), sq(-a))
    assert eq(mul(sq(a), sq(b)), mul(sq(b), sq(a)))
    assert eq(mul(sq(a), sq(b)), sq(simple_term(c, (a, b))))
    assert sq(0*a) == form(c, {})
    assert mul(sq(a), sq(a)) == form(c, {})
    assert eq(sq(Unit(c)), sq(Boundary(c, (0,))))
    assert not eq(sq(Unit(c)), form(c, {}))
    assert not eq(sq(Boundary(c, (0, 1, 2, 3))), form(c, {}))
    assert eq(sq(F(1, 3)*a), sq(a)/9)
    assert grade_part(sq(Boundary(c, (0, 1, 2))), 2).value == sq(Boundary(c, (0, 1, 2)))
    assert grade_part(sq(Boundary(c, (0, 1, 2))), 3).value == form(c, {})


def test_tetrahedron_identity_and_product():
    """T23: sum of four triangle squares equals three opposite-pair products, before metrics."""
    c = Context(tuple('abcde'))
    left = form(c, [(Boundary(c, triple), 1) for triple in combinations(range(4), 3)])
    pairs = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
    right = form(c, [(simple_term(c, (Boundary(c, a), Boundary(c, b))), 1) for a, b in pairs])
    assert eq(left, right)
    assert not eq(left, right + sq(Unit(c)))
    factor = sq(Boundary(c, (3, 4)))
    assert eq(mul(left, factor), mul(right, factor))
    # Direct ordered coordinate-pair coefficients are an independent reference
    # for the private comparison representation, not a new public constructor.
    from polyform_algebra.core.algebra import _pairs
    basis = [sp.eye(5)[:, i] for i in range(5)]
    reference = {}
    for triple in combinations(range(4), 3):
        columns = [basis[i]-basis[triple[0]] for i in triple[1:]]
        coefficients = determinant_coordinates(columns, 5)
        for a, x in coefficients.items():
            for b, y in coefficients.items():
                reference[a, b] = reference.get((a, b), 0) + x*y
    assert _pairs(right, AlgebraLimits()) == {key: value for key, value in reference.items() if value}


def test_polarization_normalization_and_cancellation():
    """T24: the polar form is u*v^T+v*u^T, with no implicit one-half."""
    c = Context(tuple('abc'))
    u = exterior(c, {(0,): 2, (1,): -1, (2,): F(1, 3)})
    v = exterior(c, {(0,): -2, (2,): 3})
    from polyform_algebra.core.algebra import _pairs
    pairs = _pairs(polar_form(u, v).value, AlgebraLimits())
    a, b = dict(u.coefficients), dict(v.coefficients)
    expected = {((i,), (j,)): a.get((i,), 0)*b.get((j,), 0)+b.get((i,), 0)*a.get((j,), 0)
                for i in range(3) for j in range(3)}
    assert pairs == {key: value for key, value in expected.items() if value}
    ab = Boundary(c, (0, 1))
    assert eq(polar_form(ab, ab).value, 2*sq(ab))
    assert polar_form(exterior(c, {}), u).value == form(c, {})
    assert sq(Atom(c, 0)).terms[0][0].grade == 1


def test_boundary_membership_after_cancellation():
    """A10/A12 must inspect the whole form: non-boundary atom terms can cancel exactly."""
    net = network(tuple('ab'), (('a', 'b', 2),))
    c = net.context
    a, b = Atom(c, 0), Atom(c, 1)
    rewritten = sq(a) + sq(b) - polar_form(a, b).value
    assert eq(rewritten, sq(Boundary(c, (0, 1))))
    assert spanning_coefficient(rewritten).value == 1
    assert evaluate_form(net, rewritten).value == F(1, 2)
    with pytest.raises(DomainError):
        evaluate_form(net, sq(a))
    with pytest.raises(DomainError):
        spanning_coefficient(sq(a))


def test_boundary_membership_with_nonboundary_terms_remaining():
    """T23/T27: an atomic multiple of a zero identity must cancel before testing closure."""
    net = network(tuple(range(5)), ((0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1)))
    c = net.context
    triangles = form(c, [(Boundary(c, t), 1) for t in combinations(range(4), 3)])
    opposite = form(c, [(simple_term(c, (Boundary(c, a), Boundary(c, b))), 1)
        for a, b in (((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2)))])
    null = mul(triangles-opposite, sq(Atom(c, 4)))
    assert null.terms and all(any(kind == 'atom' for kind, _ in t.factors) for t, _ in null.terms)
    f = sq(Boundary(c, (0, 1)))
    assert eq(null+f, f)
    assert spanning_coefficient(null+f).value == 0
    assert evaluate_form(net, null+f).value == evaluate_form(net, f).value == 1


def test_exponential_factorials_and_mixed_grades():
    """T25: a mixed-grade polynomial has exp(P)=e+P+P^2/2+..., not termwise evaluation."""
    c = Context(tuple('abcd'))
    unit = sq(Unit(c))
    p = 2*sq(Atom(c, 0)) + 3*sq(Atom(c, 1)) + 5*sq(simple_term(c, (Atom(c, 2), Atom(c, 3))))
    expected, power = unit, unit
    for k in range(1, 5):
        power = mul(power, p)
        expected = expected + power / factorial(k)
    actual = form_exp(p).value
    assert eq(actual, expected)
    assert eq(mul(actual, form_exp(-p).value), unit)
    assert eq(grade_part(actual, 4).value, 30*sq(simple_term(c, tuple(Atom(c, i) for i in range(4)))))
    zero = form(c, {})
    assert eq(form_exp(zero).value, unit)
    q = F(1, 10**12)*sq(Boundary(c, (0, 1)))
    assert eq(form_exp(q).value, unit + q)
    assert eq(grade_part(form_exp(q).value, 1).value, q)
    with pytest.raises(UnsupportedOperation):
        form_exp(p + unit)


def enumerate_trees(n, edges):
    """Enumerate edge subsets and connectivity, without cofactors or formal multiplication."""
    if n == 1:
        return F(1)
    total = F(0)
    for subset in combinations(edges, n-1):
        neighbors = [set() for _ in range(n)]
        weight = F(1)
        for a, b, c in subset:
            neighbors[a].add(b)
            neighbors[b].add(a)
            weight *= c
        visited, pending = {0}, [0]
        while pending:
            current = pending.pop()
            for nxt in neighbors[current] - visited:
                visited.add(nxt)
                pending.append(nxt)
        if len(visited) == n:
            total += weight
    return total


def test_formal_tree_coefficients_all_graphs_up_to_four_vertices():
    """T26: every simple graph on 1..4 vertices matches independent tree enumeration."""
    for n in range(1, 5):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            edges = tuple((a, b, 1) for i, (a, b) in enumerate(possible) if mask & (1 << i))
            net = network(tuple(range(n)), edges)
            m = form_exp(laplacian_form(net).value).value
            expected = enumerate_trees(n, edges)
            assert spanning_coefficient(m).value == expected
            assert spanning_tree_weight(net).value == expected


@pytest.mark.parametrize('edges', [((0, 1, 2), (1, 2, 3), (0, 2, 4)),
    ((0, 1, F(1, 3)), (1, 2, F(1, 10**12))), ((0, 1, 2),)])
def test_weighted_coefficients_and_isolates(edges):
    """T26/T29: weighted trees retain rational factors; an isolate fixes tau at zero."""
    net = network((0, 1, 2), edges)
    m = form_exp(laplacian_form(net).value).value
    assert spanning_coefficient(m).value == enumerate_trees(3, edges)
    assert m.context is net.context
    unit = sq(Unit(net.context))
    assert potential(m, unit).value == spanning_coefficient(m).value


@pytest.mark.parametrize('edges', [((0, 1, 2), (1, 2, 3), (0, 2, 4)), ((0, 1, 1), (1, 2, 1))])
def test_potential_metric_identity_and_nonmultiplicative_evaluation(edges):
    """T27: normalized potentials equal Gram evaluation, which is not multiplicative."""
    net = network((0, 1, 2), edges)
    c = net.context
    l = laplacian_form(net).value
    m = form_exp(l).value
    tau = spanning_coefficient(m).value
    p, q = sq(Boundary(c, (0, 1))), sq(Boundary(c, (0, 2)))
    product = mul(p, q)
    for f in (sq(Unit(c)), p, q, product, p-3*q+5*product):
        assert evaluate_form(net, f).value == potential(m, f).value / tau
    assert evaluate_form(net, p).value == resistance(net).value.values[0][1]
    assert potential(m, p+q).value == potential(m, p).value+potential(m, q).value
    u, v = Coordinates(c, (-1, 1, 0)), Coordinates(c, (-1, 0, 1))
    uv = inner_product(net, u, v).value
    assert evaluate_form(net, product).value == evaluate_form(net, p).value*evaluate_form(net, q).value-uv**2
    assert evaluate_form(net, product).value != evaluate_form(net, p).value*evaluate_form(net, q).value
    polar = polar_form(Boundary(c, (0, 1)), Boundary(c, (0, 2))).value
    assert evaluate_form(net, polar/2).value == uv


def test_zero_tau_potential_and_metric_domain():
    """Potential stays defined on disconnected graphs; normalization must be rejected."""
    net = network((0, 1, 2), ((0, 1, 2),))
    m = form_exp(laplacian_form(net).value).value
    f = sq(Boundary(net.context, (1, 2)))
    assert spanning_coefficient(m).value == 0
    assert potential(m, f).value == 2
    with pytest.raises(DomainError):
        evaluate_form(net, f)
    one = network(('a',), ())
    unit = sq(Unit(one.context))
    assert spanning_coefficient(unit).value == evaluate_form(one, unit).value == 1
    empty = network((), ())
    e = form_exp(laplacian_form(empty).value).value
    assert eq(e, sq(Unit(empty.context)))
    for operation in (lambda: spanning_coefficient(e), lambda: evaluate_form(empty, e)):
        with pytest.raises(DomainError):
            operation()


def test_numeric_formal_objects_are_rejected():
    """T36: the new formal layer is exact-only; legacy approximate classes remain separate."""
    c = Context(('a',), arithmetic='numeric')
    for operation in (lambda: exterior(c, {}), lambda: simple_term(c, ()),
                      lambda: form(c, {}), lambda: laplacian_form(network(('a',), (), arithmetic='numeric'))):
        with pytest.raises(UnsupportedOperation):
            operation()


def test_unsupported_square_and_polarization():
    """T36: neither arbitrary exterior sums nor higher-grade polar forms extend the grammar."""
    c = Context(tuple('abc'))
    for argument in (exterior(c, {(0,): 1, (1,): 1}), exterior(c, {}), sq(Unit(c))):
        with pytest.raises(UnsupportedOperation):
            square(argument)
    with pytest.raises(UnsupportedOperation):
        form(c, [(exterior(c, {(0,): 1}), 1)])
    with pytest.raises(UnsupportedOperation):
        polar_form(Boundary(c, (0, 1, 2)), Atom(c, 0))
    weighted = Context(tuple('ab'), epsilon=(1, 2))
    with pytest.raises(UnsupportedOperation):
        polar_form(Atom(weighted, 0), Atom(weighted, 1))
    with pytest.raises(UnsupportedOperation):
        exterior(c, {(0,): sp.Symbol('x')})


@pytest.mark.parametrize('indices', [(True,), (-1,), (3,), (0.0,)])
def test_invalid_atom_positions(indices):
    """Indices are explicit zero-based integer positions; bool and labels are not positions."""
    with pytest.raises(InvalidInput):
        exterior(Context(tuple('abc')), {indices: 1})


@pytest.mark.parametrize('coefficient', [True, 0.0, float('nan'), float('inf')])
def test_invalid_exact_coefficients(coefficient):
    """Exact algebra rejects floating and bool coefficients rather than guessing rationals."""
    c = Context(('a',))
    for operation in (lambda: exterior(c, {(0,): coefficient}),
                      lambda: form(c, [(Atom(c, 0), coefficient)])):
        with pytest.raises(InvalidInput):
            operation()


def test_formal_context_mismatch():
    """Matching labels do not identify independent spaces, including their zero elements."""
    a, b = Context(('a', 'b')), Context(('a', 'b'))
    for operation in (lambda: wedge(Atom(a, 0), Atom(b, 0)),
                      lambda: equivalent(form(a, {}), form(b, {})),
                      lambda: form_product(form(a, {}), form(b, {})),
                      lambda: simple_term(a, (Atom(b, 0),)),
                      lambda: form(a, [(Atom(b, 0), 1)]),
                      lambda: exterior(a, {}) + exterior(b, {})):
        with pytest.raises(ContextMismatch):
            operation()


def test_limits_cover_intermediates_and_keep_inputs():
    """T37: low budgets reject intermediate pair/product expansions without a partial result."""
    c = Context(tuple('abcd'))
    x = Boundary(c, (0, 1, 2))
    f = sq(x)
    p = sq(Atom(c, 0)) + sq(Atom(c, 1)) + sq(Atom(c, 2))
    original = p.terms
    for operation in (
        lambda: x.as_exterior(limits=AlgebraLimits(2)),
        lambda: equivalent(f, f, limits=AlgebraLimits(8)),
        lambda: form_product(p, p, limits=AlgebraLimits(8)),
        lambda: form_exp(p, limits=AlgebraLimits(3)),
        lambda: exterior(c, {(): 1, (0,): 1}, limits=AlgebraLimits(1)),
    ):
        with pytest.raises(ResourceLimit) as error:
            operation()
        assert error.value.code == 'RESOURCE_LIMIT'
    assert p.terms == original
    low = exterior(c, {(0,): 1, (1,): 1}, limits=AlgebraLimits(2))
    with pytest.raises(ResourceLimit):
        wedge(low, low)  # Four candidate products, even though the result cancels.
    assert wedge(low, low, limits=AlgebraLimits(4)).value.coefficients == ()


def test_low_budget_cannot_hide_cancelling_constructor_terms():
    """Intermediate storage is bounded even when later input entries would cancel it."""
    c = Context(tuple('ab'))
    with pytest.raises(ResourceLimit):
        exterior(c, [((0,), 1), ((1,), 1), ((1,), -1)], limits=AlgebraLimits(1))
    for bad in (0, -1, True, 1.5):
        with pytest.raises(InvalidInput):
            AlgebraLimits(bad)


def test_legacy_formal_api_agrees_on_triangle():
    """Compare the established Polyform API with the new exact context on a weighted triangle."""
    from polyform_algebra import Point, Polyform, quForm
    a, b, c = Point.create_list('abc')
    old = Polyform({quForm([(a, b)]): F(2), quForm([(b, c)]): F(3), quForm([(a, c)]): F(4)})
    net = network(tuple('abc'), (('a', 'b', 2), ('b', 'c', 3), ('a', 'c', 4)))
    m = form_exp(laplacian_form(net).value).value
    f = sq(Boundary(net.context, (0, 1)))
    old_m = old.exp()
    assert old_m.potential(Polyform.one()) == spanning_coefficient(m).value
    assert old_m.norm(Polyform({quForm([(a, b)]): F(1)})) == evaluate_form(net, f).value
