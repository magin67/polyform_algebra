"""T42: finite checks of notes, separate from the core API acceptance tests.

Pinned source paths, blobs and explicit numbers: reference_corpus_v1.json.
These checks do not add thematic operations or general bilinear forms to PMG.
"""
from fractions import Fraction as Q
from itertools import combinations
from math import factorial

import pytest
import sympy as sp

from polyform_algebra.core import (
    network, laplacian, resistance, spanning_tree_weight, Coordinates, gram,
    dirichlet, kron_reduce, Boundary, Unit, Form, AlgebraLimits, simple_term,
    square, form_exp, laplacian_form, form_product, grade_part, equivalent,
    potential, evaluate_form, spanning_coefficient,
)
from polyform_algebra.core.algebra import _pairs
from .test_contract_corpus import case, edges
from . import oracles


@pytest.mark.parametrize('name', ['weighted_triangle', 'weighted_tetrahedron'])
def test_note_exponential_minors_and_complementary_duality(name):
    """Each exponential grade equals forest incidence minors; duality is signed Jacobi inversion."""
    c = case(name)
    n, e = c['n'], edges(c)
    net = network(tuple(range(n)), e)
    metric = form_exp(laplacian_form(net).value).value
    l = oracles.laplacian(n, e)
    for k in range(n):
        actual = _pairs(grade_part(metric, k).value, AlgebraLimits())
        expected = oracles.forest_tensor(n, e, k)
        assert actual == expected
        minors = {(I, J): Q(l.extract(I, J).det())
                  for I in combinations(range(n), k) for J in combinations(range(n), k)}
        assert expected == {ij: c for ij, c in minors.items() if c}
    oracles.minor_duality(l[:-1, :-1])


def test_note_boundary_path_groves_and_gluing():
    """P5 has grove scale 4 and boundary Gram [[2,2],[2,4]]; two copies glue with tau=16."""
    c = case('boundary_path')
    n, e, p = c['n'], edges(c), tuple(c['selected'])
    groves = oracles.grove_weights(n, e, p)
    assert groves == {((0,), (2,), (4,)): 4, ((0, 2), (4,)): 2, ((0,), (2, 4)): 2, ((0, 2, 4),): 1}
    net = network(tuple(range(n)), e)
    reduction = kron_reduce(net, p).value
    expected_response = sp.Matrix([[sp.Rational(x) for x in row] for row in c['response']])
    assert sp.Matrix(laplacian(reduction.network).value.values) == expected_response
    vectors = [Coordinates(net.context, tuple(int(i == v)-int(i == p[0]) for i in range(n))) for v in p[1:]]
    g = sp.Matrix(gram(net, vectors).value.values)
    assert g == sp.Matrix([[sp.Rational(x) for x in row] for row in c['gram']])
    assert g.det() == Q(c['delta']) == groves[((0,), (2,), (4,))]
    assert spanning_tree_weight(reduction.network).value == Q(1, 4)
    # Distinct internal vertices in the second path, with the same three marked vertices.
    joined_edges = e + ((0, 5, 1), (5, 2, 1), (2, 6, 1), (6, 4, 1))
    tau, _ = oracles.tree_and_separating_forests(7, joined_edges)
    joined = network(tuple(range(7)), joined_edges)
    assert tau == spanning_tree_weight(joined).value == (g+g).det() == 16
    red = kron_reduce(joined, p).value.network
    assert sp.Matrix(laplacian(red).value.values) == 2*expected_response
    # Check an arbitrary non-unit basis normalization in the determinant gluing law.
    a = sp.Matrix([[2, 1], [0, 1]])
    changed = a.T*g*a
    assert (changed+changed).det()/a.det()**2 == tau


@pytest.mark.parametrize('t', [Q(-1), Q(-1, 2), Q(0), Q(1), Q(3, 2)])
def test_note_single_edge_variation(t):
    """One triangle conductance 1+t gives tau=3+2t and separating forest weights 2,2+t,2+t."""
    e = ((0, 1, 1+t), (0, 2, 1), (1, 2, 1))
    e = tuple(edge for edge in e if edge[2])
    tau, separating = oracles.tree_and_separating_forests(3, e)
    assert tau == 3+2*t
    assert (separating[0][1], separating[0][2], separating[1][2]) == (2, 2+t, 2+t)
    net = network((0, 1, 2), e)
    assert resistance(net).value.values == tuple(tuple(x/tau for x in row) for row in separating)
    changed = form_exp(laplacian_form(net).value).value
    # Reconstruct the original Laplacian in the SAME context to compare forms.
    base = Form(net.context, [(Boundary(net.context, ij), 1) for ij in combinations(range(3), 2)])
    m = form_exp(base).value
    edge_form = square(Boundary(net.context, (0, 1))).value
    expected = form_product(m, square(Unit(net.context)).value+t*edge_form).value
    assert equivalent(changed, expected).value
    assert spanning_coefficient(changed).value == tau


def test_note_higher_grades_from_complementary_forests():
    """A triangle boundary in unit C4 has norm 1/2, potential 2 and area squared 1/8."""
    e = ((0, 1, 1), (1, 2, 1), (2, 3, 1), (0, 3, 1))
    net = network((0, 1, 2, 3), e)
    # Each extra edge completes the two fixed directions to a tree or to zero.
    fixed = [sp.Matrix([-1, 1, 0]), sp.Matrix([-1, 0, 1])]
    expected = Q(0)
    for a, b, c in e:
        v = sp.Matrix([int(i == b)-int(i == a) for i in range(3)])
        expected += c*Q(sp.Matrix.hstack(*fixed, v).det())**2
    assert expected == 2
    f = square(Boundary(net.context, (0, 1, 2))).value
    m = form_exp(laplacian_form(net).value).value
    assert potential(m, f).value == expected
    assert evaluate_form(net, f).value == Q(1, 2)
    r = oracles.grounded_resistance(4, e)
    assert oracles.cayley_menger(r, (0, 1, 2)) == Q(1, 8)


def test_note_structural_rank_and_symmetric_complete_graph():
    """Tetrahedron forms have rank six; K5 forest counts differ from rewritten coefficient sums."""
    c = case('unit_k5')
    n, e = c['n'], edges(c)
    ordinary, _ = oracles.forest_statistics(n, e)
    assert ordinary == list(map(Q, c['ordinary'])) == [1, 10, 45, 110, 125]
    net = network(tuple(range(n)), e)
    m = form_exp(laplacian_form(net).value).value
    for k in range(1, n):
        symmetric = Form(net.context, [(Boundary(net.context, I), n**(k-1)) for I in combinations(range(n), k+1)])
        assert equivalent(grade_part(m, k).value, symmetric).value
    rewritten = grade_part(Form(net.context, [(Boundary(net.context, I), 5) for I in combinations(range(n), 3)]), 2).value
    assert sum(c for _, c in rewritten.terms) == 50 != ordinary[2]
    # Seven partition forms in grade two on four vertices: compute their rank
    # in an affine basis by determinants, independently of core normalization.
    basis = [sp.zeros(3, 1)] + [sp.eye(3)[:, i] for i in range(3)]
    factors = [[(I[0], I[1]), (I[0], I[2])] for I in combinations(range(4), 3)]
    factors += [[(0, 1), (2, 3)], [(0, 2), (1, 3)], [(0, 3), (1, 2)]]
    columns = []
    for pairs in factors:
        matrix = sp.Matrix.hstack(*(basis[b]-basis[a] for a, b in pairs))
        x = [matrix.extract(I, (0, 1)).det() for I in combinations(range(3), 2)]
        columns.append(sp.Matrix([a*b for a in x for b in x]))
    matrix = sp.Matrix.hstack(*columns)
    assert matrix.rank() == 6
    assert matrix*sp.Matrix([1, 1, 1, 1, -1, -1, -1]) == sp.zeros(9, 1)


@pytest.mark.parametrize('name', ['hexagons_A', 'hexagons_B'])
def test_note_three_hexagons_explicit_forests_and_spectrum(name):
    """All 65536 subsets per isomer reproduce the pinned ordinary/rooted counts and tau=204."""
    c = case(name)
    n, e = c['n'], edges(c)
    ordinary, rooted = oracles.forest_statistics(n, e)
    assert ordinary == list(map(Q, c['ordinary']))
    assert rooted == list(map(Q, c['rooted']))
    net = network(tuple(range(n)), e)
    l = sp.Matrix(laplacian(net).value.values)
    coefficients = l.charpoly().all_coeffs()
    assert coefficients == [(-1)**k*rooted[k] for k in range(n)] + [0]
    assert spanning_tree_weight(net).value == ordinary[-1] == Q(c['tau']) == 204
    r = resistance(net).value.values
    assert sum(r[i][j] for i, j in combinations(range(n), 2)) == Q(c['kirchhoff'])
    assert [sp.trace(l**k) for k in range(1, 5)] == list(map(Q, c['traces']))
    grounded = oracles.laplacian(n, e)[:-1, :-1]
    inverse = grounded.inv()
    diagonal = sp.Matrix([inverse[i, i] for i in range(n-1)])
    assert (diagonal.T*grounded*diagonal)[0]/4 == Q(c['radius'])


def test_note_isomer_polynomial_difference():
    """The first rooted difference is at grade four; full difference equals the squared C5 polynomial."""
    a, b = case('hexagons_A'), case('hexagons_B')
    assert a['ordinary'] == b['ordinary']
    assert a['rooted'][:4] == b['rooted'][:4]
    assert int(a['rooted'][4])-int(b['rooted'][4]) == 1
    x = sp.Symbol('x')
    polynomial = sum((-1)**k*(int(u)-int(v))*x**(14-k) for k, (u, v) in enumerate(zip(a['rooted'], b['rooted'])))
    assert sp.expand(polynomial-x**2*(x*x-5*x+5)**4) == 0
