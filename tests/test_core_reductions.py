"""Stage-four contracts against PMG examples and independent identities.

Sources (pinned v5 revision, sections: path, star, matrix notation, limitations):
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Selected%20Vertices%20and%20the%20Dirichlet%20Problem.md
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Vertex%20Contraction%20and%20Graph%20Gluing.md
"""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F

import numpy as np
import pytest

from polyform_algebra.core import (
    Context, Coordinates, InvalidInput, DomainError, NonuniqueProblem,
    ContextMismatch, ArithmeticMismatch, NumericalUnresolved, INFINITY,
    network, laplacian, resistance, dirichlet, extend, kron_reduce, contract,
    induced_subgraph, reorder, to_numeric,
)


@pytest.fixture(params=['exact', 'numeric'])
def mode(request):
    return request.param


def path(mode):
    return network(('a', 'x', 'y', 'b'), (('a', 'x', 2), ('x', 'y', 1), ('y', 'b', 3)), arithmetic=mode)


def equal(actual, expected, mode):
    if mode == 'exact':
        np.testing.assert_array_equal(actual, expected)
    else:
        np.testing.assert_allclose(np.asarray(actual, dtype=float), np.asarray(expected, dtype=float), atol=1e-12, rtol=1e-10)


def test_weighted_path(mode):
    """T05: solve [[3,-1],[-1,4]] B=diag(2,3); E=2*9+36+3*4=66."""
    net = path(mode)
    calculation = dirichlet(net, ('a', 'b'))
    d = calculation.value
    assert d.P == ('a', 'b') and d.I == ('x', 'y')
    equal(d.B.to_numpy(), [[F(8, 11), F(3, 11)], [F(2, 11), F(9, 11)]], mode)
    assert d.H.rows is net.context and d.B.columns is d.H.columns
    f = extend(d, {'b': 11, 'a': 0}).value
    assert f.context is net.context
    equal(f.values, [0, 3, 9, 11], mode)
    assert d.extend(Coordinates(d.H.columns, (0, 11))).value == f
    l, v = laplacian(net).value.to_numpy(), np.array(f.values)
    equal((l @ v)[[1, 2]], [0, 0], mode)
    equal(v @ l @ v, 66, mode)
    h = np.array([0, 2, -1, 0])
    equal((v+h) @ l @ (v+h), 66+h @ l @ h, mode)
    assert h @ l @ h > 0
    assert calculation.diagnostics.status == ('exact' if mode == 'exact' else 'approximate')
    if mode == 'numeric':
        assert calculation.diagnostics.numerical[0].rank == 2
        assert calculation.diagnostics.numerical[0].residual < 1e-10


def test_three_reductions_differ(mode):
    """T06: deletion gives infinity, contraction 1, elimination series R=11/6."""
    net = path(mode)
    k = kron_reduce(net, ('a', 'b')).value
    equal(k.network.edges[0][2], F(6, 11), mode)
    equal(resistance(k.network).value.values[0][1], F(11, 6), mode)
    equal(k.H.to_numpy().T @ laplacian(net).value.to_numpy() @ k.H.to_numpy(), laplacian(k.network).value.to_numpy(), mode)
    assert k.network.context is k.H.columns
    deleted = induced_subgraph(net, ('a', 'b')).value
    assert deleted.edges == ()
    assert resistance(deleted).value.values[0][1] == (INFINITY if mode == 'exact' else float('inf'))
    glued = contract(net, {'A': ('a', 'x'), 'B': ('y', 'b')}).value
    equal(resistance(glued.network).value.values[0][1], 1, mode)
    s = glued.S.to_numpy()
    equal(s.T @ laplacian(net).value.to_numpy() @ s, laplacian(glued.network).value.to_numpy(), mode)


def test_star_and_sequential_elimination(mode):
    """T07/T17: a unit star induces c=1/3; sequential Schur complements agree."""
    star = network(('x', 'a', 'b', 'c'), (('x', 'a', 1), ('x', 'b', 1), ('x', 'c', 1)), arithmetic=mode)
    k = kron_reduce(star, ('c', 'a', 'b')).value
    equal([c for _, _, c in k.network.edges], [F(1, 3)] * 3, mode)
    equal(resistance(k.network).value.to_numpy(), [[0, 2, 2], [2, 0, 2], [2, 2, 0]], mode)
    net = path(mode)
    for middle in (('a', 'y', 'b'), ('a', 'x', 'b')):
        first = kron_reduce(net, middle).value
        second = kron_reduce(first.network, ('b', 'a')).value
        direct = kron_reduce(net, ('b', 'a')).value
        equal(laplacian(second.network).value.to_numpy(), laplacian(direct.network).value.to_numpy(), mode)
        equal(first.H.to_numpy() @ second.H.to_numpy(), direct.H.to_numpy(), mode)


def test_general_network_identities(mode):
    """T16/T17: nonsymmetric topology still gives symmetric response and preserved R."""
    net = network(tuple(range(6)), ((0, 1, 2), (1, 2, 3), (2, 3, 5), (3, 4, 7),
                                   (4, 5, 11), (5, 0, 13), (1, 4, 17)), arithmetic=mode)
    p = (5, 1, 3)
    d = dirichlet(net, p).value
    k = kron_reduce(net, p).value
    l, h, b = laplacian(net).value.to_numpy(), d.H.to_numpy(), d.B.to_numpy()
    assert np.all(b >= 0)
    equal(b.sum(axis=1), [1] * 3, mode)
    equal((l @ h)[[0, 2, 4]], np.zeros((3, 3), dtype=int), mode)
    equal(h.T @ l @ h, laplacian(k.network).value.to_numpy(), mode)
    equal(resistance(k.network).value.to_numpy(), resistance(net).value.to_numpy()[np.ix_(p, p)], mode)


def test_permutations_and_immutability(mode):
    """T12: source order affects rows, selected order columns; exports are copies."""
    net = path(mode)
    shuffled = reorder(net, ('y', 'b', 'x', 'a')).value
    d = dirichlet(shuffled, ('b', 'a')).value
    assert d.I == ('y', 'x')
    equal(extend(d, {'a': 0, 'b': 11}).value.values, [9, 11, 3, 0], mode)
    exported = d.B.to_numpy()
    exported[:] = 100
    assert d.B.values[0][0] != 100
    with pytest.raises(FrozenInstanceError):
        d.P = ()
    assert net.context.labels == ('a', 'x', 'y', 'b')
    numeric_h = to_numeric(d.H).value
    assert numeric_h.rows.identity == d.H.rows.identity
    equal(numeric_h.to_numpy(), d.H.to_numpy(), 'numeric')


def test_all_selected_singleton_and_empty(mode):
    """T09: no unknowns gives a 0-by-n B; a singleton has zero response."""
    net = path(mode)
    d = dirichlet(net, ('b', 'y', 'x', 'a')).value
    assert d.B.to_numpy().shape == (0, 4)
    equal(d.H.to_numpy(), np.eye(4, dtype=int)[::-1], mode)
    equal(laplacian(kron_reduce(net, d.P).value.network).value.to_numpy(), laplacian(net).value.to_numpy()[::-1, ::-1], mode)
    one = network(('z',), (), arithmetic=mode)
    assert dirichlet(one, ('z',)).value.B.to_numpy().shape == (0, 1)
    assert kron_reduce(one, ('z',)).value.network.edges == ()
    assert induced_subgraph(net, ()).value.context.labels == ()
    empty = network((), (), arithmetic=mode)
    assert contract(empty, {}).value.S.to_numpy().shape == (0, 0)
    for operation in (dirichlet, kron_reduce):
        with pytest.raises(DomainError):
            operation(empty, ())


def test_disconnected_anchors_and_single_selected(mode):
    """Each component, including an isolate, needs an anchor; cross influence is zero."""
    net = network(('a', 'x', 'b', 'y', 'z'), (('a', 'x', 2), ('b', 'y', 3)), arithmetic=mode)
    d = dirichlet(net, ('z', 'b', 'a')).value
    equal(d.B.to_numpy(), [[0, 0, 1], [0, 1, 0]], mode)
    assert kron_reduce(net, d.P).value.network.edges == ()
    with pytest.raises(NonuniqueProblem) as error:
        dirichlet(net, ('a', 'b'))
    assert error.value.code == 'NONUNIQUE_PROBLEM' and error.value.parameters == (('z',),)
    d = dirichlet(path(mode), ('a',)).value
    equal(d.B.to_numpy(), [[1], [1], [1]], mode)
    assert kron_reduce(path(mode), ('a',)).value.network.edges == ()


def test_contraction_and_induced_degrees(mode):
    """T18: disconnected members may merge; deletion recomputes degrees from edges."""
    net = path(mode)
    first = contract(net, {'A': ('a', 'b'), 'X': ('x',), 'Y': ('y',)}).value
    second = contract(first.network, {'U': ('A', 'X'), 'V': ('Y',)}).value
    direct = contract(net, {'U': ('a', 'b', 'x'), 'V': ('y',)}).value
    equal(laplacian(second.network).value.to_numpy(), laplacian(direct.network).value.to_numpy(), mode)
    equal(first.S.to_numpy() @ second.S.to_numpy(), direct.S.to_numpy(), mode)
    sub = induced_subgraph(net, ('y', 'x')).value
    equal(laplacian(sub).value.to_numpy(), [[1, -1], [-1, 1]], mode)
    assert sub.context.identity != net.context.identity


@pytest.mark.parametrize('selection', [(), ('a', 'a'), ('unknown',), {'a', 'b'}, (True,)])
def test_invalid_boundary(selection):
    """T35: ordered nonempty unique network labels are an API requirement."""
    for operation in (dirichlet, kron_reduce):
        with pytest.raises(InvalidInput):
            operation(path('exact'), selection)


@pytest.mark.parametrize('partition', [{}, {'u': ()}, {'u': ('a', 'x')},
    {'u': ('a', 'x'), 'v': ('x', 'y', 'b')}, {'u': ('a', 'x', 'y', 'b', 'z')},
    {'u': {'a', 'x', 'y', 'b'}}, {True: ('a', 'x', 'y', 'b')}, []])
def test_invalid_partition(partition):
    """T35: contraction requires a full disjoint partition with legal labels."""
    with pytest.raises(InvalidInput):
        contract(path('exact'), partition)


def test_extension_context_and_values():
    """Boundary values cannot silently cross spaces, arithmetic modes or label sets."""
    d = dirichlet(path('exact'), ('a', 'b')).value
    for bad in ({'a': 0}, {'a': 0, 'b': 1, 'x': 2}, {'a': 0.0, 'b': 1}, [0, 1]):
        with pytest.raises(InvalidInput):
            extend(d, bad)
    with pytest.raises(ContextMismatch):
        extend(d, Coordinates(Context(('a', 'b')), (0, 1)))
    with pytest.raises(ArithmeticMismatch):
        extend(d, to_numeric(Coordinates(d.H.columns, (0, 1))).value)


@pytest.mark.parametrize('scale', [1e-200, 1.0, 1e200])
def test_numeric_scaling(scale):
    """Uniform conductance scaling leaves influence unchanged and scales Kron weights."""
    net = network(('a', 'x', 'y', 'b'), (('a', 'x', 2*scale), ('x', 'y', scale), ('y', 'b', 3*scale)), arithmetic='numeric')
    equal(dirichlet(net, ('a', 'b')).value.B.to_numpy(), [[F(8, 11), F(3, 11)], [F(2, 11), F(9, 11)]], 'numeric')
    assert kron_reduce(net, ('a', 'b')).value.network.edges[0][2] / scale == pytest.approx(6/11)


def test_numeric_independent_component_scales():
    """Componentwise grounding prevents a strong component from erasing a weak one."""
    net = network(('a', 'x', 'b', 'y'), (('a', 'x', 1e200), ('b', 'y', 1e-200)), arithmetic='numeric')
    equal(dirichlet(net, ('a', 'b')).value.B.to_numpy(), [[1, 0], [0, 1]], 'numeric')


def test_numeric_huge_internal_contraction():
    """T31: sum inter-group edges directly, preserving 1 beside internal weights 1e100."""
    net = network(('a', 'x', 'y', 'b'), (('a', 'x', 1e100), ('x', 'y', 1), ('y', 'b', 1e100)), arithmetic='numeric')
    assert contract(net, {'u': ('a', 'x'), 'v': ('y', 'b')}).value.network.edges == (('u', 'v', 1.0),)


def test_numeric_ill_conditioned_grounding():
    """A weakly anchored stiff component must fail explicitly rather than use a pseudoinverse."""
    net = network(('a', 'x', 'y'), (('a', 'x', 1e-20), ('x', 'y', 1)), arithmetic='numeric')
    for operation in (dirichlet, kron_reduce):
        with pytest.raises(NumericalUnresolved):
            operation(net, ('a',))


def test_numeric_influence_underflow():
    """A well-conditioned block may still underflow a positive influence; do not erase it."""
    net = network(('a', 'x', 'y', 'b'), (('a', 'x', 1), ('x', 'y', 1e-200),
                  ('y', 'b', 1e-200), ('y', 'a', 1)), arithmetic='numeric')
    with pytest.raises(NumericalUnresolved):
        dirichlet(net, ('a', 'b'))


def test_numeric_zero_influence_is_topological():
    """The selected vertex p shields q on x-p-q, so B_xq is exactly zero."""
    net = network(('x', 'p', 'q'), (('x', 'p', 1), ('p', 'q', 1)), arithmetic='numeric')
    assert dirichlet(net, ('p', 'q')).value.B.values == ((1.0, 0.0),)


def test_numeric_contraction_overflow():
    """Two representable cross-group weights can have an unrepresentable sum."""
    net = network(('a', 'x', 'b'), (('a', 'b', 1e308), ('x', 'b', 1e308)), arithmetic='numeric')
    with pytest.raises(NumericalUnresolved):
        contract(net, {'u': ('a', 'x'), 'v': ('b',)})


def test_legacy_matrix_api_agrees():
    """New core operations preserve the established numeric matrix algorithms on a path."""
    from polyform_algebra import Basis, Space
    from polyform_algebra.matrix import Laplacian, MatrixAxis
    net = path('numeric')
    old = Laplacian(laplacian(net).value.to_numpy(), rows=MatrixAxis(basis=Basis(list(net.context.labels), space=Space('path', 4))))
    equal(old.dirichlet(['a', 'b']).influence_matrix.to_numpy(), dirichlet(net, ('a', 'b')).value.B.to_numpy(), 'numeric')
    equal(old.eliminate(['a', 'b']).laplacian.to_numpy(), laplacian(kron_reduce(net, ('a', 'b')).value.network).value.to_numpy(), 'numeric')
    partition = {'u': ('a', 'x'), 'v': ('y', 'b')}
    equal(old.contract(partition).laplacian.to_numpy(), laplacian(contract(net, partition).value.network).value.to_numpy(), 'numeric')
    equal(old.induced_subgraph(['x', 'y']).to_numpy(), laplacian(induced_subgraph(net, ('x', 'y')).value).value.to_numpy(), 'numeric')
