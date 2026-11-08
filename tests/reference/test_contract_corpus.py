"""Independent acceptance references for core geometry, T10-T15 and T39-T40.

PMG sources and explicit case data: ../fixtures/reference_corpus_v1.json.
The oracles module never imports PMG; expected answers do not call the subject.
"""
import json
from fractions import Fraction as Q
from itertools import combinations, permutations
from pathlib import Path

import numpy as np
import pytest
import sympy as sp

from polyform_algebra.core import (
    network, laplacian, green, resistance, spanning_tree_weight, Coordinates,
    blade_inner_product, simplex_volume_squared,
)
from . import oracles


@pytest.mark.parametrize('weight', [Q(1, 10**6), Q(1, 10**14)])
def test_weak_bridge_is_resolved_or_reported(weight):
    """T32: a nonzero bridge gives series resistance 2+1/w, never disconnection."""
    from polyform_algebra.core import NumericalUnresolved
    e = ((0, 1, 1), (1, 2, weight), (2, 3, 1))
    expected = 2 + 1/weight
    assert resistance(network(range(4), e)).value.values[0][3] == expected
    try:
        actual = resistance(network(range(4), e, arithmetic='numeric')).value.values[0][3]
    except NumericalUnresolved:
        return
    assert np.isfinite(actual)
    assert actual == pytest.approx(float(expected), rel=1e-7)

FIXTURES = Path(__file__).resolve().parents[1] / 'fixtures'
CORPUS = json.loads((FIXTURES / 'reference_corpus_v1.json').read_text(encoding='utf-8'))
AUDIT = json.loads((FIXTURES / 'audit_networks.json').read_text(encoding='utf-8'))


def case(name):
    return next(c for c in CORPUS['cases'] if c['id'] == name)


def edges(c):
    return tuple((a, b, Q(w)) for a, b, w in c['edges'])


@pytest.mark.parametrize('c', AUDIT, ids=['audit-{:02}'.format(i+1) for i in range(20)])
def test_audit_answers_against_two_forest_enumeration(c):
    """T39: recompute every pinned R using separating forests, not the original grounded solver."""
    tau, two_forests = oracles.tree_and_separating_forests(c['n'], edges(c))
    assert tau == Q(c['tau'])
    expected = tuple(tuple(x/tau for x in row) for row in two_forests)
    assert expected == tuple(tuple(Q(x) for x in row) for row in c['resistance'])
    net = network(tuple(range(c['n'])), edges(c))
    assert resistance(net).value.values == expected
    assert spanning_tree_weight(net).value == tau
    numeric = network(net.context.labels, edges(c), arithmetic='numeric')
    actual = resistance(numeric).value.to_numpy()
    for i in range(c['n']):
        assert actual[i, i] == 0
        for j in range(i+1, c['n']):
            assert actual[i, j] == pytest.approx(float(expected[i][j]), rel=1e-9, abs=0)


@pytest.mark.parametrize('mode', ['exact', 'numeric'])
@pytest.mark.parametrize('scale', [Q(1, 10**12), Q(1), Q(10**12)])
def test_t15_volume_orientation_basepoint_and_scale(mode, scale):
    """T15: all 24 basepoint/orientation choices match an independent Cayley-Menger determinant."""
    c = case('weighted_tetrahedron')
    e = tuple((a, b, scale*w) for a, b, w in edges(c))
    r = oracles.grounded_resistance(4, e)
    expected = oracles.cayley_menger(r, (0, 1, 2, 3))
    net = network((0, 1, 2, 3), e, arithmetic=mode)
    for order in permutations(range(4)):
        value = simplex_volume_squared(net, order).value
        if mode == 'exact':
            assert value == expected
        else:
            assert value == pytest.approx(float(expected), rel=1e-9, abs=0)
    u, v, w = [Coordinates(net.context, tuple(int(i == j)-int(i == 0) for i in range(4))) for j in (1, 2, 3)]
    original = blade_inner_product(net, (u, v, w), (u, v, w)).value
    reversed_value = blade_inner_product(net, (v, u, w), (u, v, w)).value
    assert reversed_value == (-original if mode == 'exact' else pytest.approx(-original, rel=1e-9, abs=0))
    assert blade_inner_product(net, (u, v, u+2*v), (u, v, u+2*v)).value == 0
    # Column transformation [[1,2,0],[0,1,1],[1,0,1]] has determinant 3.
    transformed = (u+w, 2*u+v, v+w)
    value = blade_inner_product(net, transformed, transformed).value
    assert value == (9*original if mode == 'exact' else pytest.approx(9*original, rel=1e-9, abs=0))


def test_all_small_graph_components_and_projectors():
    """T10/T40: all 76 graphs have the independently traversed components and block projectors."""
    count = 0
    for n in range(5):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            e = tuple((a, b, 1) for i, (a, b) in enumerate(possible) if mask & (1 << i))
            neighbors = [set() for _ in range(n)]
            for a, b, _ in e:
                neighbors[a].add(b)
                neighbors[b].add(a)
            unseen, groups = set(range(n)), []
            while unseen:
                found, pending = set(), [min(unseen)]
                while pending:
                    v = pending.pop()
                    if v not in found:
                        found.add(v)
                        pending.extend(neighbors[v]-found)
                unseen -= found
                groups.append(tuple(sorted(found)))
            net = network(tuple(range(n)), e)
            assert net.components == tuple(groups)
            l = sp.Matrix(laplacian(net).value.to_numpy())
            assert l == oracles.laplacian(n, e)
            g = sp.Matrix(green(net).value.to_numpy())
            j = sp.zeros(n)
            for group in groups:
                for a in group:
                    for b in group:
                        j[a, b] = int(a == b)-sp.Rational(1, len(group))
            assert l*g == g*l == j
            count += 1
    assert count == 76


def test_reference_metadata_and_oracle_independence():
    """Sources, explicit cases and PythonExamples snapshots stay traceable without network access."""
    import ast
    tree = ast.parse(Path(oracles.__file__).read_text(encoding='utf-8'))
    imports = [node.module for node in ast.walk(tree) if isinstance(node, ast.ImportFrom)]
    imports += [alias.name for node in ast.walk(tree) if isinstance(node, ast.Import) for alias in node.names]
    assert not any(name and name.startswith('polyform_algebra') for name in imports)
    assert CORPUS['schema'] == 1 and len(CORPUS['notes_revision']) == 40
    assert len({c['id'] for c in CORPUS['cases']}) == len(CORPUS['cases'])
    for source in CORPUS['sources']:
        assert CORPUS['notes_revision'] in source['url'] and len(source['blob']) == 40
    pages = json.loads((FIXTURES / 'pythonexamples_v5.json').read_text(encoding='utf-8'))
    assert {p['path']: p['blob'] for p in pages} == {p['path']: p['blob'] for p in CORPUS['pythonexamples']}
