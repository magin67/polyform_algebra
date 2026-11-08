"""Exact/numeric C02/C03 M01-M04 G01-G04; reference specification T01-T15.

Source: PMG v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0,
https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Introduction
Laplacian/Green transformations; vector inner products; lengths/areas/volumes.
Independent references below use grounded current solves and tree enumeration.
"""
from fractions import Fraction as F
from itertools import combinations, permutations

import numpy as np
import pytest
import sympy as sp

from polyform_algebra.core import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, INFINITY,
    InvalidInput, DomainError, NumericalUnresolved, ContextMismatch, ArithmeticMismatch,
    network, coordinates, reorder, to_numeric, laplacian, green, resistance,
    spanning_tree_weight, network_from_laplacian, network_from_resistance,
    gram, inner_product, blade_inner_product, simplex_volume_squared,
)


def arr(calc):
    return calc.value.to_numpy()


@pytest.mark.parametrize('mode', ['exact', 'numeric'])
def test_weighted_chain_matrices(mode):
    """T02: chain weights 1,2 give tau=2 and series resistance 1+1/2."""
    net = network(['a','b','c'], [('a','b',1),('b','c',2)], arithmetic=mode)
    g = [[F(1,2),-F(1,6),-F(1,3)],[-F(1,6),F(1,6),0],[-F(1,3),0,F(1,3)]]
    r = [[0,1,F(3,2)],[1,0,F(1,2)],[F(3,2),F(1,2),0]]
    if mode == 'exact':
        assert green(net).value.values == tuple(map(tuple,g))
        assert resistance(net).value.values == tuple(map(tuple,r))
    else:
        np.testing.assert_allclose(arr(green(net)), np.array(g,dtype=float), atol=1e-14)
        np.testing.assert_allclose(arr(resistance(net)), np.array(r,dtype=float), atol=1e-14)
        assert green(net).diagnostics.numerical[0].rank == 2
    assert spanning_tree_weight(net).value == pytest.approx(2)
    l = laplacian(net).value.values
    recovered = network_from_laplacian(l, ['a','b','c'], arithmetic=mode).value
    assert recovered.edges == net.edges
    if mode == 'exact':
        assert network_from_resistance(r,['a','b','c']).value.edges == net.edges


@pytest.mark.parametrize('mode', ['exact', 'numeric'])
def test_triangle_and_gram_geometry(mode):
    """T03/T04: triangle tree sum is 26; unit-triangle area squared is 1/12."""
    weighted = network(['a','b','c'], [('a','b',2),('b','c',3),('a','c',4)], arithmetic=mode)
    assert spanning_tree_weight(weighted).value == pytest.approx(26)
    pairs = resistance(weighted, pairs=[('a','b'),('a','c'),('b','c')]).value
    assert pairs.values == pytest.approx([7/26,5/26,6/26])
    restored = network_from_resistance(resistance(weighted).value.values, ['a','b','c'], arithmetic=mode).value
    np.testing.assert_allclose(np.array(laplacian(restored).value.values,dtype=float), arr(laplacian(weighted)).astype(float), atol=1e-12)
    net = network(['a','b','c'], [('a','b',1),('b','c',1),('a','c',1)], arithmetic=mode)
    u,v = coordinates(net.context,[-1,1,0]), coordinates(net.context,[-1,0,1])
    np.testing.assert_allclose(arr(gram(net,[u,v])).astype(float), [[2/3,1/3],[1/3,2/3]], atol=1e-14)
    assert inner_product(net,u,v).value == pytest.approx(1/3)
    assert blade_inner_product(net,[u,v],[v,u]).value == pytest.approx(-1/3)
    assert blade_inner_product(net,[u,u],[u,u]).value == pytest.approx(0,abs=1e-15)
    for order in permutations(net.context.labels):
        assert simplex_volume_squared(net,list(order)).value == pytest.approx(1/12)
    assert inner_product(net,2*u-v,v).value == pytest.approx(0,abs=1e-14)


@pytest.mark.parametrize('mode', ['exact', 'numeric'])
def test_disconnected_empty_singleton_and_infinity(mode):
    """T08/T09: infinity means disconnection; a singleton contributes tree weight one."""
    net = network(['a','b','o'], [('a','b',2)], arithmetic=mode)
    r = resistance(net).value
    assert r.values[0][1] == pytest.approx(1/2)
    infinity = INFINITY if mode == 'exact' else float('inf')
    assert r.values[0][2] == infinity
    assert spanning_tree_weight(net).value == 0
    assert [x[1] for x in spanning_tree_weight(net,componentwise=True).value] == pytest.approx([2,1])
    restored = network_from_resistance(r.values,net.context.labels,arithmetic=mode).value
    assert restored.components == net.components
    with pytest.raises(DomainError):
        inner_product(net,coordinates(net.context,[1,0,-1]),coordinates(net.context,[1,0,-1]))
    with pytest.raises(DomainError):
        simplex_volume_squared(net,['a','o'])
    for labels,tau in [([],0),(['x'],1)]:
        tiny = network(labels,[],arithmetic=mode)
        assert spanning_tree_weight(tiny).value == tau
        for func in (laplacian,green,resistance):
            assert arr(func(tiny)).shape == (len(labels),len(labels))
    assert simplex_volume_squared(network(['x'],[],arithmetic=mode),['x']).value == 1
    assert blade_inner_product(network([],[],arithmetic=mode),[],[]).value == 1
    if mode == 'exact':
        assert to_numeric(r).value.values[0][2] == float('inf')


@pytest.mark.parametrize('c',[1e-12,1.,1e12])
def test_electrical_scale_and_single_edge(c):
    """T01/T28: inverse scales are handled per component; R=1/c remains finite."""
    net=network(['a','b'],[('a','b',c)],arithmetic='numeric')
    assert resistance(net).value.values[0][1] == pytest.approx(1/c,rel=1e-9,abs=0)
    np.testing.assert_allclose(arr(green(net)),np.array([[1,-1],[-1,1]])/(4*c),rtol=1e-9,atol=0)
    assert spanning_tree_weight(net).value == pytest.approx(c,rel=1e-9,abs=0)


def test_exact_scale_and_permutation_laws():
    """T10/T12/T13: scaling and relabeling commute with all matrix/geometric operations."""
    labels=['a','b','c']
    edges=[('a','b',2),('b','c',3),('a','c',4)]
    net=network(labels,edges)
    t=F(7,3)
    scaled=network(labels,[(a,b,t*c) for a,b,c in edges])
    for func,factor in [(laplacian,t),(green,1/t),(resistance,1/t)]:
        assert sp.Matrix(func(scaled).value.values) == factor*sp.Matrix(func(net).value.values)
    assert spanning_tree_weight(scaled).value == t**2*spanning_tree_weight(net).value
    assert simplex_volume_squared(scaled,labels).value == simplex_volume_squared(net,labels).value/t**2
    for order in permutations(labels):
        changed=reorder(net,order).value
        for func in (laplacian,green,resistance):
            assert func(changed).value == reorder(func(net).value,order).value
    l,g=sp.Matrix(laplacian(net).value.values),sp.Matrix(green(net).value.values)
    j=sp.eye(3)-sp.ones(3)/3
    assert l*g == g*l == j


def test_invalid_metric_arguments_and_pairs():
    """T33: contexts, grade, labels and zero-sum constraints are never guessed."""
    net=network(['a','b'],[('a','b',1)])
    u=coordinates(net.context,[-1,1])
    with pytest.raises(ContextMismatch):
        inner_product(net,u,coordinates(Context(['a','b']),[-1,1]))
    with pytest.raises(ArithmeticMismatch):
        inner_product(net,u,to_numeric(u).value)
    with pytest.raises(DomainError):
        inner_product(net,u,coordinates(net.context,[1,0]))
    with pytest.raises(DomainError):
        blade_inner_product(net,[u],[])
    for pairs in [[('a','z')],[('a','b'),('a','b')],[(True,'b')]]:
        with pytest.raises(InvalidInput):
            resistance(net,pairs=pairs)
    assert resistance(net,pairs=[]).value.values == ()


@pytest.mark.parametrize('values',[
    [[1,-1],[-2,2]], [[2,-1],[-1,1]], [[-1,1],[1,-1]], [[1,0,0]],
])
def test_invalid_laplacians(values):
    """T34: shape, symmetry, sums and conductance signs are exact input contracts."""
    with pytest.raises(InvalidInput):
        network_from_laplacian(values,['a','b'])


@pytest.mark.parametrize('values',[
    [[0,0],[0,0]], [[1,1],[1,0]], [[0,-1],[-1,0]], [[0,1],[2,0]],
    [[0,1,INFINITY],[1,0,1],[INFINITY,1,0]],
    [[0,1,9],[1,0,4],[9,4,0]],
    [[0,1,5],[1,0,2],[5,2,0]],
])
def test_invalid_resistance(values):
    """T11/T34: not every distance matrix is a connected resistive network metric."""
    with pytest.raises(InvalidInput):
        network_from_resistance(values,list(range(len(values))))


def test_numeric_corrections_and_uncertain_topology():
    """Only permitted symmetry/diagonal repairs are reported; signs are not repaired."""
    calc=network_from_laplacian([[1+1e-12,-1],[-1,1]],['a','b'],arithmetic='numeric')
    assert calc.diagnostics.corrections == (('rebuilt diagonal',0),)
    with pytest.raises(InvalidInput):
        network_from_laplacian([[-1e-20,1e-20],[1e-20,-1e-20]],['a','b'],arithmetic='numeric')
    weak=network([0,1,2,3],[(0,1,1),(1,2,1e-15),(2,3,1)],arithmetic='numeric')
    with pytest.raises(NumericalUnresolved):
        green(weak)
    assert weak.components == ((0,1,2,3),)
    path=network([0,1,2],[(0,1,1),(1,2,2)],arithmetic='numeric')
    with pytest.raises(NumericalUnresolved):
        network_from_resistance(resistance(path).value.values,[0,1,2],arithmetic='numeric')


def _tree_reference(n, edges):
    """Enumerate n-1 edge subsets, multiplying their weights when connected."""
    if n==0:return 0
    answer=0
    for chosen in combinations(edges,n-1):
        reached={0}
        for _ in range(n):
            for a,b,_ in chosen:
                if a in reached or b in reached:reached.update((a,b))
        if len(reached)==n:
            weight=1
            for _,_,c in chosen:weight*=c
            answer+=weight
    return answer


def test_all_simple_graphs_up_to_four_vertices():
    """T40: enumerate every small graph; compare trees and grounded current solves."""
    for n in range(5):
        possible=list(combinations(range(n),2))
        for mask in range(1<<len(possible)):
            edges=[(a,b,1) for k,(a,b) in enumerate(possible) if mask>>k&1]
            net=network(list(range(n)),edges)
            assert spanning_tree_weight(net).value == _tree_reference(n,edges)
            r=resistance(net).value.values
            assert network_from_resistance(r,list(range(n))).value.edges == net.edges
            if n>1 and len(net.components)==1:
                l=sp.zeros(n)
                for a,b,c in edges:
                    l[a,a]+=c;l[b,b]+=c;l[a,b]-=c;l[b,a]-=c
                for a,b in possible:
                    current=sp.zeros(n-1,1)
                    if a<n-1:current[a]+=1
                    if b<n-1:current[b]-=1
                    potential=l[:-1,:-1].inv()*current
                    expected=(potential[a] if a<n-1 else 0)-(potential[b] if b<n-1 else 0)
                    assert r[a][b] == expected


def test_twenty_pinned_audit_networks():
    """T39: explicit audit graphs and grounded rational answers survive both backends."""
    import json
    from pathlib import Path
    cases=json.loads((Path(__file__).parent/'fixtures/audit_networks.json').read_text())
    assert len(cases)==20
    for case in cases:
        expected=tuple(tuple(F(x) for x in row) for row in case['resistance'])
        for mode in ('exact','numeric'):
            net=network(list(range(case['n'])),case['edges'],arithmetic=mode)
            actual=resistance(net).value.values
            if mode=='exact':
                assert actual==expected
                assert network_from_resistance(actual,net.context.labels).value.edges==net.edges
            else:
                np.testing.assert_allclose(np.array(actual),np.array(expected,dtype=float),rtol=1e-9,atol=1e-14)
            assert spanning_tree_weight(net).value == pytest.approx(float(case['tau']),rel=1e-9)


def test_numeric_independent_component_scales():
    """A large scale in one component cannot erase conductances in another."""
    net=network([0,1,2,3,4],[(0,1,1e-12),(2,3,1e12)],arithmetic='numeric')
    result=resistance(net).value.values
    assert result[0][1] == pytest.approx(1e12,rel=1e-9)
    assert result[2][3] == pytest.approx(1e-12,rel=1e-9,abs=0)
    assert np.isinf(result[0][2])
    assert [c.rank for c in green(net).diagnostics.numerical] == [1,1,0]


def test_numeric_metric_range_and_dependence():
    """Nonzero norms cannot silently underflow; proportional directions have zero blade."""
    net=network([0,1,2],[(0,1,1),(1,2,1)],arithmetic='numeric')
    u=coordinates(net.context,[-1,1,0])
    assert blade_inner_product(net,[u,3*u],[u,3*u]).value == 0
    tiny=coordinates(net.context,[-1e-200,1e-200,0])
    with pytest.raises(NumericalUnresolved):
        inner_product(net,tiny,tiny)


def test_single_edge_exact_and_empty_metric_requests():
    """T01/T09: c=2 gives Gr entries +/-1/8; empty vector metrics are undefined."""
    net=network(['a','b'],[('a','b',2)])
    assert green(net).value.values==((F(1,8),-F(1,8)),(-F(1,8),F(1,8)))
    empty=network([],[])
    with pytest.raises(DomainError):
        gram(empty,[])
    with pytest.raises(DomainError):
        simplex_volume_squared(empty,[])


@pytest.mark.parametrize('t',[1e-12,1e12])
def test_numeric_area_scaling(t):
    """T13: scaling all triangle conductances by t scales area squared by t^-2."""
    net=network([0,1,2],[(0,1,t),(1,2,t),(0,2,t)],arithmetic='numeric')
    assert simplex_volume_squared(net,[0,1,2]).value == pytest.approx(1/(12*t*t),rel=1e-9,abs=0)


@pytest.mark.parametrize('mode',['exact','numeric'])
def test_tetrahedron_volume_normalization(mode):
    """K4 has tau=16, edge-Gram determinant 1/16 and volume squared 1/(16*36)."""
    net=network(list(range(4)),[(a,b,1) for a,b in combinations(range(4),2)],arithmetic=mode)
    assert spanning_tree_weight(net).value == pytest.approx(16)
    assert simplex_volume_squared(net,list(range(4))).value == pytest.approx(1/576)
