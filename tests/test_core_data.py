"""T12/T29/T33/T38: ordered spaces, rational coefficients and explicit conversion.

PMG v5 Core Specification, revision 27, underlying space and parallel conductors:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md
"""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F

import pytest
import sympy as sp

from polyform_algebra.core import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, InvalidInput,
    ContextMismatch, ArithmeticMismatch, NumericalUnresolved,
    coordinates, network, reorder, to_numeric,
)


def test_context_identity_and_order_are_not_names():
    """Independent spaces with identical labels cannot silently share coordinates."""
    first, second = Context(['a', 'b']), Context(['a', 'b'])
    assert first != second
    with pytest.raises(ContextMismatch):
        coordinates(first, [1, 2]) + coordinates(second, [1, 2])
    with pytest.raises(FrozenInstanceError):
        first.labels = ('b', 'a')


def test_coordinates_mapping_zero_and_arithmetic():
    """A zero keeps its space and division of rational coordinates remains exact."""
    labels = ['a', ('b', 1)]
    ctx = Context(labels, epsilon=[1, F(2, 3)])
    labels.reverse()
    x = coordinates(ctx, {('b', 1): sp.Rational(1, 3), 'a': 1})
    assert (x / 3).values == (F(1, 3), F(1, 9))
    assert (x-x).context is ctx
    assert (x-x).values == (0, 0)
    assert (2*x).values == (2, F(2, 3))
    assert ctx.labels == ('a', ('b', 1))


@pytest.mark.parametrize('labels', [{'a', 'b'}, ['a', 'a'], [True], [1, True], [('a', [])], 'ab'])
def test_reject_invalid_labels(labels):
    """Order and immutable identity are API requirements, not inferred from names."""
    with pytest.raises(InvalidInput):
        Context(labels)


@pytest.mark.parametrize('value', [1.0, True, float('nan'), float('inf'), sp.Symbol('c'), sp.sqrt(2)])
def test_exact_rejects_nonrational_inputs(value):
    """Exact calls must not silently rationalize floats, bools or symbolic expressions."""
    with pytest.raises(InvalidInput):
        coordinates(Context(['a']), [value])


@pytest.mark.parametrize('values', [{'a': 1}, {'a': 1, 'b': 2, 'c': 3}, [1], [[1], [2]], {1, 2}])
def test_coordinate_shape_and_labels(values):
    """Missing/extra labels and non-column containers are rejected, never filled."""
    with pytest.raises(InvalidInput):
        coordinates(Context(['a', 'b']), values)


def test_network_aggregation_isolates_and_reorder():
    """Parallel conductances add exactly; zero edges preserve vertices, not links."""
    edges = [('a', 'b', F(1, 3)), ('b', 'a', F(1, 6)), ('b', 'c', 0)]
    net = network(['a', 'b', 'c'], edges)
    edges.clear()
    assert net.edges == (('a', 'b', F(1, 2)),)
    assert net.components == (('a', 'b'), ('c',))
    changed = reorder(net, ['c', 'b', 'a'])
    assert changed.value.edges == (('b', 'a', F(1, 2)),)
    assert changed.value.components == (('c',), ('b', 'a'))
    assert changed.value.context.identity == net.context.identity
    assert changed.diagnostics.status == 'exact'
    assert net.context.labels == ('a', 'b', 'c')


@pytest.mark.parametrize('edges', [[('a', 'a', 1)], [('a', 'c', 1)], [('a', 'b', -1)], [('a', 'b')], [(True, 'a', 1)]])
def test_network_invalid_edges(edges):
    """Loops, unknown endpoints, negative weights and malformed edges violate C01."""
    with pytest.raises(InvalidInput):
        network(['a', 'b'], edges)


def test_parallel_rejection_and_direct_constructor_validation():
    """Reversed and zero parallel edges still count as duplicates under reject."""
    with pytest.raises(InvalidInput):
        network(['a', 'b'], [('a', 'b', 0), ('b', 'a', 1)], parallel='reject')
    with pytest.raises(InvalidInput):
        Network(Context(['a'], epsilon=[2]), ())
    with pytest.raises(InvalidInput):
        network(['a'], [], policy=NumericPolicy())
    assert network([], []).components == ()


def test_reorder_coordinates_and_matrix_preserves_values():
    """The same permutation applies to coefficients and both network matrix axes."""
    ctx = Context(['a', 'b'])
    x = coordinates(ctx, [F(1, 3), F(1, 7)])
    mat = MatrixData([[1, -1], [-1, 1]], ctx, ctx, 'laplacian')
    permuted = reorder(x, ['b', 'a']).value
    assert permuted.values == tuple(reversed(x.values))
    assert reorder(permuted, ['a', 'b']).value == x
    assert reorder(mat, ['b', 'a']).value.values == mat.values
    exported = mat.to_numpy()
    exported[0, 0] = 99
    assert mat.values[0][0] == 1
    with pytest.raises(ContextMismatch):
        x + permuted
    with pytest.raises(ContextMismatch):
        reorder(mat, ['a'])


def test_numeric_conversion_is_explicit_and_diagnosed():
    """C05 preserves space identity while changing arithmetic only by explicit call."""
    net = network(['a', 'b'], [('a', 'b', F(1, 3))])
    result = to_numeric(net)
    assert result.value.context.identity == net.context.identity
    assert result.value.edges[0][2] == float(F(1, 3))
    assert result.diagnostics.status == 'approximate'
    x = coordinates(net.context, [1, -1])
    numeric = to_numeric(x).value
    with pytest.raises(ArithmeticMismatch):
        x + numeric
    assert (numeric / 2).values == (.5, -.5)
    mat = MatrixData([[1, -1], [-1, 1]], net.context, net.context, 'laplacian')
    assert to_numeric(mat).value.to_numpy().dtype.kind == 'f'


@pytest.mark.parametrize('value', [F(10**400), F(1, 10**400)])
def test_conversion_rejects_overflow_and_underflow(value):
    """A finite rational cannot become infinity or a false zero through C05."""
    x = coordinates(Context(['a']), [value])
    with pytest.raises((NumericalUnresolved, InvalidInput)):
        to_numeric(x)
    assert x.values == (value,)


def test_matrix_empty_rectangular_and_context_validation():
    """Containers retain axes even for zero rows; they do not certify network laws."""
    ctx = Context([])
    assert MatrixData([], ctx, ctx, 'laplacian').to_numpy().shape == (0, 0)
    other = Context(['a'])
    assert MatrixData([], ctx, other, 'gram').to_numpy().shape == (0, 1)
    with pytest.raises(ContextMismatch):
        MatrixData([[0]], other, Context(['a']), 'laplacian')
    with pytest.raises(InvalidInput):
        MatrixData([[0, 0]], other, other, 'laplacian')


@pytest.mark.parametrize('bad', [-1, True, float('nan')])
def test_numeric_policy_rejects_invalid_tolerance(bad):
    """Tolerance validation is explicit and does not alter connectivity."""
    with pytest.raises(InvalidInput):
        NumericPolicy(rtol=bad)


def test_small_numeric_edge_is_never_pruned_by_policy():
    """Topology uses positive conductance even below any requested absolute tolerance."""
    net = network(['a', 'b'], [('a', 'b', 1e-200)], arithmetic='numeric',
                  policy=NumericPolicy(conductance_atol=1))
    assert net.components == (('a', 'b'),)
    assert net.edges[0][2] == 1e-200


def test_all_permutations_preserve_labeled_network_and_coordinates():
    """T12: inverse permutations recover exact values and the original space identity."""
    from itertools import permutations
    net = network(['a', 'b', 'c'], [('a', 'b', sp.Rational(2, 3)), ('b', 'c', sp.Integer(2))])
    x = coordinates(net.context, [F(1, 7), F(2, 7), -F(3, 7)])
    for order in permutations(net.context.labels):
        assert reorder(reorder(net, order).value, net.context.labels).value == net
        assert reorder(reorder(x, order).value, net.context.labels).value == x
