"""Deletion and contraction compared with PMG resistance-preserving elimination.

Sources read from v5 on 2026-10-30:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Contraction%20and%20Graph%20Gluing.md
Revision 6 (Russian revision 4), 'Three ways', equation (1), path and gluing.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
Revision 5, 'What should remain after elimination'.
"""

import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import ContractionResult, Laplacian, MatrixAxis, MatrixKind


def graph(labels, edges, **kwargs):
    axis = MatrixAxis(basis=Basis(labels, space=Space("graph", len(labels))))
    a = np.zeros((len(labels), len(labels)))
    for i, j, c in edges:
        a[i, i] += c
        a[j, j] += c
        a[i, j] -= c
        a[j, i] -= c
    return Laplacian(a, rows=axis, **kwargs)


def path():
    return graph(["a", "x", "y", "b"], [(0, 1, 2), (1, 2, 1), (2, 3, 3)])


def test_three_path_operations():
    """PMG path: deletion has R=inf, contraction R=1, elimination R=11/6."""
    lap = path()
    before = lap.to_numpy()
    deleted = lap.induced_subgraph(["a", "b"])
    contracted = lap.contract({"a": ["a", "x"], "b": ["y", "b"]})
    assert isinstance(contracted, ContractionResult)
    eliminated = lap.eliminate(["a", "b"]).laplacian
    np.testing.assert_array_equal(deleted.to_numpy(), np.zeros((2, 2)))
    assert np.isposinf(deleted.to_resistance().to_numpy()[0, 1])
    np.testing.assert_array_equal(contracted.laplacian.to_numpy(), [[1, -1], [-1, 1]])
    assert contracted.laplacian.to_resistance().to_numpy()[0, 1] == pytest.approx(1)
    assert eliminated.to_resistance().to_numpy()[0, 1] == pytest.approx(11/6)
    np.testing.assert_array_equal(lap.to_numpy(), before)


def test_induced_triangle_recomputes_degrees_and_keeps_order():
    """Deleting one unit-triangle vertex leaves degree 1, not the principal-block 2."""
    lap = graph(["a", "b", "c"], [(0, 1, 1), (1, 2, 1), (2, 0, 1)], kind=MatrixKind.OPERATOR)
    sub = lap.induced_subgraph(["b", "a"])
    assert sub.rows.labels == ("b", "a") and sub.rows is sub.columns
    assert sub.kind == lap.kind and sub.tolerance is lap.tolerance
    assert sub.rows.space is not lap.rows.space
    np.testing.assert_array_equal(sub.to_numpy(), [[1, -1], [-1, 1]])
    assert sub.to_resistance().to_numpy()[0, 1] == pytest.approx(1)


def test_induced_isolates_empty_and_all_vertices():
    """Deletion keeps selected isolates, permits empty graphs and needs no anchors."""
    lap = graph(["a", "x", "o"], [(0, 1, 2)])
    np.testing.assert_array_equal(lap.induced_subgraph(["o", "a"]).to_numpy(), np.zeros((2, 2)))
    empty = lap.induced_subgraph([])
    assert empty.to_numpy().shape == (0, 0) and empty.rows.labels == ()
    assert graph([], []).induced_subgraph([]).to_numpy().shape == (0, 0)
    np.testing.assert_array_equal(lap.induced_subgraph(["o", "x", "a"]).to_numpy(),
                                  lap.to_numpy()[::-1, ::-1])


def test_contraction_adds_couplings_discards_loops_and_preserves_energy():
    """Triangle contraction sums 2+3=5; internal weight 7 has zero constant-group energy."""
    lap = graph(["a", "b", "c"], [(0, 1, 2), (0, 2, 3), (1, 2, 7)], kind=MatrixKind.OPERATOR)
    result = lap.contract({"o": ["c", "b"], "a": ["a"]})
    s = result.membership_matrix.to_numpy()
    np.testing.assert_array_equal(s, [[0, 1], [1, 0], [1, 0]])
    assert result.membership_matrix.rows is lap.rows
    assert result.membership_matrix.columns is result.laplacian.rows
    assert result.membership_matrix.kind == MatrixKind.OPERATOR
    assert result.laplacian.kind == lap.kind and result.laplacian.tolerance is lap.tolerance
    reduced = result.laplacian.to_numpy()
    np.testing.assert_array_equal(reduced, [[5, -5], [-5, 5]])
    np.testing.assert_array_equal(reduced, s.T @ lap.to_numpy() @ s)
    f = np.array([4., 1.])
    assert (s @ f) @ lap.to_numpy() @ (s @ f) == pytest.approx(f @ reduced @ f)
    np.testing.assert_array_equal(s.sum(axis=1), 1)
    # S.T maps coefficients while preserving their affine sum.
    coefficients = np.array([0.2, 0.3, 0.5])
    assert (s.T @ coefficients).sum() == pytest.approx(coefficients.sum())


def test_nonadjacent_contraction_and_singletons():
    """PMG unit path endpoints need not be adjacent: their two couplings sum to 2."""
    lap = graph(["a", "b", "c"], [(0, 1, 1), (1, 2, 1)])
    reduced = lap.contract({"o": ["a", "c"], "b": ["b"]}).laplacian
    assert reduced.to_resistance().to_numpy()[0, 1] == pytest.approx(0.5)
    np.testing.assert_array_equal(lap.contract({"all": ["a", "b", "c"]}).laplacian.to_numpy(), [[0]])
    identity = lap.contract({"C": ["c"], "B": ["b"], "A": ["a"]})
    np.testing.assert_array_equal(identity.laplacian.to_numpy(), lap.to_numpy()[::-1, ::-1])
    assert identity.laplacian.rows.labels == ("C", "B", "A")
    assert graph([], []).contract({}).laplacian.to_numpy().shape == (0, 0)


def test_gluing_components_and_retaining_isolates():
    """Gluing endpoints of two unit edges creates a path of resistance 2, plus an isolate."""
    lap = graph(["a", "b", "c", "d", "z"], [(0, 1, 1), (2, 3, 1)])
    reduced = lap.contract({"a": ["a"], "o": ["b", "c"], "d": ["d"], "z": ["z"]}).laplacian
    resistance = reduced.to_resistance().to_numpy()
    assert resistance[0, 2] == pytest.approx(2)
    assert np.isposinf(resistance[0, 3])
    np.testing.assert_array_equal(reduced.to_numpy()[3], 0)
    separate = lap.contract({"u": ["a", "b"], "v": ["c", "d"], "z": ["z"]}).laplacian
    np.testing.assert_array_equal(separate.to_numpy(), np.zeros((3, 3)))


def test_contraction_order_and_source_permutation():
    """Successive identifications and simultaneous partition give the same weighted graph."""
    lap = path()
    expected = lap.contract({"left": ["a", "x"], "right": ["y", "b"]}).laplacian.to_numpy()
    first = lap.contract({"left": ["a", "x"], "y": ["y"], "b": ["b"]}).laplacian
    sequential = first.contract({"left": ["left"], "right": ["y", "b"]})
    np.testing.assert_array_equal(sequential.laplacian.to_numpy(), expected)
    order = ["y", "a", "b", "x"]
    reordered = lap.reordered(rows=order, columns=order)
    result = reordered.contract({"left": ["x", "a"], "right": ["b", "y"]})
    np.testing.assert_array_equal(result.laplacian.to_numpy(), expected)
    assert result.membership_matrix.rows.labels == tuple(order)


def test_large_internal_weights_do_not_cancel_external_edges():
    """Summing cross-group edges retains weight 1 beside an internal weight 1e16."""
    lap = graph(["a", "x", "b"], [(0, 1, 1e16), (1, 2, 1)])
    reduced = lap.contract({"left": ["a", "x"], "right": ["b"]}).laplacian
    np.testing.assert_array_equal(reduced.to_numpy(), [[1, -1], [-1, 1]])


def test_overflow_rejected_but_internal_weights_not_summed():
    """Two finite 1e308 cross edges overflow on aggregation; discarded internal edges do not."""
    lap = graph(range(4), [(0, 1, 1e308), (2, 3, 1e308)])
    with pytest.raises(ValueError, match="finite"):
        lap.contract({"left": [0, 2], "right": [1, 3]})
    np.testing.assert_array_equal(lap.contract({"all": [0, 1, 2, 3]}).laplacian.to_numpy(), [[0]])


@pytest.mark.parametrize("keep", [["a", "a"], ["unknown"]])
def test_invalid_induced_subset(keep):
    """Duplicate and unknown labels cannot define an induced vertex basis."""
    with pytest.raises(ValueError):
        path().induced_subgraph(keep)


@pytest.mark.parametrize("keep", ["a", {"a", "b"}])
def test_unordered_induced_subset(keep):
    """An explicit ordered iterable prevents accidental character or set ordering."""
    with pytest.raises(TypeError):
        path().induced_subgraph(keep)


@pytest.mark.parametrize("partition", [
    {}, {"u": ["a", "x"]}, {"u": []},
    {"u": ["a", "x"], "v": ["x", "y", "b"]},
    {"u": ["a", "a", "x", "y", "b"]},
    {"u": ["a", "x", "y", "b", "unknown"]},
])
def test_invalid_partition(partition):
    """A partition must cover the original vertex set exactly once in nonempty groups."""
    with pytest.raises(ValueError):
        path().contract(partition)


@pytest.mark.parametrize("partition", [[['a', 'x'], ['y', 'b']], {"u": "axyb"}, {"u": {"a", "x", "y", "b"}}])
def test_ambiguous_partition(partition):
    """Named groups and ordered member sequences prevent implicit labels or ordering."""
    with pytest.raises(TypeError):
        path().contract(partition)


def test_directed_and_invalid_sign_boundaries():
    """Directed variants are deferred; tolerated sign errors must not be silently clipped."""
    axis = MatrixAxis(basis=Basis(["a", "b"], space=Space("V", 2)))
    for values, error in [([[1, -1], [0, 0]], NotImplementedError),
                          ([[-1e-15, 1e-15], [1e-15, -1e-15]], ValueError)]:
        lap = Laplacian(values, rows=axis)
        with pytest.raises(error):
            lap.induced_subgraph(["a"])
        with pytest.raises(error):
            lap.contract({"o": ["a", "b"]})


def test_dataframe_tuple_labels_and_source_mutation():
    """Tuple labels stay atomic in exports and membership retains the validated source axis."""
    pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_dataframe
    lap = graph([("v", 1), ("v", 2), ("v", 3)], [(0, 1, 2), (1, 2, 3)])
    result = lap.contract({("group", 1): [("v", 1), ("v", 2)], ("group", 2): [("v", 3)]})
    frame = to_dataframe(result.laplacian)
    assert list(frame.index) == [("group", 1), ("group", 2)]
    sub = lap.induced_subgraph([("v", 2), ("v", 1)])
    assert list(to_dataframe(sub).columns) == [("v", 2), ("v", 1)]
    lap.rows.basis.elements.reverse()
    with pytest.raises(ValueError, match="basis changed"):
        result.membership_matrix.to_numpy()
