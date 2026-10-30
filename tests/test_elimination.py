"""Kron reduction against PMG examples and response/metric identities.

Source: PMG v5, note revision 5 (2026-10-15), read 2026-10-30.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
Equations (3)-(6), star/path examples, and 'Conditions and limitations'.
"""

import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import (
    EliminationResult, LabeledMatrix, Laplacian, MatrixAxis, MatrixKind, MatrixTolerance,
)


def graph(labels, edges, **kwargs):
    axis = MatrixAxis(basis=Basis(labels, space=Space("network", len(labels))))
    a = np.zeros((len(labels), len(labels)))
    for i, j, c in edges:
        a[i, i] += c
        a[j, j] += c
        a[i, j] -= c
        a[j, i] -= c
    return Laplacian(a, rows=axis, **kwargs)


def path():
    return graph(["a", "x", "y", "b"], [(0, 1, 2), (1, 2, 1), (2, 3, 3)])


def test_star_becomes_triangle_preserving_resistance_and_context():
    """Each induced edge is 1*1/3; parallel resistances 3 and 6 give 2."""
    lap = graph(["x", "a", "b", "c"], [(0, 1, 1), (0, 2, 1), (0, 3, 1)],
                kind=MatrixKind.OPERATOR)
    before = lap.to_numpy()
    result = lap.eliminate(keep=["c", "a", "b"])
    assert isinstance(result, EliminationResult)
    reduced = result.laplacian
    assert reduced.rows is result.dirichlet.influence_matrix.columns
    assert result.dirichlet.source_axis is lap.rows
    assert reduced.rows.labels == ("c", "a", "b")
    assert reduced.kind == lap.kind and reduced.tolerance is lap.tolerance
    np.testing.assert_allclose(reduced.to_numpy(), np.eye(3) - np.ones((3, 3))/3)
    np.testing.assert_allclose(result.dirichlet.influence_matrix.to_numpy(), [[1/3]*3])
    expected = 2 * (np.ones((3, 3)) - np.eye(3))
    np.testing.assert_allclose(reduced.to_resistance().to_numpy(), expected)
    np.testing.assert_allclose(lap.to_resistance().to_numpy()[np.ix_([3, 1, 2], [3, 1, 2])], expected)
    np.testing.assert_array_equal(lap.to_numpy(), before)


def test_path_response_energy_and_order_independence():
    """Series resistance 1/2+1+1/3=11/6; f=(0,3,9,11) gives E=66."""
    lap = path()
    result = lap.eliminate(["a", "b"])
    reduced = result.laplacian.to_numpy()
    np.testing.assert_allclose(reduced, np.array([[1, -1], [-1, 1]]) * 6/11)
    f = np.array(list(result.dirichlet.extend([0, 11]).values()))
    np.testing.assert_allclose(f, [0, 3, 9, 11])
    np.testing.assert_allclose(lap.to_numpy() @ f, [-6, 0, 0, 6], atol=1e-14)
    np.testing.assert_allclose(reduced @ [0, 11], [-6, 6])
    assert f @ lap.to_numpy() @ f == pytest.approx(66)
    assert np.array([0, 11]) @ reduced @ [0, 11] == pytest.approx(66)
    for intermediate in (["a", "y", "b"], ["a", "x", "b"]):
        sequential = lap.eliminate(intermediate).laplacian.eliminate(["a", "b"])
        np.testing.assert_allclose(sequential.laplacian.to_numpy(), reduced)


def test_existing_edges_receive_weighted_star_increment():
    """Equation (3) adds c_px*c_xq/d_x to existing couplings, not replacing them."""
    lap = graph(["x", "a", "b", "c"], [(0, 1, 2), (0, 2, 3), (0, 3, 5), (1, 2, 7)])
    adjacency = lap.eliminate(["a", "b", "c"]).laplacian.to_adjacency().to_numpy()
    np.testing.assert_allclose(adjacency, [[0, 7.6, 1], [7.6, 0, 1.5], [1, 1.5, 0]])


def test_single_retained_and_no_elimination():
    """A constant harmonic field has zero response; retaining all only reorders."""
    lap = path()
    one = lap.eliminate(["x"])
    np.testing.assert_array_equal(one.laplacian.to_numpy(), [[0]])
    assert one.dirichlet.extend([7]) == pytest.approx(dict.fromkeys(lap.rows.labels, 7))
    all_vertices = lap.eliminate(["b", "y", "x", "a"])
    np.testing.assert_array_equal(all_vertices.laplacian.to_numpy(), lap.to_numpy()[::-1, ::-1])
    assert all_vertices.dirichlet.influence_matrix.to_numpy().shape == (0, 4)
    singleton = graph(["o"], []).eliminate(["o"])
    np.testing.assert_array_equal(singleton.laplacian.to_numpy(), [[0]])


def test_components_isolates_and_unanchored_rejection():
    """Two separate edges reduce to two isolates, whose mutual resistance is infinite."""
    lap = graph(["a", "x", "b", "y", "o"], [(0, 1, 2), (2, 3, 3)])
    result = lap.eliminate(["o", "b", "a"])
    np.testing.assert_array_equal(result.laplacian.to_numpy(), np.zeros((3, 3)))
    r = result.laplacian.to_resistance().to_numpy()
    assert np.isposinf(r[~np.eye(3, dtype=bool)]).all()
    np.testing.assert_array_equal(np.diag(r), 0)
    np.testing.assert_array_equal(result.dirichlet.influence_matrix.to_numpy(), [[0, 0, 1], [0, 1, 0]])
    for keep in (["a"], ["a", "b"]):
        with pytest.raises(ValueError, match="no selected vertex"):
            lap.eliminate(keep)


def test_source_permutation_and_scale():
    """Relabeling order commutes with reduction; all conductances scale linearly."""
    lap = path()
    order = ["y", "b", "a", "x"]
    permuted = lap.reordered(rows=order, columns=order)
    np.testing.assert_allclose(permuted.eliminate(["b", "a"]).laplacian.to_numpy(),
                               lap.eliminate(["a", "b"]).laplacian.to_numpy()[::-1, ::-1])
    scaled = Laplacian(7 * lap.to_numpy(), rows=lap.rows)
    np.testing.assert_allclose(scaled.eliminate(["a", "b"]).laplacian.to_numpy(),
                               7 * lap.eliminate(["a", "b"]).laplacian.to_numpy())


def test_bipartite_special_case():
    """With diagonal internal block, reduction equals H diag(1/d) H.T off diagonal."""
    h = np.array([[2, 0], [1, 3], [0, 4.]])
    p = MatrixAxis(basis=Basis(["a", "b", "c"], space=Space("P", 3)))
    x = MatrixAxis(basis=Basis(["x", "y"], space=Space("X", 2)))
    table = LabeledMatrix(h, rows=p, columns=x, kind=MatrixKind.PAIRWISE, representation="bipartite")
    full = graph(["a", "b", "c", "x", "y"],
                 [(i, 3+j, h[i, j]) for i in range(3) for j in range(2)])
    np.testing.assert_allclose(full.eliminate(p.labels).laplacian.to_numpy(),
                               table.to_unipartite(p.space).to_numpy())


def test_large_singleton_response_has_no_cancellation_artifact():
    """One anchor has zero response even when subtracting large diagonal terms."""
    lap = path()
    large = Laplacian(lap.to_numpy() * 1e12, rows=lap.rows)
    np.testing.assert_array_equal(large.eliminate(["a"]).laplacian.to_numpy(), [[0]])


@pytest.mark.parametrize("keep", [[], ["a", "a"], ["missing"]])
def test_invalid_selection(keep):
    """The API requires a nonempty subset of distinct, known vertex labels."""
    with pytest.raises(ValueError):
        path().eliminate(keep)


@pytest.mark.parametrize("keep", ["a", {"a", "b"}])
def test_unordered_selection(keep):
    """Reject ambiguous ordering rather than silently assigning matrix axes."""
    with pytest.raises(TypeError):
        path().eliminate(keep)


def test_directed_empty_and_unresolved_cases():
    """Directed reduction is deferred; unresolved ranks never trigger a pseudoinverse."""
    axis = MatrixAxis(basis=Basis(["a", "b"], space=Space("V", 2)))
    with pytest.raises(NotImplementedError):
        Laplacian([[1, -1], [0, 0]], rows=axis).eliminate(["b"])
    with pytest.raises(ValueError):
        graph([], []).eliminate([])
    small = np.array([[1, -1], [-1, 1]]) * 1e-15
    with pytest.raises(ValueError, match="tolerance"):
        Laplacian(small, rows=axis).eliminate(["a"])
    resolved = Laplacian(small, rows=axis, tolerance=MatrixTolerance(atol=0))
    np.testing.assert_array_equal(resolved.eliminate(["a"]).laplacian.to_numpy(), [[0]])


def test_pandas_labels_and_stale_basis():
    """Response and influence export their own labeled axes; stale embeddings fail."""
    pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_dataframe
    lap = path()
    result = lap.eliminate(["b", "a"])
    frame = to_dataframe(result.laplacian)
    assert list(frame.index) == list(frame.columns) == ["b", "a"]
    assert frame.loc["a", "b"] == pytest.approx(-6/11)
    influence = to_dataframe(result.dirichlet.influence_matrix)
    assert list(influence.index) == ["x", "y"]
    lap.rows.basis.elements.reverse()
    with pytest.raises(ValueError, match="basis changed"):
        result.dirichlet.extend([11, 0])
