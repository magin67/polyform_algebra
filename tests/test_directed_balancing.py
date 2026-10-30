"""Positive source scaling of directed Laplacians.

PMG v5 revision e28771c349cb54d84040702e7dc0cd345b57fd0f:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balancing%20a%20Directed%20Graph%20by%20Spanning-Tree%20Potentials.md
Equations (1)-(5), three-participant example, and impossibility for arcs
between distinct strongly connected components.
"""

import networkx as nx
import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis, MatrixKind, MatrixTolerance
from polyform_algebra.matrix.networkx import from_networkx


def example():
    axis = MatrixAxis(basis=Basis(["a", "b", "c"], space=Space("participants", 3)))
    return Laplacian([[3, -1, -2], [-2, 3, -1], [-1, -1, 2]], rows=axis,
                      kind=MatrixKind.OPERATOR)


def test_note_example_and_independent_normalized_solve():
    """Tau=(5,4,7) gives sigma=3*tau/16 and balanced adjacency 3/16 times the note's matrix."""
    lap = example()
    result = lap.balance()
    expected = 3/16 * np.array([[0, 5, 10], [8, 0, 4], [7, 7, 0]])
    np.testing.assert_allclose(result.multipliers.to_numpy(), [15/16, 12/16, 21/16])
    np.testing.assert_allclose(result.laplacian.to_adjacency().to_numpy(), expected)
    # Equation (6) independently solves normalization without the SVD basis.
    solution = np.linalg.solve(lap.to_numpy().T + np.ones((3, 3))/3, np.ones(3))
    np.testing.assert_allclose(result.multipliers.to_numpy(), solution)
    assert sum(result.multipliers.to_dict().values()) == pytest.approx(3)
    assert result.laplacian.rows is lap.rows and result.laplacian.columns is lap.columns
    assert result.laplacian.kind == lap.kind
    assert result.laplacian.vector_type() == "vector"
    assert not result.laplacian.is_symmetric()
    np.testing.assert_array_equal(lap.to_numpy(), [[3, -1, -2], [-2, 3, -1], [-1, -1, 2]])


def test_support_allocation_recovery_and_tree_potentials():
    """Row scaling preserves allocation and can be undone; all new rooted potentials equal 315/64."""
    original = example()
    result = original.balance()
    a = original.to_adjacency().to_numpy()
    b = result.laplacian.to_adjacency().to_numpy()
    np.testing.assert_array_equal(a > 0, b > 0)
    np.testing.assert_allclose(a/a.sum(axis=1)[:, None], b/b.sum(axis=1)[:, None])
    np.testing.assert_allclose(b/result.multipliers.to_numpy()[:, None], a)
    np.testing.assert_allclose(b.sum(axis=1), [45/16, 36/16, 42/16])
    np.testing.assert_allclose(b.sum(axis=1), b.sum(axis=0))
    # The mean-one normalization fixes multipliers, not the original total 8.
    assert b.sum() == pytest.approx(123/16)
    l = result.laplacian.to_numpy()
    for i in range(3):
        assert np.linalg.det(np.delete(np.delete(l, i, axis=0), i, axis=1)) == pytest.approx(315/64)


def test_identity_and_idempotence():
    """A balanced graph needs unit source multipliers, so balancing twice changes nothing."""
    cycle = from_networkx(nx.DiGraph([(0, 1), (1, 2), (2, 0)]))
    for lap in (cycle, from_networkx(nx.path_graph(3)), example().balance().laplacian):
        result = lap.balance()
        np.testing.assert_array_equal(result.multipliers.to_numpy(), np.ones(3))
        np.testing.assert_array_equal(result.laplacian.to_numpy(), lap.to_numpy())


def test_reordering_and_common_scale():
    """Source multipliers depend on relative weights, and vertex permutations only reorder them."""
    lap = example()
    order = ["c", "a", "b"]
    reordered = lap.reordered(rows=order, columns=order).balance()
    assert reordered.multipliers.labels == tuple(order)
    assert reordered.multipliers.to_dict() == pytest.approx(lap.balance().multipliers.to_dict())
    scaled = Laplacian(7*lap.to_numpy(), rows=lap.rows).balance()
    np.testing.assert_allclose(scaled.multipliers.to_numpy(), lap.balance().multipliers.to_numpy())
    np.testing.assert_allclose(scaled.laplacian.to_numpy(), 7*lap.balance().laplacian.to_numpy())


@pytest.mark.parametrize("edges", [[(0, 1)], [(0, 1), (0, 2)],
                                  [(0, 1), (1, 0), (1, 2), (2, 3), (3, 2)]])
def test_non_strong_components_rejected_without_deleting_arcs(edges):
    """Arcs between SCCs cannot be retained by positive balancing, even with a unique left vector."""
    lap = from_networkx(nx.DiGraph(edges))
    before = lap.to_numpy()
    for componentwise in (False, True):
        with pytest.raises(ValueError, match="strongly connected"):
            lap.balance(componentwise=componentwise)
    np.testing.assert_array_equal(lap.to_numpy(), before)


def test_disconnected_components_and_isolate():
    """Independent weak components have independent mean-one scales; isolated points use multiplier one."""
    graph = nx.DiGraph()
    graph.add_nodes_from(["a", "x", "b", "y", "o"])
    graph.add_weighted_edges_from([("a", "b", 1), ("b", "a", 2),
                                   ("x", "y", 3), ("y", "x", 3)])
    lap = from_networkx(graph)
    with pytest.raises(ValueError, match="componentwise"):
        lap.balance()
    results = lap.balance(componentwise=True)
    assert [r.multipliers.labels for r in results] == [("a", "b"), ("x", "y"), ("o",)]
    # For the two-way edge: sigma_a*1=sigma_b*2, sigma_a+sigma_b=2.
    np.testing.assert_allclose(results[0].multipliers.to_numpy(), [4/3, 2/3])
    np.testing.assert_array_equal(results[1].multipliers.to_numpy(), [1, 1])
    assert results[2].multipliers.to_dict() == {"o": 1}
    for result in results:
        assert result.laplacian.vector_type() == "vector"


def test_empty_singleton_and_invalid_mode():
    """Empty default has no connected domain; singleton identity extends the note's n>=2 scope."""
    empty = from_networkx(nx.Graph())
    with pytest.raises(ValueError, match="nonempty"):
        empty.balance()
    assert empty.balance(componentwise=True) == ()
    with pytest.raises(TypeError):
        empty.balance(componentwise="yes")
    singleton = from_networkx(nx.empty_graph(1)).balance()
    assert singleton.multipliers.to_dict() == {0: 1}
    np.testing.assert_array_equal(singleton.laplacian.to_numpy(), [[0]])


def test_unresolved_kernel_and_sign_corrections_rejected():
    """Numerical rank ambiguity and negative conductances cannot silently become balancing choices."""
    graph = nx.DiGraph()
    graph.add_weighted_edges_from([(0, 1, 1e-15), (1, 0, 2e-15)])
    lap = from_networkx(graph)
    with pytest.raises(ValueError, match="tolerance"):
        lap.balance()
    precise = Laplacian(lap.to_numpy(), rows=lap.rows, tolerance=MatrixTolerance(atol=0))
    np.testing.assert_allclose(precise.balance().multipliers.to_numpy(), [4/3, 2/3])
    # Generic storage accepts these sign errors at default tolerance, but
    # support-preserving scaling must not silently clip them into a new graph.
    signed = Laplacian([[-1e-12, 1e-12], [1e-12, -1e-12]], rows=lap.rows)
    with pytest.raises(ValueError, match="sign corrections"):
        signed.balance()


def test_dataframe_balance_check():
    """Pandas tables preserve original vertex names and display equal incoming/outgoing totals."""
    pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_dataframe, to_series
    result = example().balance()
    frame = to_dataframe(result.laplacian.to_adjacency())
    assert list(frame.index) == list(frame.columns) == ["a", "b", "c"]
    np.testing.assert_allclose(frame.sum(axis=1), frame.sum(axis=0))
    assert to_series(result.multipliers).name == "balancing_multiplier"
