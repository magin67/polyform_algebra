"""PMG adjacency eigenvector centrality, not an eigenvector of the Laplacian.

Source: v5 revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Centroid,%20Spanning-Tree,%20and%20Eigenvector%20Centralities.md
Sections 'Eigenvector centrality and a common normalization' and
'The complete graph, a path, and a star'.
"""

import networkx as nx
import numpy as np
import pytest

from polyform_algebra.matrix import Laplacian, MatrixTolerance
from polyform_algebra.matrix.networkx import from_networkx


@pytest.mark.parametrize("graph, expected, eigenvalue", [
    (nx.path_graph(5), np.array([1, np.sqrt(3), 2, np.sqrt(3), 1]) / (4 + 2*np.sqrt(3)), np.sqrt(3)),
    (nx.star_graph(4), np.array([2, 1, 1, 1, 1]) / 6, 2),
    (nx.complete_graph(5), np.ones(5) / 5, 4),
    (nx.path_graph(2), np.ones(2) / 2, 1),
])
def test_pmg_exact_eigenpairs(graph, expected, eigenvalue):
    """Known Perron vectors solve Av=lambda*v; bipartite -lambda must not be selected."""
    lap = from_networkx(graph)
    result = lap.eigenvector_centrality()
    assert result.eigenvalue == pytest.approx(eigenvalue)
    np.testing.assert_allclose(result.scores.to_numpy(), expected)
    assert result.scores.axis is lap.rows
    assert result.scores.labels == tuple(graph.nodes)
    assert result.scores.to_numpy().sum() == pytest.approx(1)
    assert np.all(result.scores.to_numpy() > 0)
    np.testing.assert_allclose(lap.to_adjacency().to_numpy() @ result.scores.to_numpy(),
                               result.eigenvalue * result.scores.to_numpy())


@pytest.mark.parametrize("scale", [1e-100, 1, 1e100])
def test_weighted_path_and_scale_invariance(scale):
    """Weights 3 and 4 give lambda=5 and v=(3,5,4); common scaling affects lambda only."""
    graph = nx.Graph()
    graph.add_weighted_edges_from([("a", "x", 3 * scale), ("x", "b", 4 * scale)])
    lap = from_networkx(graph, tolerance=MatrixTolerance(atol=0))
    result = lap.eigenvector_centrality()
    np.testing.assert_allclose(result.scores.to_numpy(), np.array([3, 5, 4]) / 12)
    assert result.eigenvalue / scale == pytest.approx(5)
    reordered = lap.reordered(rows=["b", "a", "x"], columns=["b", "a", "x"])
    assert reordered.eigenvector_centrality().scores.to_dict() == pytest.approx(
        {"a": 3/12, "x": 5/12, "b": 4/12})
    assert reordered.eigenvector_centrality().scores.labels == ("b", "a", "x")


def test_disconnected_normalization_and_isolate():
    """Component scores each sum to one even with unequal eigenvalues; isolates get score one."""
    graph = nx.Graph()
    graph.add_nodes_from(["a", "b", "x", "y", "isolated"])
    graph.add_weighted_edges_from([("a", "x", 2), ("b", "y", 5)])
    lap = from_networkx(graph)
    with pytest.raises(ValueError, match="componentwise"):
        lap.eigenvector_centrality()
    results = lap.eigenvector_centrality(componentwise=True)
    assert [r.scores.labels for r in results] == [("a", "x"), ("b", "y"), ("isolated",)]
    assert [r.eigenvalue for r in results] == pytest.approx([2, 5, 0])
    for result in results:
        assert sum(result.scores.to_dict().values()) == pytest.approx(1)
    assert results[-1].scores.to_dict() == {"isolated": 1}


def test_singleton_empty_and_directed_boundaries():
    """Singleton normalization is unique; empty default and directed implementation are explicit errors."""
    result = from_networkx(nx.empty_graph(1)).eigenvector_centrality()
    assert result.eigenvalue == 0 and result.scores.to_dict() == {0: 1}
    empty = from_networkx(nx.Graph())
    with pytest.raises(ValueError):
        empty.eigenvector_centrality()
    assert empty.eigenvector_centrality(componentwise=True) == ()
    with pytest.raises(TypeError):
        empty.eigenvector_centrality(componentwise="yes")
    with pytest.raises(NotImplementedError):
        from_networkx(nx.DiGraph([(0, 1), (1, 2), (2, 0)])).eigenvector_centrality()


def test_unresolved_perron_gap():
    """Two equal edges joined weakly have a small top gap; do not choose an arbitrary mixture."""
    graph = nx.Graph()
    graph.add_weighted_edges_from([(0, 1, 1), (2, 3, 1), (1, 2, 1e-6)])
    lap = from_networkx(graph)
    loose = Laplacian(lap.to_numpy(), rows=lap.rows, tolerance=MatrixTolerance(atol=1e-4))
    with pytest.raises(ValueError, match="gap"):
        loose.eigenvector_centrality()
    # Inversion symmetry gives equal pairs; the default tolerance resolves the gap.
    scores = lap.eigenvector_centrality().scores.to_numpy()
    np.testing.assert_allclose(scores, scores[::-1], rtol=1e-8)
    assert scores.sum() == pytest.approx(1)


def test_tuple_labels_and_pandas_column():
    """The existing Series adapter adds an independent named column without losing tuple labels."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_series
    graph = nx.Graph()
    graph.add_edge(("a", 1), ("b", 2))
    lap = from_networkx(graph)
    result = lap.eigenvector_centrality()
    series = to_series(result.scores)
    assert series.name == "eigenvector_centrality"
    assert not isinstance(series.index, pd.MultiIndex)
    frame = pd.concat([to_series(lap.spanning_tree_centrality().scores), series], axis=1)
    assert list(frame.columns) == ["spanning_tree_centrality", "eigenvector_centrality"]
    np.testing.assert_allclose(frame.to_numpy(), .5)
    series.iloc[0] = 99
    assert result.scores[("a", 1)] == pytest.approx(.5)
