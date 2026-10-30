"""Circumspheres and weighted spanning-tree centrality from PMG v5.

Source revision: d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Centroid,%20Spanning-Tree,%20and%20Eigenvector%20Centralities.md
Equations (4)-(5), the complete graph, P5 and K1,4 examples.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Metric%20parameters/Resistance%20Radius%20and%20Diameter%20of%20a%20Graph.md
Equations (3)-(4); singleton and componentwise centrality are explicit API policies.
"""

from dataclasses import FrozenInstanceError
from fractions import Fraction
from itertools import combinations

import networkx as nx
import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis, MatrixTolerance, VertexValues
from polyform_algebra.matrix.networkx import from_networkx


def weighted_triangle():
    graph = nx.Graph()
    graph.add_weighted_edges_from([("a", "b", 2), ("b", "c", 3), ("a", "c", 4)])
    return graph


@pytest.mark.parametrize("graph, degrees, coordinates", [
    (nx.path_graph(5), [1, 2, 2, 2, 1], [.5, 0, 0, 0, .5]),
    (nx.star_graph(4), [4, 1, 1, 1, 1], [-1, .5, .5, .5, .5]),
])
def test_path_and_star(graph, degrees, coordinates):
    """Both unit trees have rho=4/4=1, but distinct degrees and center coordinates."""
    lap = from_networkx(graph)
    sphere = lap.circumsphere()
    centrality = lap.spanning_tree_centrality()
    assert sphere.resistance_radius == pytest.approx(1)
    np.testing.assert_allclose(sphere.center_coordinates.to_numpy(), coordinates, atol=1e-14)
    np.testing.assert_allclose(centrality.expected_degrees.to_numpy(), degrees)
    np.testing.assert_allclose(centrality.scores.to_numpy(), np.array(degrees) / 8)
    assert sphere.center_coordinates.axis is lap.rows
    assert sphere.center_coordinates.labels == tuple(graph.nodes)
    # Explicit Euclidean realizations with orthogonal unit edge vectors have
    # center (1/2,1/2,1/2,1/2) and squared distance one to every vertex.
    points = np.vstack([np.zeros(4), np.eye(4)])
    if nx.is_isomorphic(graph, nx.path_graph(5)):
        points = np.vstack([np.zeros(4), np.cumsum(np.eye(4), axis=0)])
    center = sphere.center_coordinates.to_numpy() @ points
    np.testing.assert_allclose(center, np.full(4, .5))
    np.testing.assert_allclose(np.sum((points - center) ** 2, axis=1), sphere.resistance_radius)


@pytest.mark.parametrize("n", [2, 3, 5])
def test_complete_graph(n):
    """K_n has R_ij=2/n, uniform center and centrality, and rho=(n-1)/n^2."""
    lap = from_networkx(nx.complete_graph(n))
    sphere = lap.circumsphere()
    centrality = lap.spanning_tree_centrality()
    assert sphere.resistance_radius == pytest.approx((n - 1) / n ** 2)
    np.testing.assert_allclose(sphere.center_coordinates.to_numpy(), np.full(n, 1 / n))
    np.testing.assert_allclose(centrality.expected_degrees.to_numpy(), np.full(n, 2 * (n - 1) / n))
    np.testing.assert_allclose(centrality.scores.to_numpy(), np.full(n, 1 / n))


def test_weighted_triangle_exact_values():
    """Tree weight sum 2*3+2*4+3*4=26 gives delta=(17,16,19)/13 and rho=105/1352."""
    lap = from_networkx(weighted_triangle())
    sphere = lap.circumsphere()
    result = lap.spanning_tree_centrality()
    np.testing.assert_allclose(result.expected_degrees.to_numpy(), np.array([17, 16, 19]) / 13)
    np.testing.assert_allclose(sphere.center_coordinates.to_numpy(), np.array([9, 10, 7]) / 26)
    # R_ab=7/26, R_ac=5/26, R_bc=6/26; R*alpha=105/676 * 1.
    assert sphere.resistance_radius == pytest.approx(Fraction(105, 1352))


@pytest.mark.parametrize("graph", [weighted_triangle(), nx.cycle_graph(4)])
def test_expected_degrees_by_tree_enumeration(graph):
    """Independent enumeration averages ordinary degrees with product-of-conductance weights."""
    labels = list(graph.nodes)
    totals = dict.fromkeys(labels, Fraction(0))
    total_weight = Fraction(0)
    for edges in combinations(graph.edges(data=True), len(labels) - 1):
        tree = nx.Graph()
        tree.add_nodes_from(labels)
        tree.add_edges_from((u, v) for u, v, _ in edges)
        if not nx.is_tree(tree):
            continue
        weight = Fraction(1)
        for _, _, data in edges:
            weight *= Fraction(data.get("weight", 1))
        total_weight += weight
        for label in labels:
            totals[label] += weight * tree.degree(label)
    result = from_networkx(graph).spanning_tree_centrality()
    expected = [float(totals[label] / total_weight) for label in labels]
    np.testing.assert_allclose(result.expected_degrees.to_numpy(), expected)
    assert sum(result.expected_degrees.to_dict().values()) == pytest.approx(2 * (len(labels) - 1))
    assert sum(result.scores.to_dict().values()) == pytest.approx(1)


def test_center_by_independent_resistance_solve():
    """A grounded inverse gives R independently; solving Rv=1 recovers alpha and rho."""
    lap = from_networkx(weighted_triangle())
    a = lap.to_numpy()
    grounded = np.linalg.solve(a[1:, 1:], np.eye(2))
    gram = np.zeros_like(a)
    gram[1:, 1:] = grounded
    diagonal = np.diag(gram)
    resistance = diagonal[:, None] + diagonal[None, :] - 2 * gram
    v = np.linalg.solve(resistance, np.ones(3))
    sphere = lap.circumsphere()
    alpha = sphere.center_coordinates.to_numpy()
    np.testing.assert_allclose(alpha, v / v.sum())
    assert sphere.resistance_radius == pytest.approx(1 / (2 * v.sum()))
    np.testing.assert_allclose(resistance @ alpha, np.full(3, 2 * sphere.resistance_radius))


def test_scaling_and_reordering():
    """Scaling conductances by t divides rho by t but preserves alpha and normalized centrality."""
    lap = from_networkx(weighted_triangle())
    scaled = Laplacian(7 * lap.to_numpy(), rows=lap.rows)
    assert scaled.circumsphere().resistance_radius == pytest.approx(lap.circumsphere().resistance_radius / 7)
    np.testing.assert_allclose(scaled.circumsphere().center_coordinates.to_numpy(),
                               lap.circumsphere().center_coordinates.to_numpy())
    np.testing.assert_allclose(scaled.spanning_tree_centrality().scores.to_numpy(),
                               lap.spanning_tree_centrality().scores.to_numpy())
    order = ["c", "a", "b"]
    reordered = lap.reordered(rows=order, columns=order)
    assert reordered.circumsphere().center_coordinates.labels == tuple(order)
    assert reordered.circumsphere().center_coordinates.to_dict() == pytest.approx(
        lap.circumsphere().center_coordinates.to_dict())


def test_disconnected_components_and_isolate():
    """Two edges have separate spheres/scores; an isolate has rho=0, alpha=1, delta=0, scores=None."""
    graph = nx.Graph()
    graph.add_nodes_from(["a", "b", "x", "y", "o"])
    graph.add_weighted_edges_from([("a", "x", 2), ("b", "y", 4)])
    lap = from_networkx(graph)
    for method in (lap.circumsphere, lap.spanning_tree_centrality):
        with pytest.raises(ValueError, match="componentwise"):
            method()
    spheres = lap.circumsphere(componentwise=True)
    centralities = lap.spanning_tree_centrality(componentwise=True)
    assert [r.center_coordinates.labels for r in spheres] == [("a", "x"), ("b", "y"), ("o",)]
    # A single edge of conductance c has endpoint resistance 1/c and rho=1/(4c).
    assert [r.resistance_radius for r in spheres] == pytest.approx([1/8, 1/16, 0])
    for result in centralities[:2]:
        np.testing.assert_allclose(result.scores.to_numpy(), [.5, .5])
    assert centralities[2].scores is None
    assert centralities[2].expected_degrees.to_dict() == {"o": 0}
    assert spheres[2].center_coordinates.to_dict() == {"o": 1}


def test_singleton_empty_and_explicit_component_mode():
    """Singleton geometry is a point; empty default is undefined and empty component results are ()."""
    lap = from_networkx(nx.empty_graph(1))
    assert lap.circumsphere().resistance_radius == 0
    assert lap.circumsphere().center_coordinates[0] == 1
    assert lap.spanning_tree_centrality().scores is None
    assert len(lap.circumsphere(componentwise=True)) == 1
    empty = from_networkx(nx.Graph())
    for method in (empty.circumsphere, empty.spanning_tree_centrality):
        with pytest.raises(ValueError, match="nonempty"):
            method()
        assert method(componentwise=True) == ()
        with pytest.raises(TypeError):
            method(componentwise="yes")


def test_directed_and_numerically_unresolved_domains():
    """Directed geometry is deferred; weak unresolved components are never silently split."""
    graph = nx.DiGraph([(0, 1), (1, 2), (2, 0)])
    for method in (from_networkx(graph).circumsphere, from_networkx(graph).spanning_tree_centrality):
        with pytest.raises(NotImplementedError):
            method()
    axis = MatrixAxis(basis=Basis(["a", "b"], space=Space("weak edge", 2)))
    data = 1e-15 * np.array([[1, -1], [-1, 1]])
    lap = Laplacian(data, rows=axis)
    with pytest.raises(ValueError, match="tolerance"):
        lap.circumsphere()
    precise = Laplacian(data, rows=axis, tolerance=MatrixTolerance(atol=0))
    assert precise.circumsphere().resistance_radius == pytest.approx(1 / (4e-15))


def test_labeled_values_and_detached_pandas_export():
    """Tuple/integer labels stay scalar identities; exported arrays and Series cannot mutate results."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_series
    graph = nx.Graph()
    graph.add_edge(("a", 1), ("b", 2), weight=2)
    lap = from_networkx(graph)
    coords = lap.circumsphere().center_coordinates
    series = to_series(coords)
    assert not isinstance(series.index, pd.MultiIndex)
    assert series.name == "center_coordinate"
    assert coords[("a", 1)] == pytest.approx(.5)
    series.iloc[0] = 99
    array = coords.to_numpy()
    array[0] = 99
    assert coords[("a", 1)] == pytest.approx(.5)
    with pytest.raises(FrozenInstanceError):
        coords.name = "other"
    with pytest.raises(KeyError):
        coords[0]
    with pytest.raises(TypeError):
        to_series(None)
    lap.rows.basis.elements.reverse()
    with pytest.raises(ValueError, match="basis changed"):
        coords.to_dict()


@pytest.mark.parametrize("values", [[np.nan], [[1]], [1, 2]])
def test_vertex_values_reject_bad_shape_or_nonfinite(values):
    """Named numerical results require one finite scalar per vertex, without broadcasting."""
    axis = MatrixAxis(basis=Basis(["o"], space=Space("point", 1)))
    with pytest.raises(ValueError):
        VertexValues(values, axis=axis, name="test")
