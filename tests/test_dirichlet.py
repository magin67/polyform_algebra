"""Undirected Dirichlet projections and Laplacian compatibility.

Source: PMG v5, revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Selected%20Vertices%20and%20the%20Dirichlet%20Problem.md
Weighted path example, equations (3), (7), (8), and 'Disconnected graphs'.
"""

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import (
    LabeledMatrix, Laplacian, MatrixAxis, MatrixKind, MatrixTolerance,
)


def path_laplacian():
    axis = MatrixAxis(basis=Basis(["a", "x", "y", "b"], space=Space("path", 4)))
    # Conductances 2, 1, 3 on a-x-y-b give these weighted degrees.
    return Laplacian([[2, -2, 0, 0], [-2, 3, -1, 0],
                      [0, -1, 4, -3], [0, 0, -3, 3]], rows=axis)


def test_weighted_path_projection_coordinates_and_extension():
    """Solving [[3,-1],[-1,4]] B=diag(2,3) gives [[8,3],[2,9]]/11."""
    lap = path_laplacian()
    result = lap.dirichlet(["a", "b"])
    assert result.selected == ("a", "b") and result.internal == ("x", "y")
    assert result.source_axis is lap.rows
    b = result.influence_matrix
    assert b.kind == MatrixKind.OPERATOR and b.representation == "generic"
    np.testing.assert_allclose(b.to_numpy(), np.array([[8, 3], [2, 9]]) / 11)
    assert result.extend({"b": 11, "a": 0}) == pytest.approx({"a": 0, "x": 3, "y": 9, "b": 11})
    assert list(result.extend([0, 11])) == ["a", "x", "y", "b"]


def test_projection_residual_is_orthogonal():
    """In resistance geometry, x-q_x and y-q_y are orthogonal to b-a."""
    lap = path_laplacian()
    b = lap.dirichlet(["a", "b"]).influence_matrix.to_numpy()
    green = lap.to_greenian().to_numpy()
    direction = np.array([-1, 0, 0, 1])
    for row, i in enumerate([1, 2]):
        residual = np.eye(4)[i]
        residual[[0, 3]] -= b[row]
        assert residual @ green @ direction == pytest.approx(0, abs=1e-14)
        # The note independently obtains residual norms 4/11 and 3/11.
        assert residual @ green @ residual == pytest.approx([4/11, 3/11][row])


def test_harmonicity_maximum_principle_and_energy():
    """The harmonic extension minimizes energy: perturbations add h.T L h > 0."""
    lap = path_laplacian()
    result = lap.dirichlet(["a", "b"])
    b = result.influence_matrix.to_numpy()
    assert np.all(b >= 0)
    np.testing.assert_allclose(b.sum(axis=1), 1)
    f = np.array(list(result.extend([0, 11]).values()))
    a = lap.to_numpy()
    np.testing.assert_allclose((a @ f)[1:3], 0, atol=1e-14)
    assert np.all((f >= 0) & (f <= 11))
    # 2*3^2 + (9-3)^2 + 3*(11-9)^2 = 66.
    assert f @ a @ f == pytest.approx(66)
    h = np.array([0, 2, -1, 0])
    assert (f + h) @ a @ (f + h) == pytest.approx(66 + h @ a @ h)
    assert h @ a @ h > 0


def test_selection_order_and_source_permutation():
    """Selected order controls columns; source order controls internal rows and extension."""
    lap = path_laplacian().reordered(rows=["y", "b", "x", "a"], columns=["y", "b", "x", "a"])
    assert isinstance(lap, Laplacian)
    result = lap.dirichlet(["b", "a"])
    assert result.internal == ("y", "x")
    np.testing.assert_allclose(result.influence_matrix.to_numpy(), np.array([[9, 2], [3, 8]]) / 11)
    assert list(result.extend([11, 0])) == ["y", "b", "x", "a"]


def test_one_or_all_selected():
    """One selected vertex yields constant extension; all selected leaves zero internal rows."""
    lap = path_laplacian()
    one = lap.dirichlet(["x"])
    np.testing.assert_allclose(one.influence_matrix.to_numpy(), np.ones((3, 1)))
    assert one.extend([7]) == pytest.approx(dict.fromkeys(lap.rows.labels, 7))
    all_selected = lap.dirichlet(["b", "y", "x", "a"])
    assert all_selected.influence_matrix.to_numpy().shape == (0, 4)
    assert all_selected.extend([11, 9, 3, 0]) == {"a": 0, "x": 3, "y": 9, "b": 11}


def test_components_and_isolated_vertex():
    """Every component needs an anchor; cross-component influence is exactly zero."""
    axis = MatrixAxis(basis=Basis(range(5), space=Space("disconnected", 5)))
    a = np.zeros((5, 5))
    a[:2, :2] = [[2, -2], [-2, 2]]
    a[2:4, 2:4] = [[3, -3], [-3, 3]]
    lap = Laplacian(a, rows=axis)
    for selected in ([0], [0, 2]):
        with pytest.raises(ValueError, match="not unique"):
            lap.dirichlet(selected)
    result = lap.dirichlet([4, 2, 0])
    np.testing.assert_array_equal(result.influence_matrix.to_numpy(), [[0, 0, 1], [0, 1, 0]])
    assert result.extend([8, 5, 3]) == {0: 3, 1: 3, 2: 5, 3: 5, 4: 8}


def test_singleton():
    """An isolated selected point has no unknowns and preserves its assigned value."""
    axis = MatrixAxis(basis=Basis(["o"], space=Space("point", 1)))
    result = Laplacian([[0]], rows=axis).dirichlet(["o"])
    assert result.influence_matrix.to_numpy().shape == (0, 1)
    assert result.extend({"o": 4}) == {"o": 4}


@pytest.mark.parametrize("selected", [[], ["a", "a"], ["unknown"]])
def test_invalid_selection(selected):
    """API selection must be a nonempty unique subset, not missing or repeated labels."""
    with pytest.raises(ValueError):
        path_laplacian().dirichlet(selected)


@pytest.mark.parametrize("selected", ["a", {"a", "b"}])
def test_unordered_or_bare_string_selection(selected):
    """An explicit ordered iterable prevents accidental character or set ordering."""
    with pytest.raises(TypeError):
        path_laplacian().dirichlet(selected)


@pytest.mark.parametrize("values", [{"a": 1}, {"a": 1, "b": 2, "x": 3}, [1], [[1, 2]], [np.nan, 2]])
def test_invalid_prescribed_values(values):
    """Missing labels, dimensions and nonfinite data must not silently align or propagate."""
    with pytest.raises(ValueError):
        path_laplacian().dirichlet(["a", "b"]).extend(values)


def test_directed_domain_and_unresolved_rank():
    """Directed solves are deferred; unresolved blocks cannot use a pseudoinverse fallback."""
    axis = MatrixAxis(basis=Basis(["a", "b"], space=Space("V", 2)))
    with pytest.raises(NotImplementedError):
        Laplacian([[1, -1], [0, 0]], rows=axis).dirichlet(["b"])
    values = np.array([[1, -1], [-1, 1]]) * 1e-15
    with pytest.raises(ValueError, match="tolerance"):
        Laplacian(values, rows=axis).dirichlet(["a"])
    lap = Laplacian(values, rows=axis, tolerance=MatrixTolerance(atol=0))
    np.testing.assert_allclose(lap.dirichlet(["a"]).influence_matrix.to_numpy(), [[1]])


def test_laplacian_invariants_and_legacy_conversions():
    """Subclass preserves matrix APIs; semantic conversion upgrades legacy instances."""
    lap = path_laplacian()
    assert isinstance(lap, LabeledMatrix)
    legacy = LabeledMatrix(lap.to_numpy(), rows=lap.rows, columns=lap.columns,
                            kind=MatrixKind.OPERATOR, representation="laplacian")
    assert isinstance(legacy.to_laplacian(), Laplacian)
    assert isinstance(lap.to_adjacency().to_laplacian(), Laplacian)
    assert isinstance(lap.to_resistance().to_laplacian(), Laplacian)
    for other in (lap.to_greenian(), lap.to_adjacency(), lap.pseudoinverse(),
                  lap.with_representation("generic")):
        assert not isinstance(other, Laplacian)
    with pytest.raises(ValueError):
        Laplacian(lap.to_numpy(), rows=lap.rows, representation="generic")
    with pytest.raises(ValueError):
        Laplacian(np.eye(4), rows=lap.rows)
    with pytest.raises(FrozenInstanceError):
        lap.representation = "generic"


def test_stale_source_basis_is_detected():
    """A result cannot silently remap its embedding after legacy basis-list mutation."""
    lap = path_laplacian()
    result = lap.dirichlet(["a", "b"])
    lap.rows.basis.elements.reverse()
    with pytest.raises(ValueError, match="basis changed"):
        result.extend([0, 11])


def test_pandas_and_networkx_adapters():
    """Adapters produce Laplacian and preserve labels when exporting projection coordinates."""
    pd = pytest.importorskip("pandas")
    import networkx as nx
    from polyform_algebra.matrix.networkx import from_networkx
    from polyform_algebra.matrix.pandas import from_dataframe, to_dataframe
    graph = nx.Graph()
    graph.add_weighted_edges_from([("a", "x", 2), ("x", "y", 1), ("y", "b", 3)])
    lap = from_networkx(graph)
    assert isinstance(lap, Laplacian)
    imported = from_dataframe(to_dataframe(lap), rows=lap.rows, columns=lap.columns,
                              kind=lap.kind, representation="laplacian")
    assert isinstance(imported, Laplacian)
    result = imported.dirichlet(["a", "b"])
    frame = to_dataframe(result.influence_matrix)
    assert list(frame.index) == ["x", "y"] and list(frame.columns) == ["a", "b"]
    assert frame.loc["x", "b"] == pytest.approx(3/11)
    with pytest.raises(TypeError):
        result.extend(pd.Series({"b": 11, "a": 0}))
    assert result.extend(pd.Series({"b": 11, "a": 0}).to_dict())["x"] == pytest.approx(3)
