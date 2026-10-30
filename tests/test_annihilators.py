"""Directed Laplacian kernels, with independent minor and nullspace checks.

PMG v5 revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balanced%20Graphs%20and%20Spanning-Tree%20Potentials.md
'A balance criterion using spanning-tree potentials' and 'Why positivity matters'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Balanced%20graphs/Balancing%20a%20Directed%20Graph%20by%20Spanning-Tree%20Potentials.md
'A three-participant example', tau=(5,4,7).
"""

import numpy as np
import pytest

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis


def laplacian(values, labels=None):
    values = np.asarray(values, dtype=float)
    labels = list(range(len(values))) if labels is None else labels
    axis = MatrixAxis(basis=Basis(labels, space=Space("vertices", len(labels))))
    return Laplacian(values, rows=axis)


@pytest.mark.parametrize("values", [
    [[1, -1, 0], [0, 1, -1], [-1, 0, 1]],
    [[1, -1, 0], [-1, 2, -1], [0, -1, 1]],
])
def test_balanced_and_symmetric_normalized_annihilators(values):
    """A connected balanced graph has both kernels spanned by ones, even without symmetry."""
    lap = laplacian(values, ["a", "b", "c"])
    result = lap.annihilators()
    assert result.left_dimension == result.right_dimension == 1
    np.testing.assert_allclose(result.left.to_numpy(), np.ones(3))
    np.testing.assert_array_equal(result.right.to_numpy(), np.ones(3))
    assert result.left_basis.rows is lap.rows
    assert result.right_basis.rows is lap.columns
    for basis in (result.left_basis, result.right_basis):
        np.testing.assert_allclose(basis.to_numpy().T @ basis.to_numpy(), [[1]])


def test_unbalanced_example_against_principal_minors():
    """Principal minors (5,4,7) give the independent left vector 3*tau/sum(tau)."""
    lap = laplacian([[3, -1, -2], [-2, 3, -1], [-1, -1, 2]], ["a", "b", "c"])
    result = lap.annihilators()
    a = lap.to_numpy()
    minors = [np.linalg.det(np.delete(np.delete(a, i, axis=0), i, axis=1)) for i in range(3)]
    np.testing.assert_allclose(minors, [5, 4, 7])
    np.testing.assert_allclose(result.left.to_numpy(), 3 * np.array(minors) / sum(minors))
    np.testing.assert_allclose(result.left.to_numpy(), [15/16, 12/16, 21/16])
    np.testing.assert_allclose(result.left.to_numpy() @ a, 0, atol=1e-14)
    np.testing.assert_allclose(a @ result.right.to_numpy(), 0, atol=1e-14)
    assert result.left.to_numpy().sum() == pytest.approx(3)
    assert result.right.to_numpy().sum() == pytest.approx(3)
    # The normalized vector is not the raw spanning-tree potential itself.
    assert not np.allclose(result.left.to_numpy(), minors)


def test_multidimensional_kernels_have_no_preferred_vector():
    """For a->b,a->c, left=(0,s,t), right=((s+t)/2,s,t); both dimensions are two."""
    lap = laplacian([[2, -1, -1], [0, 0, 0], [0, 0, 0]])
    result = lap.annihilators()
    assert result.left_dimension == result.right_dimension == 2
    assert result.left is result.right is None
    expected_left = np.array([[0, 0], [1, 0], [0, 1.]])
    expected_right = np.array([[.5, .5], [1, 0], [0, 1.]])
    for actual, expected in ((result.left_basis, expected_left), (result.right_basis, expected_right)):
        q = actual.to_numpy()
        # Compare subspaces, not the arbitrary SVD rotation or column signs.
        projector = expected @ np.linalg.solve(expected.T @ expected, expected.T)
        np.testing.assert_allclose(q @ q.T, projector, atol=1e-14)
        np.testing.assert_allclose(q.T @ q, np.eye(2), atol=1e-14)
    np.testing.assert_allclose(result.left_basis.to_numpy().T @ lap.to_numpy(), 0, atol=1e-14)
    np.testing.assert_allclose(lap.to_numpy() @ result.right_basis.to_numpy(), 0, atol=1e-14)


def test_one_dimensional_reducible_graph_can_have_zero_left_entries():
    """For a->b, only the sink has a rooted tree; normalized left=(0,2), right=(1,1)."""
    result = laplacian([[1, -1], [0, 0]]).annihilators()
    np.testing.assert_allclose(result.left.to_numpy(), [0, 2], atol=1e-14)
    np.testing.assert_array_equal(result.right.to_numpy(), [1, 1])


@pytest.mark.parametrize("values, dimension", [
    (np.zeros((0, 0)), 0), (np.zeros((1, 1)), 1), (np.zeros((3, 3)), 3),
    ([[1, -1, 0], [-1, 1, 0], [0, 0, 0]], 2),
])
def test_empty_zero_and_disconnected_graphs(values, dimension):
    """Nullity includes isolates; zero and empty matrices retain every explicitly supplied vertex."""
    result = laplacian(values).annihilators()
    assert result.left_dimension == result.right_dimension == dimension
    assert result.left_basis.to_numpy().shape == (len(values), dimension)
    if dimension == 1:
        assert result.left.to_dict() == result.right.to_dict() == {0: 1}
    else:
        assert result.left is result.right is None


def test_tolerance_is_explicit_and_does_not_mutate_matrix():
    """Numerical nullity follows the singular cutoff; lowering atol resolves a weak edge."""
    lap = laplacian(np.array([[1, -1], [-1, 1]]) * 1e-15)
    assert lap.annihilators().right_dimension == 2
    resolved = lap.annihilators(atol=0)
    assert resolved.right_dimension == 1
    assert resolved.tolerance.atol == 0 and lap.tolerance.atol == 1e-10
    np.testing.assert_allclose(resolved.left.to_numpy(), [1, 1])


def test_reordering_and_scaling():
    """Changing vertex order permutes coordinates; a common positive scale does not change kernels."""
    lap = laplacian([[3, -1, -2], [-2, 3, -1], [-1, -1, 2]], ["a", "b", "c"])
    order = ["c", "a", "b"]
    reordered = lap.reordered(rows=order, columns=order)
    result = reordered.annihilators()
    assert result.left.labels == tuple(order)
    assert result.left.to_dict() == pytest.approx(lap.annihilators().left.to_dict())
    large = Laplacian(lap.to_numpy() * 1e200, rows=lap.rows).annihilators()
    np.testing.assert_allclose(large.left.to_numpy(), lap.annihilators().left.to_numpy())


def test_dataframe_basis_and_normalized_series():
    """Vertex rows and coefficient-space columns remain explicit through pandas adapters."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import to_dataframe, to_series
    lap = laplacian([[1, -1], [0, 0]], [("a", 1), ("b", 2)])
    result = lap.annihilators()
    frame = to_dataframe(result.left_basis)
    assert not isinstance(frame.index, pd.MultiIndex)
    assert tuple(frame.index) == lap.rows.labels
    assert list(frame.columns) == ["left_0"]
    assert to_series(result.left).name == "left_annihilator"
    frame.iloc[0, 0] = 99
    assert result.left_basis.to_numpy()[0, 0] == pytest.approx(0)
