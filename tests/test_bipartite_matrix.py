"""Conductance conversions and Kron reduction of a participation table.

PMG v5 revision f4e80ba6d7bec8d3e32774b957fe921bf76de264:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Couplings%20Within%20One%20Part%20of%20a%20Bipartite%20Graph.md
Equations (3), (6), 'Two induced graphs', and 'Calculation ... limitations'.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Introduction/Laplacian%20-%20graph,%20electrical%20network%20and%20quadratic%20form.md
Section 'The Laplacian matrix': L=diag(A1)-A.
"""

import numpy as np
import pytest
from polyform_algebra import Basis, Space
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def table(values, representation="bipartite"):
    values = np.asarray(values, dtype=float)
    rows = MatrixAxis(basis=Basis(range(values.shape[0]), space=Space("objects", values.shape[0])))
    columns = (rows if representation == "adjacency" else
               MatrixAxis(basis=Basis(range(values.shape[1]), space=Space("objects", values.shape[1]))))
    return LabeledMatrix(values, rows=rows, columns=columns,
                         kind=MatrixKind.PAIRWISE, representation=representation)


def test_pmg_employee_example():
    """Equation (3) gives ab=bc=1/2+1/3 and ac=1/3; (6) fixes the diagonal."""
    h = table([[1, 1, 0, 0], [1, 1, 1, 0], [0, 1, 1, 1]])
    lap = h.to_unipartite(h.rows.space)
    expected = np.array([[7, -5, -2], [-5, 10, -5], [-2, -5, 7]]) / 6
    np.testing.assert_allclose(lap.to_numpy(), expected)
    assert lap.rows is h.rows and lap.columns is h.rows
    assert lap.vector_type() == "vector"
    # Triangle formula R_ab=(c_ac+c_bc)/(c_ab*c_ac+c_ab*c_bc+c_ac*c_bc).
    np.testing.assert_allclose(lap.to_resistance().to_numpy(),
                               [[0, 14/15, 4/3], [14/15, 0, 14/15], [4/3, 14/15, 0]], atol=1e-14)
    adjacency = h.to_unipartite(h.rows.space, representation="adjacency")
    np.testing.assert_allclose(adjacency.to_laplacian().to_numpy(), expected)
    np.testing.assert_array_equal(h.to_numpy().sum(axis=1), [2, 3, 3])


@pytest.mark.parametrize("side", ["rows", "columns"])
def test_schur_complement_and_full_graph_resistance(side):
    """Elimination equals D-H E^-1 H.T and preserves full-network resistances."""
    h = table([[1, 2, 0], [3, 0, 4]])
    values = h.to_numpy() if side == "rows" else h.to_numpy().T
    selected = getattr(h, side)
    lap = h.to_unipartite(selected.space)
    expected = np.diag(values.sum(axis=1)) - values @ np.linalg.solve(
        np.diag(values.sum(axis=0)), values.T)
    np.testing.assert_allclose(lap.to_numpy(), expected)
    n, m = values.shape
    full = table(np.block([[np.zeros((n, n)), values],
                           [values.T, np.zeros((m, m))]]), "adjacency")
    np.testing.assert_allclose(lap.to_resistance().to_numpy(),
                               full.to_resistance().to_numpy()[:n, :n], atol=1e-13)


def test_two_parts_have_different_couplings():
    """The note's weighted example yields 1*2/3 for employees and 1*1/2 for documents."""
    h = table([[1, 1], [2, 0]])
    assert h.to_unipartite(h.rows.space).to_adjacency().to_numpy()[0, 1] == pytest.approx(2/3)
    assert h.to_unipartite(h.columns.space).to_adjacency().to_numpy()[0, 1] == pytest.approx(1/2)


def test_isolates_zero_columns_and_disconnected_components():
    """A zero column contributes nothing; lone participants and empty rows stay isolates."""
    h = table([[2, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 0]])
    lap = h.to_unipartite(h.rows.space)
    np.testing.assert_array_equal(lap.to_numpy(), [[1, -1, 0, 0], [-1, 1, 0, 0],
                                                 [0, 0, 0, 0], [0, 0, 0, 0]])
    assert lap.components == ((0, 1), (2,), (3,))
    assert np.isinf(lap.to_resistance().to_numpy()[0, 3])


@pytest.mark.parametrize("shape", [(0, 0), (0, 3), (3, 0), (1, 4), (3, 2)])
def test_zero_tables_and_empty_parts(shape):
    """No original couplings means a zero Laplacian on every retained vertex."""
    h = table(np.zeros(shape))
    for axis in (h.rows, h.columns):
        np.testing.assert_array_equal(h.to_unipartite(axis.space).to_numpy(),
                                      np.zeros((axis.space.dimension,) * 2))


def test_scaling_and_permutation():
    """Each coefficient is homogeneous of degree one; labels determine ordering only."""
    h = table([[1, 2], [3, 4], [0, 5]])
    result = h.to_unipartite(h.rows.space).to_numpy()
    scaled = table(7 * h.to_numpy())
    np.testing.assert_allclose(scaled.to_unipartite(scaled.rows.space).to_numpy(), 7 * result)
    reordered = h.reordered(rows=[2, 0, 1], columns=[1, 0])
    np.testing.assert_allclose(reordered.to_unipartite(h.rows.space).to_numpy(),
                               result[np.ix_([2, 0, 1], [2, 0, 1])])


@pytest.mark.parametrize("weight", [1e-100, 1e308])
def test_finite_extreme_conductances(weight):
    """Two equal conductances in series yield c/2, without degree/product overflow."""
    h = table([[weight], [weight]])
    np.testing.assert_allclose(h.to_unipartite(h.rows.space).to_adjacency().to_numpy() / weight,
                               [[0, .5], [.5, 0]])


@pytest.mark.parametrize("values", [[[0, 1, 0], [1, 0, 2], [0, 2, 0]],
                                    [[0, 3, 0], [0, 0, 0], [0, 0, 0]], np.zeros((0, 0))])
def test_adjacency_roundtrip(values):
    """Outgoing L=diag(A1)-A is inverted off the diagonal, including directed arcs."""
    a = table(values, "adjacency")
    lap = a.to_laplacian()
    np.testing.assert_array_equal(lap.to_numpy(), np.diag(a.to_numpy().sum(axis=1)) - a.to_numpy())
    np.testing.assert_array_equal(lap.to_adjacency().to_numpy(), a.to_numpy())
    if a.is_symmetric():
        np.testing.assert_allclose(a.to_greenian().to_adjacency().to_numpy(), a.to_numpy(), atol=1e-13)
    else:
        with pytest.raises(NotImplementedError):
            a.to_resistance()


@pytest.mark.parametrize("values", [[[-1]], [[np.inf]], [[np.nan]]])
def test_invalid_conductances(values):
    """The supported electrical model requires finite nonnegative conductances."""
    with pytest.raises(ValueError):
        table(values)


def test_context_and_dispatch_restrictions():
    """Distinct parts require Space identity, explicit meaning, and a supported output."""
    h = table([[1, 2], [3, 4]])
    with pytest.raises(ValueError):
        h.to_unipartite(Space("objects", 2))
    with pytest.raises(ValueError):
        h.to_unipartite(h.rows.space, representation="greenian")
    for conversion in (h.to_greenian, h.to_laplacian, h.to_adjacency, h.to_resistance):
        with pytest.raises(ValueError):
            conversion()
    with pytest.raises(ValueError):
        LabeledMatrix(h.to_numpy(), rows=h.rows, columns=h.rows,
                      kind=MatrixKind.PAIRWISE, representation="bipartite")
    with pytest.raises(ValueError):
        table([[1, 0], [0, 0]], "adjacency")
    with pytest.raises(ValueError):
        h.with_representation("generic").to_unipartite(h.rows.space)


def test_dataframe_named_parts():
    """The pandas adapter retains employee/project labels through reduction and export."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import from_dataframe, to_dataframe
    rows = MatrixAxis(basis=Basis(["Alice", "Bob"], space=Space("employees", 2)))
    cols = MatrixAxis(basis=Basis(["Project"], space=Space("projects", 1)))
    h = from_dataframe(pd.DataFrame([[2], [6]], index=rows.labels, columns=cols.labels),
                       rows=rows, columns=cols, kind=MatrixKind.PAIRWISE, representation="bipartite")
    result = to_dataframe(h.to_unipartite(rows.space, representation="adjacency"))
    assert list(result.index) == list(result.columns) == ["Alice", "Bob"]
    assert result.loc["Alice", "Bob"] == pytest.approx(2 * 6 / 8)
