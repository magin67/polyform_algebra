"""Storage safeguards; these tests do not assert new PMG identities."""

import numpy as np
import pytest

from polyform_algebra import Basis
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def axis(labels=("a", "b"), space="vertex-values"):
    return MatrixAxis(space=space, basis=Basis(labels))


def test_context_preserves_order_and_zero_rows():
    """Storage keeps isolated labels, rather than deriving support from entries."""
    vertices = axis(("b", "a", "isolated"))
    matrix = LabeledMatrix(np.zeros((3, 3)), rows=vertices, columns=vertices,
                           kind=MatrixKind.FORM_COEFFICIENTS)
    assert matrix.rows.labels == ("b", "a", "isolated")
    assert matrix.to_numpy().shape == (3, 3)


def test_rectangular_matrix_and_detached_buffers():
    """Editing input/output buffers must not silently change stored data."""
    values = np.array([[1., 2.]])
    matrix = LabeledMatrix(values, rows=axis(("out",), "output"),
                           columns=axis(), kind=MatrixKind.OPERATOR)
    values[:] = 0
    exported = matrix.to_numpy()
    exported[:] = 5
    np.testing.assert_array_equal(matrix.to_numpy(), [[1, 2]])


@pytest.mark.parametrize("values, error", [
    ([[1, 2]], ValueError),
    ([[1, np.nan], [0, 1]], ValueError),
    ([[1, np.inf], [0, 1]], ValueError),
    ([[1j, 0], [0, 1]], TypeError),
    ([["1", "0"], ["0", "1"]], TypeError),
])
def test_invalid_storage(values, error):
    """Invalid storage fails before a future numerical algorithm can consume it."""
    vertices = axis()
    with pytest.raises(error):
        LabeledMatrix(values, rows=vertices, columns=vertices, kind=MatrixKind.METRIC)


def test_axes_require_basis_unique_labels_and_named_space():
    """These are context requirements, not claims of linear independence."""
    with pytest.raises(TypeError):
        MatrixAxis(space="V", basis=["a"])
    with pytest.raises(ValueError):
        axis(("a", "a"))
    with pytest.raises(ValueError):
        axis(space=" ")


@pytest.mark.parametrize("mutation", ["add", "remove", "reorder"])
def test_mutable_basis_invalidates_axis(mutation):
    """Direct mutations of the legacy list must invalidate numerical context."""
    vertices = axis()
    matrix = LabeledMatrix(np.eye(2), rows=vertices, columns=vertices,
                           kind=MatrixKind.METRIC)
    if mutation == "add":
        vertices.basis.elements.append("c")
    elif mutation == "remove":
        vertices.basis.elements.remove("a")
    else:
        vertices.basis.elements.reverse()
    with pytest.raises(ValueError, match="basis changed"):
        matrix.to_numpy()


def test_empty_axes_and_explicit_kind():
    """Zero-dimensional storage is allowed; roles cannot be guessed from values."""
    empty = axis(())
    assert LabeledMatrix(np.empty((0, 0)), rows=empty, columns=empty,
                         kind=MatrixKind.OPERATOR).to_numpy().shape == (0, 0)
    with pytest.raises(TypeError):
        LabeledMatrix(np.empty((0, 0)), rows=empty, columns=empty, kind="operator")


def test_pandas_roundtrip_tuple_labels_and_context():
    """Optional tables retain tuple IDs and never guess spaces from labels."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import from_dataframe, to_dataframe

    vertices = axis((("a", 1), ("b", 2)))
    matrix = LabeledMatrix(np.eye(2), rows=vertices, columns=vertices,
                           kind=MatrixKind.METRIC)
    frame = to_dataframe(matrix)
    assert not isinstance(frame.index, pd.MultiIndex)
    restored = from_dataframe(frame, rows=vertices, columns=vertices,
                              kind=matrix.kind)
    assert restored.rows is vertices
    np.testing.assert_array_equal(restored.to_numpy(), matrix.to_numpy())
    frame.iloc[0, 0] = 7
    assert restored.to_numpy()[0, 0] == 1


def test_pandas_refuses_implicit_alignment_and_duplicates():
    """Label alignment must be explicit to prevent plausible but wrong coordinates."""
    pd = pytest.importorskip("pandas")
    from polyform_algebra.matrix.pandas import from_dataframe

    vertices = axis()
    for labels in [("b", "a"), ("a", "a"), ("a", "c")]:
        frame = pd.DataFrame(np.eye(2), index=labels, columns=("a", "b"))
        with pytest.raises(ValueError):
            from_dataframe(frame, rows=vertices, columns=vertices,
                           kind=MatrixKind.METRIC)
