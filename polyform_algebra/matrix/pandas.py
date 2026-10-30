"""Optional DataFrame adapters with strict ordering and no implicit alignment."""

from .context import LabeledMatrix
from .laplacian import Laplacian
from .results import VertexValues


def from_dataframe(frame, *, rows, columns, kind, representation="generic", tolerance=None):
    """Import explicit contexts; reject missing, duplicate or reordered labels.

    A DataFrame alone does not identify spaces or tensor roles. Reindex it
    explicitly before import if its order differs from the required basis.
    """
    import pandas as pd

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    if not frame.index.is_unique or not frame.columns.is_unique:
        raise ValueError("DataFrame labels must be unique")
    if tuple(frame.index) != rows.labels or tuple(frame.columns) != columns.labels:
        raise ValueError("DataFrame labels and order must match axis contexts")
    matrix_class = Laplacian if representation == "laplacian" else LabeledMatrix
    return matrix_class(frame.to_numpy(), rows=rows, columns=columns, kind=kind,
                         representation=representation, tolerance=tolerance)


def to_dataframe(matrix):
    """Export a detached table; spaces and mathematical role stay in the matrix."""
    import pandas as pd

    return pd.DataFrame(
        matrix.to_numpy(),
        index=pd.Index(matrix.rows.labels, tupleize_cols=False),
        columns=pd.Index(matrix.columns.labels, tupleize_cols=False),
    )


def to_series(values):
    """Export detached vertex values with their scalar labels and result name.

    Undefined singleton centrality has scores=None; handle that explicitly
    instead of silently inventing a value or interpreting None as an empty series.
    """
    import pandas as pd

    if not isinstance(values, VertexValues):
        raise TypeError("values must be VertexValues; undefined scores must be handled explicitly")
    return pd.Series(values.to_numpy(), index=pd.Index(values.labels, tupleize_cols=False),
                     name=values.name)
