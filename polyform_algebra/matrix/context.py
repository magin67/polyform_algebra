"""Numerical storage and explicit axis context; no implicit tensor operations."""

from dataclasses import dataclass, field
from enum import Enum

import numpy as np

from .._core.basis import Basis


class MatrixKind(Enum):
    """Mathematical roles, not graph representation names or validity claims."""

    OPERATOR = "operator"
    METRIC = "metric"
    FORM_COEFFICIENTS = "form_coefficients"
    PAIRWISE = "pairwise"


@dataclass(frozen=True, eq=False)
class MatrixAxis:
    """An ordered Basis in an explicitly named space.

    The basis is reused, not reimplemented. Its ordered elements are captured
    at construction; later basis mutations invalidate this context. Space
    names are caller-controlled identifiers, not inferred from element labels.
    An axis does not certify linear independence of arbitrary basis elements.
    """

    space: str
    basis: Basis
    _labels: tuple = field(init=False, repr=False)

    def __post_init__(self):
        if not isinstance(self.space, str) or not self.space.strip():
            raise ValueError("space must be a nonempty string")
        if not isinstance(self.basis, Basis):
            raise TypeError("basis must be a Basis")
        labels = tuple(self.basis.elements)
        if len(set(labels)) != len(labels):
            raise ValueError("axis basis elements must be distinct")
        object.__setattr__(self, "_labels", labels)

    @property
    def labels(self):
        """Return the captured ordering, rejecting stale mutable bases."""
        if tuple(self.basis.elements) != self._labels:
            raise ValueError("axis basis changed; construct a new matrix context")
        return self._labels


@dataclass(frozen=True, eq=False, init=False)
class LabeledMatrix:
    """A finite real dense array with row/column contexts and an explicit role.

    Construction validates storage, not properties such as symmetry or the
    Laplacian sign convention. Specialized constructors will provide those
    guarantees. No arithmetic or automatic label alignment is implied.
    """

    rows: MatrixAxis
    columns: MatrixAxis
    kind: MatrixKind
    _values: np.ndarray = field(repr=False)

    def __init__(self, values, *, rows, columns, kind):
        if not isinstance(rows, MatrixAxis) or not isinstance(columns, MatrixAxis):
            raise TypeError("rows and columns must be MatrixAxis instances")
        if not isinstance(kind, MatrixKind):
            raise TypeError("kind must be a MatrixKind")
        array = np.asarray(values)
        if array.dtype.kind not in "iuf":
            raise TypeError("matrix entries must be real numeric values")
        if array.shape != (len(rows.labels), len(columns.labels)):
            raise ValueError("matrix shape must match both axis contexts")
        if not np.isfinite(array).all():
            raise ValueError("matrix entries must be finite")
        # Copy caller-owned buffers so edits cannot silently change a matrix.
        array = np.array(array, dtype=float, copy=True)
        if not np.isfinite(array).all():
            raise ValueError("matrix entries must fit finite float storage")
        array.setflags(write=False)
        object.__setattr__(self, "rows", rows)
        object.__setattr__(self, "columns", columns)
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "_values", array)

    def to_numpy(self):
        """Return an independent array in the declared axis order."""
        self.rows.labels
        self.columns.labels
        return self._values.copy()
