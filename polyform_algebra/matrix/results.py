"""Immutable labeled numerical results; pandas stays an optional adapter."""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from .context import LabeledMatrix, MatrixAxis
from ._numeric import MatrixTolerance


@dataclass(frozen=True, eq=False, init=False)
class VertexValues:
    """Scalar values in an explicit vertex order, not a new formal algebra type.

    Indexing uses vertex labels (including integer labels), never positions.
    NumPy/dict exports are detached. The axis detects legacy Basis mutation.
    """

    axis: MatrixAxis
    name: str
    _values: tuple = field(repr=False)

    def __init__(self, values, *, axis, name):
        if not isinstance(axis, MatrixAxis):
            raise TypeError("axis must be a MatrixAxis")
        if not isinstance(name, str) or not name:
            raise ValueError("name must be a nonempty string")
        raw = np.asarray(values)
        if raw.dtype.kind not in "iuf":
            raise TypeError("vertex values must be real numbers")
        array = np.asarray(raw, dtype=float)
        if array.shape != (len(axis.labels),) or not np.isfinite(array).all():
            raise ValueError("vertex values must be finite and match the axis")
        object.__setattr__(self, "axis", axis)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "_values", tuple(array.tolist()))

    @property
    def labels(self):
        return self.axis.labels

    def __getitem__(self, label):
        try:
            position = self.labels.index(label)
        except ValueError:
            raise KeyError(label) from None
        return self._values[position]

    def to_numpy(self):
        self.labels
        return np.array(self._values, dtype=float)

    def to_dict(self):
        return dict(zip(self.labels, self._values))


@dataclass(frozen=True, eq=False)
class CircumsphereResult:
    """Resistance radius (squared linear radius) and affine center coordinates."""

    resistance_radius: float
    center_coordinates: VertexValues


@dataclass(frozen=True, eq=False)
class SpanningTreeCentralityResult:
    """Expected tree degrees and scores normalized within one component.

    For a singleton, expected degree is zero and scores is None: choosing a
    random tree edge end is undefined. No artificial unit score is assigned.
    """

    expected_degrees: VertexValues
    scores: Optional[VertexValues]


@dataclass(frozen=True, eq=False)
class EigenvectorCentralityResult:
    """Largest adjacency eigenvalue and positive vertex scores summing to one."""

    eigenvalue: float
    scores: VertexValues


@dataclass(frozen=True, eq=False)
class AnnihilatorsResult:
    """Numerical left/right kernel bases and unique sum-n representatives.

    Basis vectors are columns: left_basis.T @ L=0 and L @ right_basis=0.
    They are Euclidean-coordinate orthonormal, not orthonormal in a PMG metric.
    left/right are None unless the corresponding numerical kernel is 1-D.
    """

    left_basis: LabeledMatrix
    right_basis: LabeledMatrix
    left: Optional[VertexValues]
    right: Optional[VertexValues]
    tolerance: MatrixTolerance

    @property
    def left_dimension(self):
        return len(self.left_basis.columns.labels)

    @property
    def right_dimension(self):
        return len(self.right_basis.columns.labels)
