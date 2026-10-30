"""Context-aware numerical matrices, independent of polyform exponentials."""

from .context import MatrixAxis, MatrixKind, MatrixRepresentation, LabeledMatrix
from ._numeric import MatrixTolerance

__all__ = ["MatrixAxis", "MatrixKind", "MatrixRepresentation", "MatrixTolerance", "LabeledMatrix"]
