"""Context-aware numerical matrices, independent of polyform exponentials."""

from .context import MatrixAxis, MatrixKind, MatrixRepresentation, LabeledMatrix
from ._numeric import MatrixTolerance
from .laplacian import Laplacian, DirichletResult

__all__ = ["MatrixAxis", "MatrixKind", "MatrixRepresentation", "MatrixTolerance",
           "LabeledMatrix", "Laplacian", "DirichletResult"]
