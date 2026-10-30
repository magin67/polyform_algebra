"""Context-aware numerical matrices, independent of polyform exponentials."""

from .context import MatrixAxis, MatrixKind, MatrixRepresentation, LabeledMatrix
from ._numeric import MatrixTolerance
from .laplacian import Laplacian, DirichletResult
from .results import (
    VertexValues, CircumsphereResult, SpanningTreeCentralityResult,
    EigenvectorCentralityResult, AnnihilatorsResult,
)

__all__ = ["MatrixAxis", "MatrixKind", "MatrixRepresentation", "MatrixTolerance",
           "LabeledMatrix", "Laplacian", "DirichletResult", "VertexValues",
           "CircumsphereResult", "SpanningTreeCentralityResult", "EigenvectorCentralityResult",
           "AnnihilatorsResult"]
