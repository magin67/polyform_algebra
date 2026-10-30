"""Context-aware numerical matrices, independent of polyform exponentials."""

from .context import MatrixAxis, MatrixKind, MatrixRepresentation, LabeledMatrix
from ._numeric import MatrixTolerance
from .laplacian import Laplacian, DirichletResult, EliminationResult
from .results import (
    VertexValues, CircumsphereResult, SpanningTreeCentralityResult,
    EigenvectorCentralityResult, AnnihilatorsResult, BalancingResult,
)

__all__ = ["MatrixAxis", "MatrixKind", "MatrixRepresentation", "MatrixTolerance",
           "LabeledMatrix", "Laplacian", "DirichletResult", "EliminationResult", "VertexValues",
           "CircumsphereResult", "SpanningTreeCentralityResult", "EigenvectorCentralityResult",
           "AnnihilatorsResult", "BalancingResult"]
