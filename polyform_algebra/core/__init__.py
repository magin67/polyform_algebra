"""Strict reference-core data, network geometry and reductions (stages two-four).

PMG v5 Core Specification, revision 27: underlying space, affine weight,
and aggregation of parallel conductances; Introduction matrix/geometry notes.
"""
from .data import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, Calculation,
    Diagnostics, PMGError, InvalidInput, ContextMismatch, ArithmeticMismatch,
    NumericalUnresolved, DomainError, NonuniqueProblem, INFINITY, coordinates, network, reorder, to_numeric,
)
from .network_math import (
    PairValues, SolveCheck, laplacian, green, resistance, spanning_tree_weight,
    network_from_laplacian, network_from_resistance, inner_product, gram,
    blade_inner_product, simplex_volume_squared,
)
from .reductions import (
    DirichletResult, ReductionResult, ContractionResult,
    dirichlet, extend, kron_reduce, contract, induced_subgraph,
)

__all__ = [
    'Context', 'Coordinates', 'MatrixData', 'Network', 'NumericPolicy',
    'Calculation', 'Diagnostics', 'PMGError', 'InvalidInput', 'ContextMismatch',
    'ArithmeticMismatch', 'NumericalUnresolved', 'coordinates', 'network',
    'reorder', 'to_numeric',
    'DomainError', 'INFINITY', 'PairValues', 'SolveCheck', 'laplacian', 'green',
    'resistance', 'spanning_tree_weight', 'network_from_laplacian',
    'network_from_resistance', 'inner_product', 'gram', 'blade_inner_product',
    'simplex_volume_squared',
    'NonuniqueProblem', 'DirichletResult', 'ReductionResult', 'ContractionResult',
    'dirichlet', 'extend', 'kron_reduce', 'contract', 'induced_subgraph',
]
