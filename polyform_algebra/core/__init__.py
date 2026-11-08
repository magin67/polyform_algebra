"""Strict reference-core data and undirected network geometry (stages two/three).

PMG v5 Core Specification, revision 27: underlying space, affine weight,
and aggregation of parallel conductances; Introduction matrix/geometry notes.
"""
from .data import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, Calculation,
    Diagnostics, PMGError, InvalidInput, ContextMismatch, ArithmeticMismatch,
    NumericalUnresolved, DomainError, INFINITY, coordinates, network, reorder, to_numeric,
)
from .network_math import (
    PairValues, SolveCheck, laplacian, green, resistance, spanning_tree_weight,
    network_from_laplacian, network_from_resistance, inner_product, gram,
    blade_inner_product, simplex_volume_squared,
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
]
