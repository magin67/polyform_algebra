"""Strict reference-core data contracts (stage two).

PMG v5 Core Specification, revision 27: underlying space, affine weight,
and aggregation of parallel conductances. No metric solver is exposed here.
"""
from .data import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, Calculation,
    Diagnostics, PMGError, InvalidInput, ContextMismatch, ArithmeticMismatch,
    NumericalUnresolved, coordinates, network, reorder, to_numeric,
)

__all__ = [
    'Context', 'Coordinates', 'MatrixData', 'Network', 'NumericPolicy',
    'Calculation', 'Diagnostics', 'PMGError', 'InvalidInput', 'ContextMismatch',
    'ArithmeticMismatch', 'NumericalUnresolved', 'coordinates', 'network',
    'reorder', 'to_numeric',
]
