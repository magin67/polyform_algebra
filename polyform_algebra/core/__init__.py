"""Strict reference-core data, network geometry, reductions and formal algebra.

PMG v5 Core Specification, revision 27: underlying space, affine weight,
and aggregation of parallel conductances; Introduction matrix/geometry notes.
"""
from .data import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, Calculation,
    Diagnostics, PMGError, InvalidInput, ContextMismatch, ArithmeticMismatch,
    NumericalUnresolved, DomainError, NonuniqueProblem, INFINITY, coordinates, network, reorder, to_numeric,
    UnsupportedOperation, ResourceLimit,
)
from .network_math import (
    PairValues, SolveCheck, CircumsphereResult, LaplacianCircumsphereResult,
    laplacian, green, resistance, circumsphere, degree_potential,
    laplacian_circumsphere, spanning_tree_weight,
    network_from_laplacian, network_from_resistance, inner_product, gram,
    blade_inner_product, simplex_volume_squared,
)
from .reductions import (
    DirichletResult, ReductionResult, ContractionResult,
    dirichlet, extend, kron_reduce, contract, induced_subgraph,
)
from .isotropic import IsotropicExtension
from .algebra import (
    AlgebraLimits, Exterior, SimpleTerm, Form, exterior, simple_term, form,
    Atom, Boundary, Unit, wedge, boundary, grade_part, square, form_product,
    partition_form, form_coefficients, equivalent, polar_form, form_exp, laplacian_form, spanning_coefficient,
    potential, evaluate_form,
)

__all__ = [
    'Context', 'Coordinates', 'MatrixData', 'Network', 'NumericPolicy',
    'Calculation', 'Diagnostics', 'PMGError', 'InvalidInput', 'ContextMismatch',
    'ArithmeticMismatch', 'NumericalUnresolved', 'coordinates', 'network',
    'reorder', 'to_numeric',
    'DomainError', 'INFINITY', 'PairValues', 'SolveCheck', 'laplacian', 'green',
    'resistance', 'CircumsphereResult', 'circumsphere', 'spanning_tree_weight', 'network_from_laplacian',
    'LaplacianCircumsphereResult', 'degree_potential', 'laplacian_circumsphere',
    'network_from_resistance', 'inner_product', 'gram', 'blade_inner_product',
    'simplex_volume_squared',
    'NonuniqueProblem', 'DirichletResult', 'ReductionResult', 'ContractionResult',
    'dirichlet', 'extend', 'kron_reduce', 'contract', 'induced_subgraph',
    'IsotropicExtension',
    'UnsupportedOperation', 'ResourceLimit', 'AlgebraLimits', 'Exterior',
    'SimpleTerm', 'Form', 'exterior', 'simple_term', 'form', 'Atom', 'Boundary',
    'Unit', 'wedge', 'boundary', 'grade_part', 'square', 'form_product',
    'partition_form', 'form_coefficients', 'equivalent', 'polar_form', 'form_exp', 'laplacian_form',
    'spanning_coefficient', 'potential', 'evaluate_form',
]
