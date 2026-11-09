"""Compatibility imports; implementations live in polyform_algebra.adapters.matrix_networkx.

Existing callers and published PythonExamples retain the same functions.
"""
from ..adapters.matrix_networkx import from_networkx, to_networkx

__all__ = ['from_networkx', 'to_networkx']
