"""Compatibility imports; implementations live in polyform_algebra.adapters.pandas.

Existing callers and published PythonExamples retain the same functions.
"""
from ..adapters.pandas import from_dataframe, to_dataframe, to_series

__all__ = ['from_dataframe', 'to_dataframe', 'to_series']
