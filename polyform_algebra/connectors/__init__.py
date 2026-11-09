"""External data exchange contracts, independent of graph mathematics.

Concrete transports (files, HTTP, 1C) implement these contracts separately.
Importing this package neither opens connections nor loads optional libraries.
"""
from .base import DataSource, ResultSink

__all__ = ['DataSource', 'ResultSink']
