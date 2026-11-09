"""Minimal extension points for external reads and result publication.

Payload schemas belong to the applied scenario. Connector implementations own
authentication, transport and source-specific identifiers, not graph weights.
"""
from abc import ABC, abstractmethod
from typing import Any, Mapping


class DataSource(ABC):
    """Read-only source; implementations must report incomplete snapshots explicitly."""

    @abstractmethod
    def read_snapshot(self, query: Mapping[str, Any]) -> Mapping[str, Any]:
        """Read a scenario-defined snapshot with stable object identifiers.

        Query and snapshot schemas are agreed by the source and application.
        A failed read must not masquerade as a successful empty dataset.
        """
        raise NotImplementedError


class ResultSink(ABC):
    """Explicit publication capability, separate from read-only data access."""

    @abstractmethod
    def publish(self, result: Mapping[str, Any], *, idempotency_key: str) -> Mapping[str, Any]:
        """Publish results and return a transport receipt.

        A concrete sink must handle retries for the same key without duplicate
        effects and reject reuse with a conflicting payload. The abstract
        contract alone does not implement persistence or idempotency.
        """
        raise NotImplementedError
