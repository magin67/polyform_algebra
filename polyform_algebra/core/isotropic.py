"""Exact isotropic extension of the resistance geometry of a connected graph.

The basis consists of the local graph vertices followed by a distinguished
null vector. The vertex block of the Gram matrix is minus half the resistance
matrix, and each local vertex pairs with the null vector to give one.
"""

from dataclasses import InitVar, dataclass, field

from .data import (
    Context, Coordinates, DomainError, InvalidInput, MatrixData, Network,
    UnsupportedOperation, _compatible, _label,
)
from .network_math import _require, resistance


@dataclass(frozen=True)
class IsotropicExtension:
    """Extend an exact connected network, preserving its vertex order.

    Use the resulting context with Coordinates or the exterior algebra. The
    null label must be distinct from the network's vertex labels. Norms are
    signed quadratic values, and a local center is defined for affine weight one.
    """

    network: InitVar[Network]
    null_label: object = 'z'
    context: Context = field(init=False)
    gram: MatrixData = field(init=False)

    def __post_init__(self, network):
        _require(network)
        if network.context.arithmetic != 'exact':
            raise UnsupportedOperation('isotropic extension requires exact arithmetic')
        if len(network.components) != 1:
            raise DomainError('isotropic extension requires a nonempty connected network')

        labels = network.context.labels
        context = Context(labels + (self.null_label,), epsilon=(1,) * len(labels) + (0,))
        resistances = resistance(network).value.values
        values = tuple(tuple(-value / 2 for value in row) + (1,) for row in resistances)
        values += ((1,) * len(labels) + (0,),)
        object.__setattr__(self, 'context', context)
        object.__setattr__(self, 'gram', MatrixData(values, context, context, 'gram'))

    @property
    def z(self):
        """Return the distinguished null vector of affine weight zero."""
        return Coordinates(self.context, (0,) * (len(self.context.labels) - 1) + (1,))

    def vertex(self, label):
        """Return the local point corresponding to a network vertex."""
        labels = self.context.labels
        if not _label(label) or label == self.null_label or label not in labels:
            raise InvalidInput('expected a network vertex label')
        return Coordinates(self.context, tuple(int(value == label) for value in labels))

    def _check_coordinates(self, value):
        if not isinstance(value, Coordinates):
            raise InvalidInput('expected Coordinates')
        _compatible(self.context, value.context)

    def inner_product(self, left, right):
        """Return the bilinear product in the extended Gram matrix."""
        self._check_coordinates(left)
        self._check_coordinates(right)
        return sum(
            coefficient * sum(entry * value for entry, value in zip(row, right.values))
            for coefficient, row in zip(left.values, self.gram.values)
        )

    def affine_weight(self, value):
        """Return the pairing with z, which sums the vertex coefficients."""
        self._check_coordinates(value)
        return sum(value.values[:-1])

    def norm(self, value):
        """Return the signed quadratic norm of a point, element, or vector."""
        return self.inner_product(value, value)

    def center(self, element):
        """Return the unique local point on the element's line parallel to z."""
        if self.affine_weight(element) != 1:
            raise DomainError('center requires an element of affine weight one')
        return element - (self.norm(element) / 2) * self.z
