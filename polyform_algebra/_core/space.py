"""Identity and dimension of a real linear space, independent of a metric."""

from dataclasses import dataclass
from numbers import Integral


@dataclass(frozen=True, eq=False)
class Space:
    """A space has identity: equal display names do not identify two spaces."""

    name: str
    dimension: int

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("space name must be a nonempty string")
        if isinstance(self.dimension, bool) or not isinstance(self.dimension, Integral):
            raise TypeError("space dimension must be an integer")
        if self.dimension < 0:
            raise ValueError("space dimension must be nonnegative")
