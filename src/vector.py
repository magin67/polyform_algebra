# @title Вектор

from .monoid import Monoid
from .affine_vector import  AffineVector

class Vector(Monoid):
    _cache = {}
    def __new__(cls, coeffs=None):
        if coeffs is None:
            obj = super().__new__(cls)
            obj.components = ()
            obj.multiplicity = -1
            return obj
        key = frozenset((sym, round(coeff, 10)) for sym, coeff in coeffs.items() if abs(coeff) > 1e-10)
        if key in cls._cache:
            return cls._cache[key]
        obj = super().__new__(cls)
        obj.components = ()
        obj._coeffs = dict(coeffs)
        obj._key = key
        obj.multiplicity = 0  # здесь устанавливаем
        cls._cache[key] = obj
        return obj

    def __init__(self, coeffs=None): # ничего не делаем, чтобы не перезаписывать multiplicity
        pass

    def dot(self, other):
        if not isinstance(other, (Vector, AffineVector)):
            return NotImplemented
        total = 0
        for sym, c in self._coeffs.items():
            total += c * other._coeffs.get(sym, 0)
        return total

    def __hash__(self):
        return hash(self._key)

    def __eq__(self, other):
        return isinstance(other, Vector) and self._key == other._key

    def __repr__(self):
        return f"Vector({self._coeffs})"