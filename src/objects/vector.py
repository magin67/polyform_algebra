# @title Вектор

from .._core.monoid import Monoid
from ..combinations.affine_vector import  AffineVector

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
        if isinstance(other, Vector):
            total = 0
            for sym, c in self._coeffs.items():
                total += c * other._coeffs.get(sym, 0)
            return total
        elif isinstance(other, AffineVector):
            total = 0
            for term, coeff in other.terms.items():
                # term — Simplex([sym])
                sym = term.components[0]
                total += coeff * self._coeffs.get(sym, 0)
            return total
        else:
            return NotImplemented

    def __hash__(self):
        return hash(self._key)

    def __eq__(self, other):
        return isinstance(other, Vector) and self._key == other._key

    def __repr__(self):
        return f"Vector({self._coeffs})"