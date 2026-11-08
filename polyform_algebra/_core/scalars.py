"""Shared rational arithmetic; legacy floating-point operations remain available."""
from fractions import Fraction
from numbers import Rational, Real


def is_scalar(value):
    return isinstance(value, Real)


def exact(value):
    """Normalize Python and SymPy rational numbers without a float conversion."""
    return Fraction(int(value.numerator), int(value.denominator))


def quotient(left, right):
    """Keep rational division exact, including integer division and factorials."""
    if isinstance(left, Rational) and isinstance(right, Rational):
        return exact(left) / exact(right)
    return left / right


def significant(value, tolerance):
    """Rational nonzero coefficients never disappear under a float tolerance."""
    return value != 0 if isinstance(value, Rational) else abs(value) > tolerance
