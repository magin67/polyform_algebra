"""The four-vertex path in PMG's KVP Decomposition of Pair Characteristics."""

from sympy import Matrix, Rational
from polyform_algebra.kvp import decompose


def main():
    table = Matrix([[0, 1, 0, 0], [1, 0, 1, 0], [0, 1, 0, 1], [0, 0, 1, 0]])
    parts = decompose(table)
    assert parts.K[0, 1] == Rational(1, 2)
    assert parts.K + parts.V + parts.P == table
    print(parts.P)


if __name__ == "__main__":
    main()
