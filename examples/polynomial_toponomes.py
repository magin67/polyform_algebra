"""Vector and vertex polynomials and prescribed roots in PMG v5.

https://magin67.github.io/pmg-notes/en/polyform-algebra-and-combinatorics/polynomial-toponomes/
"""

from polyform_algebra.core import network
from polyform_algebra.toponomes import (
    coefficient_table, forest_toponome, laplacian_polynomial,
    rooted_forest_weight, vertex_toponome,
)


def main():
    # Forest Coefficients and the Laplacian Polynomial: a complete C4 example.
    cycle = network(tuple("abcd"), (("a", "b", 1), ("b", "c", 1), ("c", "d", 1), ("a", "d", 1)))
    table = coefficient_table(vertex_toponome(forest_toponome(cycle)))
    coefficients = tuple(sum(value for (k, _), value in table.items() if k == grade)
                         for grade in range(4))
    assert coefficients == laplacian_polynomial(cycle) == (1, 8, 20, 16)
    assert rooted_forest_weight(cycle, ("a", "c")) == 4
    print("P(C4) = 1 + 8t + 20t² + 16t³; roots a, c: 4")


if __name__ == "__main__":
    main()
