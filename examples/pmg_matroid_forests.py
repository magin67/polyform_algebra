"""PMG Matroid Parameters series, the weighted tree and resistance radius examples.

Source: https://github.com/magin67/pmg-notes/tree/8f30a5b9440c20ed841b50c9e0ce40241dd65e02/pmg-library/en/Invariants%2C%20Decompositions%2C%20and%20Spectra/Matroid%20parameters%20and%20forest%20characteristics%20of%20graphs
"""

from fractions import Fraction

from polyform_algebra.core import circumsphere, network
from polyform_algebra.matroid import (
    forest_coefficients, forest_power_sums, forest_resistance_radius,
    normalized_forest_coefficients,
)


def main():
    tree = network(tuple("abcd"), [("a", "b", 1), ("b", "c", 2), ("c", "d", 3)])
    coefficients = forest_coefficients(tree)
    powers = forest_power_sums(normalized_forest_coefficients(tree), 3)
    radius = forest_resistance_radius(tree)
    assert coefficients == (1, 6, 11, 6)
    assert radius == circumsphere(tree).value.resistance_radius == Fraction(11, 24)
    print("Forest coefficients:", tuple(map(str, coefficients)))
    print("Power sums:", tuple(map(str, powers)))
    print("Resistance radius:", radius)


if __name__ == "__main__":
    main()
