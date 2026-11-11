"""PMG Metric Distributions, the overlap of two paths in a weighted tree.

https://github.com/magin67/pmg-notes/tree/7cebaa55451d801c44be39ca1149c69bd00dd6fd/pmg-library/en/Resistance%20Geometry%20and%20Duality/Metric%20Distributions%20and%20Edge%20Influence
"""

from fractions import Fraction as Q

from polyform_algebra.core import network
from polyform_algebra.metric_distributions import metric_distribution, pair_blade


def main():
    path = network(tuple("abcde"), [("a", "b", 3), ("b", "c", 2), ("c", "d", Q(1, 2)), ("d", "e", 4)])
    data = metric_distribution(path)
    left = pair_blade(path, (("a", "d"),))
    right = pair_blade(path, (("b", "e"),))
    coefficients = data.coefficients(left, right)
    assert sum(coefficients.values()) == data.evaluate(left, right) == Q(5, 2)
    print("Overlap coefficients:", {a+b: str(value) for ((a, b),), value in coefficients.items()})


if __name__ == "__main__":
    main()
