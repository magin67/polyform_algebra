"""The four-vertex path from Laplacian and Greenian in PMG v5."""

from fractions import Fraction

from polyform_algebra.core import form_exp, laplacian_form, network, partition_form
from polyform_algebra.duality import affine_coefficients, dual_coefficients, greenian_form


def main():
    graph = network(tuple("abcd"), (("a", "b", 1), ("b", "c", 1), ("c", "d", 1)))
    direct = form_exp(laplacian_form(graph).value).value
    inverse = form_exp(greenian_form(graph).value).value
    top = partition_form(graph.context, (graph.context.labels,)).value
    assert dual_coefficients(direct, top * Fraction(1, 16)) == affine_coefficients(inverse)
    print("P4: exp(Gr) is the complementary image of exp(L)")


if __name__ == "__main__":
    main()
