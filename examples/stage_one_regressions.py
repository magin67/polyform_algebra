"""Stage-one regression examples; PMG Core Specification v5, revision 27.

https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md
Sections: products of simplex boundaries and the resistance metric.
"""
from polyform_algebra import Basis, Element, Simplex, Space, get_version
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def main():
    print("Version:", get_version())
    axis = MatrixAxis(basis=Basis(["a", "b"], space=Space("V", 2)))
    lap = LabeledMatrix([[1e12, -1e12], [-1e12, 1e12]], rows=axis,
                        columns=axis, kind=MatrixKind.FORM_COEFFICIENTS,
                        representation="laplacian")
    print("Resistance (expected 1e-12):", lap.to_resistance().to_numpy()[0, 1])
    a, b, c = Element("a"), Element("b", affine_weight=2), Element("c")
    product = Simplex([a, b]).boundary() * Simplex([b, c]).boundary()
    print("Weighted boundary identity:", product == 2 * Simplex([a, b, c]).boundary())


if __name__ == "__main__":
    main()
