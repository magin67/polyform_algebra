"""PMG Spectra and Diagonalization, the path P3 example.

https://github.com/magin67/pmg-notes/tree/0f478b36dcc35c20a06c8d6c3041493996742910/pmg-library/ru/Инварианты%2C%20разложения%20и%20спектры
"""

from sympy import Matrix

from polyform_algebra.core import green, laplacian, network
from polyform_algebra.spectral import laplacian_projectors


def main():
    path = network(tuple("abc"), [("a", "b", 1), ("b", "c", 1)])
    blocks = laplacian_projectors(path)
    assert blocks[1] + 3 * blocks[3] == Matrix(laplacian(path).value.values)
    assert blocks[1] + blocks[3] / 3 == Matrix(green(path).value.values)
    print("Path eigenvalues:", sorted(blocks))


if __name__ == "__main__":
    main()
