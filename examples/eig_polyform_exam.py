"""Compare the endpoint norm of a path in vertex and eigenvector coordinates.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from math import isclose
from polyform_algebra import Point, Polyform, quForm


def main():
    """A unit path of three edges has endpoint resistance three in either basis."""
    a,b,c,d = Point.create_list("abcd")
    laplacian = quForm([(a,b)]) + quForm([(b,c)]) + quForm([(c,d)])
    diagonal, eigenbasis = laplacian.to_eigenbasis()
    points = eigenbasis.create_basis_from_frames()
    eigenbasis.link_to(points)
    vector = d-a
    coordinates = vector.transform(points, eigenbasis)
    # The metric is diagonal here; mixed terms do not contribute to this norm.
    diagonal_square = coordinates.to_diagonal_quform(eigenbasis)
    original_norm = laplacian.exp().norm(vector @ vector)
    eigen_norm = diagonal.exp().norm(diagonal_square)
    print("Eigenbasis Laplacian:", diagonal)
    print("Expanded Laplacian:", diagonal.expand())
    print("Endpoint vector coordinates:", coordinates)
    print("Diagonal square:", diagonal_square)
    print("Vertex-basis norm:", original_norm)
    print("Eigenbasis norm:", eigen_norm)
    assert isclose(original_norm, 3)
    assert isclose(eigen_norm, original_norm)


if __name__ == "__main__":
    main()
