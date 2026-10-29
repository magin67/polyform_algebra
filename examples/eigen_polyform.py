"""Cycle C5: potentials and norms before and after an eigenbasis change.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from math import isclose
import networkx as nx
from polyform_algebra import Polyform


def main():
    """Replace the obsolete PolyBForm/AffineVector example with the current API."""
    laplacian, vertices = Polyform.from_networkx(nx.cycle_graph(5), return_points="list")
    diagonal, eigenbasis = laplacian.to_eigenbasis()
    point_basis = eigenbasis.create_basis_from_frames()
    eigenbasis.link_to(point_basis)
    metric, eigenmetric = laplacian.exp(), diagonal.exp()
    vector = -vertices[0] + 2*vertices[1] - vertices[3]
    square = vector @ vector
    coordinates = vector.transform(point_basis, eigenbasis)
    diagonal_square = coordinates.to_diagonal_quform(eigenbasis)
    print("Eigenbasis Laplacian:", diagonal)
    print("Grade sizes:", metric.sizes(), eigenmetric.sizes())
    # Potentials depend on the top-form normalization; normalized norms agree.
    print("Potentials:", metric.potential(square), eigenmetric.potential(diagonal_square))
    original_norm = metric.norm(square)
    eigen_norm = eigenmetric.norm(diagonal_square)
    print("Norms:", original_norm, eigen_norm)
    assert isclose(original_norm, eigen_norm)


if __name__ == "__main__":
    main()
