"""The Petersen graph metric and its spanning-tree coefficient.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from math import isclose
import networkx as nx
import numpy as np
from polyform_algebra import Polyform, quForm


def main():
    """Compare the exponential's spanning coefficient with a Laplacian cofactor."""
    laplacian = Polyform.from_networkx(nx.petersen_graph())
    metric = laplacian.exp()
    matrix, _ = laplacian.to_matrix()
    trees = metric.potential(quForm.one())
    # Kirchhoff's theorem supplies an independent check, without SciPy.
    expected = np.linalg.det(matrix[1:,1:])
    assert isclose(trees, expected)
    print("Grade sizes:", metric.sizes())
    print("Spanning-tree coefficient:", trees)


if __name__ == "__main__":
    main()
