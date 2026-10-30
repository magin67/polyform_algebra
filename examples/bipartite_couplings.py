"""Reproduce the PMG employee/document example using labeled numeric matrices.

Source: https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Couplings%20Within%20One%20Part%20of%20a%20Bipartite%20Graph.md
Equations (3), (6), revision f4e80ba6d7bec8d3e32774b957fe921bf76de264.
Run: python -m examples.bipartite_couplings
"""

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def main():
    employees = Space("employees", 3)
    documents = Space("documents", 4)
    rows = MatrixAxis(basis=Basis(["Andrey", "Boris", "Vladimir"], space=employees))
    columns = MatrixAxis(basis=Basis(["x1", "x2", "x3", "x4"], space=documents))
    participation = LabeledMatrix(
        [[1, 1, 0, 0], [1, 1, 1, 0], [0, 1, 1, 1]],
        rows=rows, columns=columns, kind=MatrixKind.PAIRWISE,
        representation="bipartite",
    )
    employee_laplacian = participation.to_unipartite(employees)
    document_adjacency = participation.to_unipartite(documents, representation="adjacency")
    # The shared documents contribute ab=bc=1/2+1/3 and ac=1/3.
    print("Employees:", employee_laplacian.rows.labels)
    print("Employee Laplacian:\n", employee_laplacian.to_numpy())
    print("Employee adjacency:\n", employee_laplacian.to_adjacency().to_numpy())
    print("Employee resistance:\n", employee_laplacian.to_resistance().to_numpy())
    print("Documents:", document_adjacency.rows.labels)
    print("Document adjacency:\n", document_adjacency.to_numpy())
    # Activity is measured in the input, not by the reduced graph degrees.
    print("Original employee activity:", participation.to_numpy().sum(axis=1))


if __name__ == "__main__":
    main()
