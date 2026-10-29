"""Keep the explicit vertex order of the PMG weighted three-vertex chain.

Source: pmg-notes v5, Introduction/Laplacian, Green Matrix, and Effective
Resistance Matrix, section 'Example: a three-vertex path', revision
1e2a1221a4323b22e7fdbc96a4cb40254076cb5d.
https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Introduction
This example stores the stated Laplacian; it does not compute transformations.
"""

from polyform_algebra import Basis
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def main():
    """Store L=diag(A1)-A for chain conductances 1 and 2."""
    vertices = MatrixAxis(space="chain:vertex-values", basis=Basis([1, 2, 3]))
    laplacian = LabeledMatrix(
        [[1, -1, 0], [-1, 3, -2], [0, -2, 2]],
        rows=vertices, columns=vertices, kind=MatrixKind.FORM_COEFFICIENTS,
    )
    print("Vertex order:", laplacian.rows.labels)
    print(laplacian.to_numpy())


if __name__ == "__main__":
    main()
