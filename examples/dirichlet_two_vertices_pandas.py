"""Bi-coordinates of projections onto two selected vertices, displayed in pandas.

Run: python -m examples.dirichlet_two_vertices_pandas
Optional dependency: python -m pip install pandas
Source: PMG v5 revision d3b9cb30f8a990ebdc4c46b9cf4deb7f34d7a661,
'Selected Vertices and the Dirichlet Problem', weighted path and equation (7).
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Selected%20Vertices%20and%20the%20Dirichlet%20Problem.md
"""

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis
from polyform_algebra.matrix.pandas import to_dataframe


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    vertices = Space("weighted path", 4)
    axis = MatrixAxis(basis=Basis(["a", "x", "y", "b"], space=vertices))
    lap = Laplacian(
        [[2, -2, 0, 0], [-2, 3, -1, 0],
         [0, -1, 4, -3], [0, 0, -3, 3]],
        rows=axis,
    )
    print("Path: a --2-- x --1-- y --3-- b (conductances)")
    print("\nOriginal Laplacian:")
    print(to_dataframe(lap).to_string())

    problem = lap.dirichlet(selected=["a", "b"])
    coordinates = to_dataframe(problem.influence_matrix)
    print("\nSelected vertices: a, b")
    print("Bi-coordinates of the projections of x and y:")
    print(coordinates.to_string(float_format=lambda value: f"{value:.6f}"))
    # Solve [[3,-1],[-1,4]] B = diag(2,3); determinant 11 fixes the denominator.
    print("\nExact rows: x -> (8/11, 3/11); y -> (2/11, 9/11).")
    print("q_x = (8*a + 3*b)/11; q_y = (2*a + 9*b)/11.")
    print("These are projections onto the line a-b, not the original vertices.")

    prescribed = pd.Series({"a": 0.0, "b": 11.0}, name="prescribed")
    extended = pd.Series(problem.extend(prescribed.to_dict()), name="value")
    print("\nHarmonic extension of a=0, b=11:")
    print(extended.to_frame().to_string(float_format=lambda value: f"{value:.6f}"))
    print("At x: 2*(3-0) + (3-9) = 0; at y: (9-3) + 3*(9-11) = 0.")


if __name__ == "__main__":
    main()
