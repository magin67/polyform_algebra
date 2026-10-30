"""Eliminate a star's center and compare the resulting resistance geometry.

Source: PMG v5, note revision 5 (2026-10-15), read 2026-10-30.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Vertex%20Elimination%20and%20Induced%20Couplings.md
Section 'A star becomes a triangle', equations (3)-(6).
Run: python -m examples.elimination_star_pandas
"""

import numpy as np

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import Laplacian, MatrixAxis
from polyform_algebra.matrix.pandas import to_dataframe


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    axis = MatrixAxis(basis=Basis(["x", "a", "b", "c"], space=Space("star", 4)))
    star = Laplacian([[3, -1, -1, -1], [-1, 1, 0, 0],
                      [-1, 0, 1, 0], [-1, 0, 0, 1]], rows=axis)
    leaves = ["a", "b", "c"]
    result = star.eliminate(keep=leaves)
    triangle = result.laplacian
    print("Original star Laplacian:")
    print(to_dataframe(star).round(6))
    print("\nAfter eliminating x: triangle Laplacian:")
    print(to_dataframe(triangle).round(6))
    print("\nInduced adjacency (each edge has conductance 1/3):")
    print(to_dataframe(triangle.to_adjacency()).round(6))
    print("\nBi-coordinates of the projection of x onto the retained vertices:")
    print(to_dataframe(result.dirichlet.influence_matrix).round(6))

    before = to_dataframe(star.to_resistance()).loc[leaves, leaves]
    after = to_dataframe(triangle.to_resistance())
    print("\nLeaf resistances in the original star:")
    print(before.round(6))
    print("\nLeaf resistances after elimination:")
    print(after.round(6))
    # In the star each pair uses two unit resistances. In the triangle the
    # resistance-3 edge is parallel to a resistance-6 two-edge path.
    expected = 2 * (np.ones((3, 3)) - np.eye(3))
    np.testing.assert_allclose(before.to_numpy(), expected)
    np.testing.assert_allclose(after.to_numpy(), expected)
    np.testing.assert_allclose(triangle.to_adjacency().to_numpy(),
                               (np.ones((3, 3)) - np.eye(3)) / 3)

    prescribed = np.array([0., 3., 6.])
    extended = result.dirichlet.extend(dict(zip(leaves, prescribed)))
    full = np.array(list(extended.values()))
    original_response = star.to_numpy() @ full
    reduced_response = triangle.to_numpy() @ prescribed
    print("\nHarmonic extension:")
    print(pd.Series(extended, name="potential").to_frame())
    print("\nResponse at retained vertices:")
    print(pd.DataFrame({"original": original_response[1:],
                        "reduced": reduced_response}, index=leaves).round(6))
    np.testing.assert_allclose(original_response[0], 0, atol=1e-12)
    np.testing.assert_allclose(original_response[1:], reduced_response, atol=1e-12)
    # The extension is (3,0,3,6); energy is 9+0+9=18 in either graph.
    np.testing.assert_allclose(full @ star.to_numpy() @ full, 18)
    np.testing.assert_allclose(prescribed @ triangle.to_numpy() @ prescribed, 18)
    print("\nChecks passed: resistance, response and energy are preserved.")
    print("Edge conductance 1/3 is not 1/R=1/2: the third vertex provides another path.")


if __name__ == "__main__":
    main()
