"""Turn a pandas employee/project comment table into employee couplings.

Run: python -m examples.bipartite_comments_pandas
Optional dependency: python -m pip install pandas

The illustrative comment counts come from the accompanying discussion.
The reduction follows PMG, 'Couplings Within One Part of a Bipartite Graph',
equations (3), (6), revision f4e80ba6d7bec8d3e32774b957fe921bf76de264:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Spaces/Couplings%20Within%20One%20Part%20of%20a%20Bipartite%20Graph.md
"""

from polyform_algebra import Basis, Space
from polyform_algebra.matrix import MatrixAxis, MatrixKind
from polyform_algebra.matrix.pandas import from_dataframe, to_dataframe


def main():
    try:
        import pandas as pd
    except ImportError:
        print("This example requires pandas: python -m pip install pandas")
        return

    comments = pd.DataFrame(
        {"Проект X": [4, 2, 0], "Проект Y": [2, 0, 2]},
        index=["Анна", "Борис", "Вера"],
    )
    employees = Space("employees", len(comments.index))
    projects = Space("projects", len(comments.columns))
    rows = MatrixAxis(basis=Basis(comments.index, space=employees))
    columns = MatrixAxis(basis=Basis(comments.columns, space=projects))
    participation = from_dataframe(
        comments, rows=rows, columns=columns,
        kind=MatrixKind.PAIRWISE, representation="bipartite",
    )

    # Counts are modeled as conductances. Eliminate project vertices while
    # preserving the original network's effective employee resistances.
    laplacian = participation.to_unipartite(employees)
    connections = to_dataframe(laplacian.to_adjacency())
    laplacian_frame = to_dataframe(laplacian)

    print("Original comments by employee and project:")
    print(comments.to_string())
    print("\nEmployee couplings (adjacency):")
    # Formatting affects only the display, not the stored calculation.
    print(connections.to_string(float_format=lambda value: f"{value:.3f}"))
    print("\nEmployee Laplacian:")
    print(laplacian_frame.to_string(float_format=lambda value: f"{value:.3f}"))

    print("\nEach project contributes H[i,p] * H[j,p] / sum(H[:,p]).")
    print("Anna - Boris: 4 * 2 / (4 + 2) = 4/3.")
    print("Anna - Vera:  2 * 2 / (2 + 2) = 1.")
    print("Boris - Vera: no shared project; direct coupling is zero.")
    print("These weights are conductances, not comment counts or probabilities.")
    print("Indirect connections remain even when a direct coupling is zero.")
    print("\nOriginal activity totals (distinct from reduced graph degrees):")
    print(comments.sum(axis=1).to_string())


if __name__ == "__main__":
    main()
