"""PMG note 10: grade expansions, tree coefficients, and resistances for P5 and C4.

Source: https://github.com/magin67/pmg-notes/blob/33a873c5e62ec9d5e5afe94086009223c129cfcb/pmg-library/en/Introduction/Laplacian%20Exponential.md
Sections: five-vertex path, four-cycle, and calculation of potentials and norms.
"""
from itertools import combinations

import networkx as nx

from polyform_algebra import Polyform, quForm


def calculate():
    """Return both unit-conductance examples, retaining their forest representations."""
    graphs = {
        "P5": nx.path_graph(range(1, 6)),
        "C4": nx.cycle_graph(range(1, 5)),
    }
    results = {}
    for name, graph in graphs.items():
        laplacian, points = Polyform.from_networkx(graph, return_points="dict")
        metric = laplacian.exp()
        # Both graphs have four edges; include grade 4 to display the cycle's zero.
        grades = [metric.extract_grade(k) for k in range(5)]
        factors = Polyform.one()
        for i, j in graph.edges():
            factors = factors * (Polyform.one() + quForm([(points[i], points[j])]))
        pairs = list(combinations(graph.nodes(), 2))
        forms = {pair: quForm([(points[pair[0]], points[pair[1]])]) for pair in pairs}
        results[name] = {
            "laplacian": laplacian,
            "metric": metric,
            "points": points,
            "grades": grades,
            # These sums refer to the retained forest expansion, not an arbitrary basis.
            "forest_weight_sums": [sum(part.terms.values()) for part in grades],
            "stored_term_counts": [len(part.terms) for part in grades],
            "tau": metric.potential(quForm.one()),
            "potentials": {pair: metric.potential(form) for pair, form in forms.items()},
            "resistances": {pair: metric.norm(form) for pair, form in forms.items()},
            "factorization_matches": metric.is_equivalent(factors),
        }
    return results


def main():
    """Print every grade, forest weight sum, and pair resistance for the two graphs."""
    for name, result in calculate().items():
        print(f"{name}: unit conductances")
        for k, part in enumerate(result["grades"]):
            print(f"  M{k} = {part}")
        for key in ("forest_weight_sums", "stored_term_counts", "tau",
                    "potentials", "resistances", "factorization_matches"):
            print(f"  {key}: {result[key]}")


if __name__ == "__main__":
    main()
