"""Delete a cycle edge and a bridge using inverse edge factors.

Source: https://github.com/magin67/pmg-notes/blob/90a92a61044c79686925d0e32f09f150dfdd2c43/pmg-library/en/Introduction/Laplacian%20Exponential.md
Sections: adding and deleting an edge, disconnected graphs, and potentials.
"""
import networkx as nx

from polyform_algebra import Polyform, quForm


def calculate():
    """Keep the full vertex set when deletion isolates a vertex."""
    cycle = nx.cycle_graph(range(1, 5))
    chain = nx.Graph()
    chain.add_weighted_edges_from([(1, 2, 1), (2, 3, 2)])
    results = {}
    for name, graph, edge in [("cycle", cycle, (1, 2)), ("bridge", chain, (2, 3))]:
        laplacian, points = Polyform.from_networkx(graph, return_points="dict")
        metric = laplacian.exp()
        i, j = edge
        conductance = graph[i][j].get("weight", 1)
        removed_form = quForm([(points[i], points[j])])
        # A negative factor coefficient cancels an existing positive conductance.
        updated = metric * (Polyform.one() - conductance*removed_form)
        remaining = graph.copy()
        remaining.remove_edge(i, j)  # Preserve isolated vertices.
        rebuilt = Polyform.from_networkx(remaining).exp()
        tau_after = updated.potential(quForm.one())
        endpoint_form = quForm([(points[1], points[max(points)])])
        result = {
            "before": metric,
            "after": updated,
            "points": points,
            "removed_form": removed_form,
            "conductance": conductance,
            "tau_before": metric.potential(quForm.one()),
            "tau_after": tau_after,
            "removed_pair_R_before": metric.norm(removed_form),
            "deletion_ratio": 1-conductance*metric.norm(removed_form),
            "top_grade_after": updated.grade(),
            "endpoint_potential_after": updated.potential(endpoint_form),
            "rebuild_matches": updated.is_equivalent(rebuilt),
            "restoration_matches": (updated*(Polyform.one()+conductance*removed_form)).is_equivalent(metric),
        }
        # Do not catch arbitrary exceptions or report a disconnected norm as zero.
        if tau_after:
            result["removed_pair_R_after"] = updated.norm(removed_form)
        else:
            result["removed_pair_R_after"] = None
        results[name] = result
    return results


def main():
    """Show finite connected resistance and explain the unavailable disconnected normalization."""
    for name, result in calculate().items():
        print(f"{name}: delete an edge")
        for key, value in result.items():
            if key not in ("points", "removed_form"):
                print(f"  {key}: {value}")
        if result["tau_after"] == 0:
            print("  Global norm is undefined: zero spanning coefficient on the full vertex set.")


if __name__ == "__main__":
    main()
