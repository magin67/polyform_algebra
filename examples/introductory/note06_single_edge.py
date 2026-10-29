"""PMG note 6: strengthen edge 12 of a unit triangle from conductance 1 to 2.

Source: https://github.com/magin67/pmg-notes/blob/v5/pmg-library/ru/Введение/Вариация%20одной%20связи.md
Sections: triangle example, spanning coefficient, edge probability and forest weight.
Source checked on 2026-10-29. All calculations stay in the quadratic algebra.
"""
from polyform_algebra import Point, Polyform, quForm


def calculate():
    """Compute the note's numbers by multiplying the metric by one edge factor."""
    a, b, c = Point.create_list(["1", "2", "3"])
    f12 = quForm([(a, b)])
    f13 = quForm([(a, c)])
    f23 = quForm([(b, c)])
    laplacian = f12 + f13 + f23
    metric = laplacian.exp()
    delta = 1

    # f12*f12=0, so exp(L+delta*f12)=exp(L)*(1+delta*f12).
    updated = metric * (Polyform.one() + delta*f12)
    recomputed = (laplacian + delta*f12).exp()
    tau = metric.potential(quForm.one())
    tau_updated = updated.potential(quForm.one())
    r12 = metric.norm(f12)
    r13 = metric.norm(f13)
    coupling = metric.norm((b-a) @ (c-a))
    return {
        "tau_before": tau,
        "R12_before": r12,
        "R13_before": r13,
        "R23_before": metric.norm(f23),
        "coupling": coupling,
        "tau_after": tau_updated,
        "R12_after": updated.norm(f12),
        "R13_after": updated.norm(f13),
        "R23_after": updated.norm(f23),
        # P(edge in T)=c*R: the changed edge now has conductance two.
        "edge_probability_before": r12,
        "edge_probability_after": (1+delta)*updated.norm(f12),
        # The unnormalized potential counts weighted spanning 2-forests.
        "forest13_before": metric.potential(f13),
        "forest13_after": updated.potential(f13),
        "forest13_change": updated.potential(f13)-metric.potential(f13),
        # Compare with the closed formulas in the note, using the old metric.
        "tau_formula": tau*(1+delta*r12),
        "R13_formula": r13-delta*coupling**2/(1+delta*r12),
        "forest13_change_formula": delta*tau*(r12*r13-coupling**2),
        "factor_update_matches_exp": updated.is_equivalent(recomputed),
    }


def main():
    """Print the before/after quantities; regression tests check every result."""
    print("PMG note 6: strengthen edge 12 in K3 (conductance 1 -> 2)")
    for name, value in calculate().items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
