"""PMG note 7: strengthen edges 12 and 13 of a unit triangle by one each.

Source: https://github.com/magin67/pmg-notes/blob/v5/pmg-library/ru/Введение/Совместная%20вариация%20нескольких%20связей.md
Sections: triangle example, joint spanning coefficient and resistance between 2 and 3.
Source checked on 2026-10-29. NumPy only evaluates the note's comparison formulas.
"""
import numpy as np
from polyform_algebra import Point, Polyform, quForm


def calculate():
    """Compare two edge factors, sequential updates and the Gram-matrix formula."""
    a, b, c = Point.create_list(["1", "2", "3"])
    f12 = quForm([(a, b)])
    f13 = quForm([(a, c)])
    f23 = quForm([(b, c)])
    laplacian = f12 + f13 + f23
    metric = laplacian.exp()
    first_factor = Polyform.one() + f12
    second_factor = Polyform.one() + f13
    after_first = metric * first_factor
    updated = after_first * second_factor
    recomputed = (laplacian + f12 + f13).exp()

    # Both chosen edges point away from vertex 1; the measured edge is 2 -> 3.
    u, v, w = b-a, c-a, c-b
    gram = np.array([[metric.norm(u@u), metric.norm(u@v)],
                     [metric.norm(v@u), metric.norm(v@v)]])
    q = np.array([metric.norm(u@w), metric.norm(v@w)])
    tau = metric.potential(quForm.one())
    tau_after_first = after_first.potential(quForm.one())
    interaction = np.eye(2) + gram
    # The mixed coefficient is tau*det(G)=1, not the product of first derivatives.
    mixed_coefficient = metric.potential(f12*f13)
    return {
        "gram": gram,
        "q": q,
        "tau_before": tau,
        "tau_after_first": tau_after_first,
        "R13_after_first": after_first.norm(f13),
        "tau_after": updated.potential(quForm.one()),
        "R23_after": updated.norm(f23),
        "tau_ratio_formula": float(np.linalg.det(interaction)),
        "tau_formula": tau*float(np.linalg.det(interaction)),
        "R23_formula": metric.norm(f23)-float(q @ np.linalg.solve(interaction,q)),
        "tau_sequential_formula": tau_after_first*(1+after_first.norm(f13)),
        "mixed_coefficient": mixed_coefficient,
        "mixed_coefficient_formula": tau*float(np.linalg.det(gram)),
        # Both distinct unit edges are present in exactly one of three trees.
        "joint_probability_before": metric.norm(f12*f13),
        "factor_update_matches_exp": updated.is_equivalent(recomputed),
        "reverse_order_matches": updated.is_equivalent(metric*second_factor*first_factor),
    }


def main():
    """Print the joint variation and its independent matrix-formula comparison."""
    print("PMG note 7: strengthen edges 12 and 13 in K3 (both 1 -> 2)")
    for name, value in calculate().items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
