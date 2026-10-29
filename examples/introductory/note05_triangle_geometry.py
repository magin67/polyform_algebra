"""Geometry of the resistance triangle of unit-conductance K3.

PMG notes 5 and 9: K3 Gram determinant, boundary norm, and simplex area.
https://github.com/magin67/pmg-notes/blob/048491a28a3774ef315afd417168f8638b51e74d/pmg-library/en/Introduction/From%20Lengths%20to%20Areas%20and%20Volumes.md
https://github.com/magin67/pmg-notes/blob/048491a28a3774ef315afd417168f8638b51e74d/pmg-library/en/Introduction/Metric%20of%20Higher-Grade%20Objects.md
"""
from math import sqrt

import numpy as np

from polyform_algebra import Point, quForm


def calculate():
    """Evaluate first-grade products and a second-grade quadratic boundary form."""
    a, b, c = Point.create_list(["1", "2", "3"])
    f12, f13, f23 = (quForm([pair]) for pair in ((a, b), (a, c), (b, c)))
    metric = (f12 + f13 + f23).exp()
    u, v, w = b-a, c-a, c-b

    # @ constructs a scalar polyform; norm supplies its numerical metric value.
    uv = metric.norm(u @ v)
    gram = np.array([[metric.norm(u @ u), uv],
                     [metric.norm(v @ u), metric.norm(v @ v)]])
    # Evaluate the supported quadratic form directly, without higher-grade @.
    boundary_norm = metric.norm(quForm([(a, b, c)]))
    area_squared = boundary_norm / 4
    return {
        "tau": metric.potential(quForm.one()),
        "R12": metric.norm(f12),
        "R13": metric.norm(f13),
        "R23": metric.norm(f23),
        "uv": uv,
        "uw": metric.norm(u @ w),
        "vw": metric.norm(v @ w),
        "reversed_uv": metric.norm((a-b) @ v),
        "gram": gram,
        "gram_determinant": float(np.linalg.det(gram)),
        "boundary_norm": boundary_norm,
        "boundary_linear_norm": sqrt(boundary_norm),
        # The boundary's linear norm is the parallelogram area, twice the triangle area.
        "area_squared": area_squared,
        "area": sqrt(area_squared),
        "product_form_norm": metric.norm(f12*f13),
        "product_of_norms": metric.norm(f12)*metric.norm(f13),
    }


def main():
    """Print the geometry; unit conductances do not mean unit geometric sides."""
    print("PMG: geometry of the unit-conductance triangle K3")
    for name, value in calculate().items():
        print(f"{name}: {value}")


if __name__ == "__main__":
    main()
