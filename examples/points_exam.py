"""Points, affine vectors and boundaries from the PMG introductory series.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from polyform_algebra import Point, Simplex, Vector


def main():
    """Demonstrate stored frames, boundary signs and first-grade polarization."""
    a, b, c = Point.create_list("abc")
    midpoint = Point("d", frame=(a + b) / 2)
    v = Vector("v", frame=b - a)
    u = Vector("u", frame=2*c - b - a)
    print("Point product:", a*b)
    print("Midpoint frame:", midpoint[0])
    print("Midpoint times a:", midpoint[0]*a[0])
    print("Boundary of ab:", Simplex([a, b]).boundary())
    # The vector has affine weight zero, so d(av)=v.
    print("Boundary of av:", Simplex([a, v]).boundary())
    print("Boundary times endpoint:", Simplex([(a, b), a]))
    print("Expanded boundary product:", Simplex([(a, b), c]).to_polysimplex())
    print("Squared v:", v[0] @ v[0])
    print("Squared u:", u[0] @ u[0])
    print("Scalar polyform v,u:", v[0] @ u[0])
    print("Squared midpoint:", midpoint[0] @ midpoint[0])
    print("Scalar polyform v,midpoint:", v[0] @ midpoint[0])


if __name__ == "__main__":
    main()
