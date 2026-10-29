"""Arithmetic means and affine weights in point and mixed bases.

Sources: PMG v5, terminology ('Centroid', 'Affine weight') and
Spaces/Projections, Bi-coordinates, and Di-coordinates ('Bi-coordinates and frames').
https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Spaces
"""

from polyform_algebra import Basis, Point, Vector


def main():
    """Show epsilon(mean)=sum epsilon(h_i)/n without choosing a metric."""
    a, b, c = Point.create_list("abc")
    centroid, frame = Basis([a, b, c]).centroid(name="h")
    print("Point centroid:", centroid, "=", frame)
    print("Affine weight:", centroid.affine_weight)
    mean, mixed_frame = Basis([a, Vector("u")]).centroid()
    print("Mixed mean:", mixed_frame)
    print("Affine weight:", mean.affine_weight)


if __name__ == "__main__":
    main()
