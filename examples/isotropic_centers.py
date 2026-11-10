"""Reproduce the nonlocal-vertices example from the PMG Spaces series.

Source: PMG v5, Nonlocal Elements, Spheres, and Their Centers,
"Example: nonlocal vertices and an unchanged boundary norm".
https://github.com/magin67/pmg-notes/blob/6a80e80026413c7842e2d37de508b67c7ad7b687/pmg-library/en/Resistance%20Geometry%20and%20Duality/Spaces/Nonlocal%20Elements%2C%20Spheres%2C%20and%20Their%20Centers.md
"""

from polyform_algebra.core import IsotropicExtension, network


def triangle_boundary_norm(extension, a, b, c):
    """Compute the quadratic norm of the boundary via its Gram determinant."""
    u, v = b - a, c - a
    return extension.norm(u) * extension.norm(v) - extension.inner_product(u, v)**2


def main():
    path = network(('a', 'b', 'c'), (('a', 'b', 1), ('b', 'c', 1)))
    extension = IsotropicExtension(path)
    a, b, c = (extension.vertex(label) for label in path.context.labels)
    z = extension.z
    a_shifted, b_shifted, c_shifted = a + z, b - 2*z, c + 3*z

    for label, vertex, shifted in zip(path.context.labels, (a, b, c), (a_shifted, b_shifted, c_shifted)):
        center = extension.center(shifted)
        assert center == vertex and extension.norm(center) == 0
        print(f"{label}': norm = {extension.norm(shifted)}, center = {label}")

    shifted_norm = triangle_boundary_norm(extension, a_shifted, b_shifted, c_shifted)
    center_norm = triangle_boundary_norm(extension, a, b, c)
    assert shifted_norm == center_norm
    print(f'Boundary norm before and after centering: {shifted_norm}')


if __name__ == '__main__':
    main()
