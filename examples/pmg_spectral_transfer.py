"""PMG Spectra of Isomers, transferring a passport to a new boundary.

Source: https://github.com/magin67/pmg-notes/tree/a55218dafc004ce1d8bd30dca60ea417d6eb5a74/pmg-library/en/Graph%20Transformations%20and%20Gluing/Spectra%20of%20Isomers
"""

from sympy import Matrix, Symbol, expand, factor

from polyform_algebra.boundary_data import spectral_pair_passport, spectral_pair_transfer
from polyform_algebra.core import network


def main():
    x = Symbol("x")
    entry, exit = ("a", "b"), ("d", "e")
    left_edges = [("a", "c", 1), ("c", "b", 1)]
    right_edges = [("a", "d", 1), ("d", "e", 1), ("e", "b", 1)]
    old = network(tuple("abc"), left_edges)
    kernel = network(tuple("abde"), right_edges)
    glued = network(tuple("abcde"), left_edges + right_edges)
    transfer = spectral_pair_transfer(kernel, entry, exit, x)
    result = (transfer * Matrix(spectral_pair_passport(old, entry, x))).applyfunc(expand)
    assert tuple(result) == spectral_pair_passport(glued, exit, x)
    print("Characteristic polynomial:", factor(result[0]))


if __name__ == "__main__":
    main()
