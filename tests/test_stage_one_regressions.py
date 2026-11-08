"""Regressions for PMG reference-core stage one.

Source: v5 Core Specification, revision 27 (2026-11-08), sections
'Products of simplex boundaries' and Laplacian/Greenian/resistance:
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md
"""
from importlib import metadata

import numpy as np
import pytest

import polyform_algebra as package
from polyform_algebra import Basis, Element, Simplex, Space
from polyform_algebra.matrix import LabeledMatrix, MatrixAxis, MatrixKind


def matrix(values, representation="laplacian"):
    axis = MatrixAxis(basis=Basis(range(len(values)), space=Space("V", len(values))))
    kind = MatrixKind.PAIRWISE if representation == "resistance" else MatrixKind.FORM_COEFFICIENTS
    return LabeledMatrix(values, rows=axis, columns=axis, kind=kind, representation=representation)


@pytest.mark.parametrize("c", [1., 1e6, 1e12])
def test_large_conductance_keeps_finite_resistance(c):
    """A single edge has R=1/c; small Gr must never become two isolates."""
    lap = matrix(c * np.array([[1., -1.], [-1., 1.]]))
    green = lap.to_greenian()
    assert green.components == ((0, 1),)
    resistance = green.to_resistance()
    assert resistance.to_numpy()[0, 1] == pytest.approx(1/c, rel=1e-10, abs=0)
    np.testing.assert_allclose(resistance.to_laplacian().to_numpy(), lap.to_numpy(), rtol=1e-10, atol=0)


def test_disconnected_greenian_uses_separate_scales():
    """Edges of different scales and an isolate retain exactly three components."""
    lap = np.zeros((5, 5))
    lap[:2, :2] = [[1, -1], [-1, 1]]
    lap[2:4, 2:4] = 1e12 * np.array([[1, -1], [-1, 1]])
    result = matrix(lap).to_resistance()
    assert result.components == ((0, 1), (2, 3), (4,))
    assert result.to_numpy()[2, 3] == pytest.approx(1e-12, rel=1e-10, abs=0)
    assert np.isinf(result.to_numpy()[0, 2])
    np.testing.assert_allclose(result.to_laplacian().to_numpy(), lap, rtol=1e-10, atol=0)


def test_unresolved_greenian_rejects_instead_of_splitting():
    """A weak nonzero mode is an unresolved connected block, not new components."""
    u = np.array([1., -1., 0.]) / np.sqrt(2)
    v = np.array([1., 1., -2.]) / np.sqrt(6)
    with pytest.raises(ValueError, match="unresolved"):
        matrix(np.outer(u, u) + 1e-15 * np.outer(v, v), "greenian")


@pytest.mark.parametrize("weight", [0, 2, -1])
def test_compact_boundary_rejects_nonunit_weights(weight):
    """Compact merging is restricted to unit weights, not a general epsilon rule."""
    a, b = Element("a"), Element("b", affine_weight=weight)
    with pytest.raises(ValueError, match="affine weight one"):
        Simplex([(a, b)])
    with pytest.raises(ValueError, match="affine weight one"):
        Simplex([(b,)])


def test_weighted_coordinate_boundary_keeps_pivot_factor():
    """For epsilon(b)=2, d(ab)d(bc)=2d(abc), preventing the missing factor two."""
    a, b, c = Element("a"), Element("b", affine_weight=2), Element("c")
    left = Simplex([a, b]).boundary() * Simplex([b, c]).boundary()
    assert left == 2 * Simplex([a, b, c]).boundary()


def test_version_matches_installed_metadata():
    """Both public entry points read the same installed distribution metadata."""
    assert package.get_version() == package.__version__ == metadata.version("polyform_algebra")


def test_missing_metadata_has_no_fabricated_version(monkeypatch):
    """An uninstalled checkout reports missing metadata without hiding its classes."""
    def missing(name):
        raise metadata.PackageNotFoundError(name)
    monkeypatch.setattr(metadata, "version", missing)
    assert package.Point is not None
    with pytest.raises(metadata.PackageNotFoundError):
        package.get_version()
    with pytest.raises(metadata.PackageNotFoundError):
        _ = package.__version__
