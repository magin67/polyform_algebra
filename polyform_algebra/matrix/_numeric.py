"""Shared tolerances and dense coordinate algorithms (no polyform operations)."""

from dataclasses import dataclass
from numbers import Real
import numpy as np
from .._core.topology import connected_components


@dataclass(frozen=True)
class MatrixTolerance:
    """Zero threshold = atol + rtol * explicitly specified reference scale."""

    atol: float = 1e-10
    rtol: float = 1e-12

    def __post_init__(self):
        for value in (self.atol, self.rtol):
            if isinstance(value, bool) or not isinstance(value, Real):
                raise ValueError("tolerances must be finite nonnegative real numbers")
            if not np.isfinite(value) or value < 0:
                raise ValueError("tolerances must be finite nonnegative real numbers")

    def threshold(self, scale=0):
        return self.atol + self.rtol * scale


def scale(a):
    return float(np.max(np.abs(a))) if a.size else 0.0


def symmetric(a, tol):
    return a.shape[0] == a.shape[1] and bool(
        np.all(np.abs(a - a.T) <= tol.threshold(scale(a))))


def vector_type(a, tol):
    # Compare each sum against its own absolute sum, avoiding a large unrelated
    # row masking a small but meaningful imbalance elsewhere in the matrix.
    def zero_sums(axis):
        # Normalize before summation: finite inputs can otherwise overflow
        # both the residual and its threshold, making inf <= inf look valid.
        unit = np.max(np.abs(a), axis=axis, initial=0, keepdims=True)
        unit = np.where(unit == 0, 1.0, unit)
        scaled = a / unit
        with np.errstate(over="ignore"):
            threshold = tol.atol / unit.squeeze(axis) + tol.rtol * np.abs(scaled).sum(axis=axis)
        return np.all(np.abs(scaled.sum(axis=axis)) <= threshold)

    left = zero_sums(0)
    right = zero_sums(1)
    if left and right:
        return "vector"
    if left:
        return "left_vector"
    if right:
        return "right_vector"
    return "non_vector"


def pinv(a, tol):
    """Moore-Penrose inverse using the same singular cutoff as rank/kernel."""
    u, s, vh = np.linalg.svd(a, full_matrices=False)
    inverse = np.zeros_like(s)
    keep = s > tol.threshold(scale(s))
    inverse[keep] = 1 / s[keep]
    return (vh.T * inverse) @ u.T


def deviation(a):
    """Double center without constructing the dense centering projector."""
    if not a.size:
        return a.copy()
    return -0.5 * (a - a.mean(axis=0)[None, :]
                   - a.mean(axis=1)[:, None] + a.mean())


def distance(a):
    """Symmetric pair data; the skew part of a nonsymmetric input is lost."""
    d = np.diag(a)
    return d[:, None] + d[None, :] - a - a.T


def components(adjacency):
    """Deterministic connected components of a Boolean undirected relation."""
    return connected_components(tuple(range(len(adjacency))),
                                ((i, int(j)) for i, row in enumerate(adjacency)
                                 for j in np.flatnonzero(row)))


def validate_laplacian(a, tol):
    if not np.isfinite(a).all():
        raise ValueError("laplacian must be finite")
    if vector_type(a, tol) not in ("vector", "right_vector"):
        raise ValueError("outgoing Laplacian must have zero row sums")
    off = a.copy()
    np.fill_diagonal(off, 0)
    if np.any(off > tol.threshold(scale(a))) or np.any(np.diag(a) < -tol.threshold(scale(a))):
        raise ValueError("Laplacian requires nonpositive off-diagonal conductances")


def validate_greenian(a, tol):
    """Validate inverse networks without applying conductance units to Gr.

    A connected centered Greenian has rank n-1, so its nonzero support
    cannot split into smaller centered blocks. Preserve every nonzero link;
    reject unresolved rank instead of inventing disconnected components.
    Each block is normalized separately to allow different electrical scales.
    """
    if not np.isfinite(a).all():
        raise ValueError("greenian must be finite")
    groups = components((a != 0) | (a.T != 0))
    for group in groups:
        block = a[np.ix_(group, group)]
        unit = scale(block)
        relative = MatrixTolerance(0, max(tol.rtol, 100 * len(group) * np.finfo(float).eps))
        normalized = block / unit if unit else block
        if not symmetric(normalized, relative) or vector_type(normalized, relative) != "vector":
            raise ValueError("greenian must be symmetric and vector type")
        if np.any(np.linalg.eigvalsh(normalized) < -relative.threshold(scale(normalized))):
            raise ValueError("greenian must be positive semidefinite")
        s = np.linalg.svd(normalized, compute_uv=False)
        if np.count_nonzero(s > relative.threshold(scale(s))) != len(group) - 1:
            raise ValueError("greenian component is unresolved at this rank tolerance")
        validate_laplacian(pinv(normalized, relative), relative)
    return groups


def inverse_greenian_block(a, tol):
    """Invert a validated Gr block in its own units, then restore conductance."""
    unit = scale(a)
    if not unit:
        return np.zeros_like(a)
    relative = MatrixTolerance(0, max(tol.rtol, 100 * len(a) * np.finfo(float).eps))
    result = pinv(a / unit, relative) / unit
    if not np.isfinite(result).all():
        raise ValueError("inverse Greenian is outside the finite numerical range")
    return result


def validate_resistance(a, tol):
    """Finite blocks are networks; +inf entries must separate entire components."""
    if np.isnan(a).any() or np.isneginf(a).any():
        raise ValueError("resistance permits only finite entries or positive infinity")
    finite = np.isfinite(a)
    if not np.array_equal(finite, finite.T):
        raise ValueError("resistance infinity pattern must be symmetric")
    if not np.isfinite(np.diag(a)).all():
        raise ValueError("resistance diagonal must be finite and zero")
    groups = components(finite)
    for group in groups:
        block = a[np.ix_(group, group)]
        if not np.isfinite(block).all():
            raise ValueError("finite resistance pairs must form complete component blocks")
        eps = tol.threshold(scale(block))
        if not symmetric(block, tol) or np.any(np.abs(np.diag(block)) > eps) or np.any(block < -eps):
            raise ValueError("resistance must be symmetric, nonnegative and zero-diagonal")
        inner = validate_greenian(deviation(block), tol)
        if len(inner) != 1:
            raise ValueError("each finite resistance block must describe a connected network")
    return groups
