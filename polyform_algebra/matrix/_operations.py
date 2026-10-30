"""Implementation methods for the single public LabeledMatrix class."""

import numpy as np
from . import _numeric as num


class MatrixOperations:
    """Internal organization only; not a second public matrix type."""

    def vector_type(self, atol=None, rtol=None):
        """Left means 1.T @ A=0; right means A @ 1=0, in current coordinates."""
        return num.vector_type(self._array(), self._tol(atol, rtol))

    def is_symmetric(self, atol=None, rtol=None):
        a = self._array(finite=False)
        tol = self._tol(atol, rtol)
        if a.shape[0] != a.shape[1]:
            return False
        finite = np.isfinite(a)
        if not np.array_equal(finite, finite.T):
            return False
        return num.symmetric(np.where(finite, a, 0), tol)

    def is_skew_symmetric(self, atol=None, rtol=None):
        a = self._array()
        tol = self._tol(atol, rtol)
        return a.shape[0] == a.shape[1] and bool(np.all(np.abs(a + a.T) <= tol.threshold(num.scale(a))))

    def rank(self, atol=None, rtol=None):
        s = np.linalg.svd(self._array(), compute_uv=False)
        return int(np.count_nonzero(s > self._tol(atol, rtol).threshold(num.scale(s))))

    def is_singular(self, atol=None, rtol=None):
        a = self._array()
        if a.shape[0] != a.shape[1]:
            raise ValueError("singularity requires a square matrix; use rank() for rectangular data")
        return self.rank(atol, rtol) < len(a)

    def kernel(self, side="right", atol=None, rtol=None):
        """Euclidean-orthonormal coordinate columns; left kernel satisfies K.T @ A=0.

        This makes no orthonormality claim in a PMG metric.
        """
        if side not in ("left", "right"):
            raise ValueError("side must be left or right")
        a = self._array()
        if side == "left":
            a = a.T
        _, s, vh = np.linalg.svd(a, full_matrices=True)
        rank = np.count_nonzero(s > self._tol(atol, rtol).threshold(num.scale(s)))
        return vh[rank:].T.copy()

    def is_positive_semidefinite(self, atol=None, rtol=None):
        a = self._array()
        tol = self._tol(atol, rtol)
        if not num.symmetric(a, tol):
            raise ValueError("positive semidefiniteness requires a symmetric matrix")
        return bool(np.all(np.linalg.eigvalsh(a) >= -tol.threshold(num.scale(a))))

    def is_positive_definite(self, atol=None, rtol=None):
        a = self._array()
        tol = self._tol(atol, rtol)
        if not num.symmetric(a, tol):
            raise ValueError("positive definiteness requires a symmetric matrix")
        return bool(np.all(np.linalg.eigvalsh(a) > tol.threshold(num.scale(a))))

    def is_laplacian(self, atol=None, rtol=None):
        """Test the outgoing, nonnegative-conductance Laplacian convention."""
        a = self._array()
        tol = self._tol(atol, rtol)
        if a.shape[0] != a.shape[1]:
            return False
        try:
            num.validate_laplacian(a, tol)
        except ValueError:
            return False
        return True

    def pseudoinverse(self, atol=None, rtol=None):
        """Moore-Penrose inverse in Euclidean coordinates; a generic operator.

        Axes are exchanged. This does not silently resolve an ambiguous
        metric projection or assert a directed Greenian convention.
        """
        from .context import MatrixKind
        return self._new(num.pinv(self._array(), self._tol(atol, rtol)),
                         kind=MatrixKind.OPERATOR, rows=self.columns, columns=self.rows)

    def distance_transform(self):
        """D(A)_ij=A_ii+A_jj-A_ij-A_ji, only for vector-type matrices.

        Nonsymmetric input loses its skew part. Output is generic pair data,
        not a declaration of directed resistance or network realizability.
        """
        from .context import MatrixKind
        self._require_same_axes()
        a = self._array()
        if self.vector_type() != "vector":
            raise ValueError("distance transformation requires vector type on both sides")
        return self._new(num.distance(a), kind=MatrixKind.PAIRWISE)

    def deviation_transform(self):
        """-J R J/2 for finite symmetric zero-diagonal data.

        Inverts distance_transform on symmetric vector matrices only.
        """
        from .context import MatrixKind
        self._require_same_axes()
        a = self._array()
        if not self.is_symmetric() or np.any(np.abs(np.diag(a)) > self.tolerance.threshold(num.scale(a))):
            raise ValueError("deviation transformation requires symmetric zero-diagonal data")
        return self._new(num.deviation(a), kind=MatrixKind.FORM_COEFFICIENTS)

    def with_representation(self, representation, *, kind=None):
        """Declare meaning explicitly and validate it."""
        return self._new(self._array(finite=False), kind=kind, representation=representation)

    def _require_network(self):
        if self.representation == "generic":
            raise ValueError("declare a network representation before semantic conversion")
        if not self.is_symmetric():
            raise NotImplementedError("directed Greenian/resistance conventions are deferred; see MATRIX_DEFERRED.md")

    def _blockwise(self, function, fill=0.0):
        a = self._array(finite=False)
        result = np.full(a.shape, fill)
        for group in self._components:
            ix = np.ix_(group, group)
            result[ix] = function(a[ix])
        return result

    def to_greenian(self):
        """Dispatch by representation, with a separate centroid gauge per component."""
        from .context import MatrixKind
        self._require_network()
        if self.representation == "laplacian":
            for group in self._components:
                block = self._array()[np.ix_(group, group)]
                s = np.linalg.svd(block, compute_uv=False)
                if np.count_nonzero(s > self.tolerance.threshold(num.scale(s))) != len(group) - 1:
                    raise ValueError("network component is unresolved at this rank tolerance; choose a smaller tolerance")
            values = self._blockwise(lambda block: num.pinv(block, self.tolerance))
        elif self.representation == "resistance":
            values = self._blockwise(num.deviation)
        else:
            values = self._array()
        return self._new(values, kind=MatrixKind.FORM_COEFFICIENTS, representation="greenian")

    def to_laplacian(self):
        """Recover a symmetric network and validate its off-diagonal signs."""
        from .context import MatrixKind
        if self.representation == "laplacian":
            values = self._array()
        else:
            self._require_network()
            green = self.to_greenian()
            values = green._blockwise(lambda block: num.pinv(block, self.tolerance))
        return self._new(values, kind=MatrixKind.FORM_COEFFICIENTS, representation="laplacian")

    def to_resistance(self):
        """Within-component D(Gr), +inf between components, zero self-resistance."""
        from .context import MatrixKind
        self._require_network()
        if self.representation == "resistance":
            values = self._array(finite=False)
        else:
            values = self.to_greenian()._blockwise(num.distance, fill=np.inf)
        return self._new(values, kind=MatrixKind.PAIRWISE, representation="resistance")
