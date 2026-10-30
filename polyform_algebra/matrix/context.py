"""Shared matrix storage with explicit roles, representations and axis context."""

from dataclasses import dataclass, field
from enum import Enum
import numpy as np

from .._core.basis import Basis
from .._core.space import Space
from . import _numeric as num
from ._numeric import MatrixTolerance
from ._operations import MatrixOperations


class MatrixKind(Enum):
    """Tensor roles determine basis-change rules, independently of representation."""
    OPERATOR = "operator"
    METRIC = "metric"
    FORM_COEFFICIENTS = "form_coefficients"
    PAIRWISE = "pairwise"


class MatrixRepresentation(str, Enum):
    """Declared meanings, validated rather than inferred from numerical coincidences."""
    GENERIC = "generic"
    ADJACENCY = "adjacency"
    BIPARTITE = "bipartite"
    LAPLACIAN = "laplacian"
    GREENIAN = "greenian"
    RESISTANCE = "resistance"


@dataclass(frozen=True, eq=False, init=False)
class MatrixAxis:
    """Snapshot of an existing Basis; the Basis owns the Space.

    Legacy space= binds an unbound Basis once. A string creates a Space;
    reuse that Space instance for other bases of the same space.
    """
    basis: Basis
    _labels: tuple = field(repr=False)

    def __init__(self, space=None, basis=None):
        if not isinstance(basis, Basis):
            raise TypeError("basis must be a Basis")
        if isinstance(space, str):
            if not space.strip():
                raise ValueError("space name must be nonempty")
            if basis.space is not None:
                if basis.space.name != space:
                    raise ValueError("axis space disagrees with basis space")
                space = basis.space
            else:
                space = Space(space, len(basis))
        if space is not None:
            basis.bind_space(space)
        if basis.space is None:
            raise ValueError("matrix axis requires a Basis with an explicit Space")
        labels = tuple(basis.elements)
        if len(set(labels)) != len(labels):
            raise ValueError("axis basis elements must be distinct")
        if len(labels) != basis.space.dimension:
            raise ValueError("basis size differs from space dimension")
        object.__setattr__(self, "basis", basis)
        object.__setattr__(self, "_labels", labels)

    @property
    def space(self):
        return self.basis.space

    @property
    def labels(self):
        if tuple(self.basis.elements) != self._labels:
            raise ValueError("axis basis changed; construct a new matrix context")
        return self._labels


@dataclass(frozen=True, eq=False, init=False)
class LabeledMatrix(MatrixOperations):
    """Dense real storage with optional validated network meaning.

    Generic data must be finite. Only resistance matrices allow +inf,
    separating connected components. Operations return new objects.
    """
    rows: MatrixAxis
    columns: MatrixAxis
    kind: MatrixKind
    representation: MatrixRepresentation
    tolerance: MatrixTolerance
    _values: np.ndarray = field(repr=False)
    _components: tuple = field(repr=False)

    def __init__(self, values, *, rows, columns, kind,
                 representation="generic", tolerance=None):
        if not isinstance(rows, MatrixAxis) or not isinstance(columns, MatrixAxis):
            raise TypeError("rows and columns must be MatrixAxis instances")
        if not isinstance(kind, MatrixKind):
            raise TypeError("kind must be a MatrixKind")
        representation = MatrixRepresentation(representation)
        tolerance = MatrixTolerance() if tolerance is None else tolerance
        if not isinstance(tolerance, MatrixTolerance):
            raise TypeError("tolerance must be a MatrixTolerance")
        array = np.asarray(values)
        if array.dtype.kind not in "iuf":
            raise TypeError("matrix entries must be real numeric values")
        if array.shape != (len(rows.labels), len(columns.labels)):
            raise ValueError("matrix shape must match both axis contexts")
        array = np.array(array, dtype=float, copy=True)
        if representation != "resistance" and not np.isfinite(array).all():
            raise ValueError("matrix entries must be finite")
        for key, value in (("rows", rows), ("columns", columns), ("kind", kind),
                           ("representation", representation), ("tolerance", tolerance)):
            object.__setattr__(self, key, value)
        groups = None
        if representation != "generic":
            if representation == "bipartite":
                if rows.space is columns.space:
                    raise ValueError("bipartite axes require distinct Space instances")
            else:
                self._require_same_axes()
            for label in rows.labels + columns.labels:
                if hasattr(label, "affine_weight") and label.affine_weight != 1:
                    raise ValueError("network representations require vertex labels of affine weight one")
            if representation in ("adjacency", "bipartite"):
                if kind != MatrixKind.PAIRWISE:
                    raise ValueError("adjacency and bipartite data require pairwise role")
                if np.any(array < 0):
                    raise ValueError("coupling coefficients must be nonnegative")
                if representation == "adjacency":
                    if np.any(np.diag(array) != 0):
                        raise ValueError("adjacency requires zero diagonal; loops are unsupported")
                    if num.symmetric(array, tolerance):
                        groups = num.components((array > 0) | (array.T > 0))
            elif representation == "laplacian":
                if kind not in (MatrixKind.OPERATOR, MatrixKind.FORM_COEFFICIENTS):
                    raise ValueError("laplacian requires operator or form-coefficient role")
                num.validate_laplacian(array, tolerance)
                if num.symmetric(array, tolerance):
                    # Do not delete a weak edge just because it is below a
                    # numerical rank tolerance; inversion will reject it if unresolved.
                    groups = num.components((array < 0) | (array.T < 0))
            elif representation == "greenian":
                if kind not in (MatrixKind.METRIC, MatrixKind.FORM_COEFFICIENTS):
                    raise ValueError("greenian requires metric or form-coefficient role")
                groups = num.validate_greenian(array, tolerance)
            else:
                if kind != MatrixKind.PAIRWISE:
                    raise ValueError("resistance requires pairwise role")
                groups = num.validate_resistance(array, tolerance)
        # Bytes-backed storage cannot be made writable through setflags().
        stored = np.frombuffer(array.tobytes(), dtype=float).reshape(array.shape)
        object.__setattr__(self, "_values", stored)
        object.__setattr__(self, "_components", groups)

    def _require_same_axes(self):
        if self.rows.space is not self.columns.space or self.rows.labels != self.columns.labels:
            raise ValueError("operation requires the same space and ordered basis on both axes")

    def _array(self, finite=True):
        self.rows.labels
        self.columns.labels
        if finite and not np.isfinite(self._values).all():
            raise ValueError("operation requires finite data; use componentwise network conversions")
        return self._values

    def _tol(self, atol=None, rtol=None):
        return MatrixTolerance(self.tolerance.atol if atol is None else atol,
                               self.tolerance.rtol if rtol is None else rtol)

    def _new(self, values, *, kind=None, representation="generic", rows=None, columns=None):
        matrix_class = LabeledMatrix
        if representation == "laplacian":
            from .laplacian import Laplacian
            matrix_class = Laplacian
        return matrix_class(values, rows=self.rows if rows is None else rows,
                             columns=self.columns if columns is None else columns,
                             kind=self.kind if kind is None else kind,
                             representation=representation, tolerance=self.tolerance)

    def to_numpy(self):
        """Return independent storage in the declared axis order."""
        return self._array(finite=False).copy()

    @property
    def components(self):
        """Label tuples of symmetric network components, including isolated vertices."""
        self._array(finite=False)
        if self._components is None:
            raise ValueError("components require a symmetric network representation")
        return tuple(tuple(self.rows.labels[i] for i in group) for group in self._components)

    def reordered(self, *, rows=None, columns=None):
        """Full permutations; bipartite axes may be reordered independently."""
        self._array(finite=False)
        rb = self.rows.basis.reordered(self.rows.labels if rows is None else rows)
        cb = self.columns.basis.reordered(self.columns.labels if columns is None else columns)
        ri = [self.rows.labels.index(x) for x in rb.elements]
        ci = [self.columns.labels.index(x) for x in cb.elements]
        return self._new(self._values[np.ix_(ri, ci)], rows=MatrixAxis(basis=rb),
                         columns=MatrixAxis(basis=cb), representation=self.representation)

    def change_basis(self, *, rows, columns, row_transition, column_transition):
        """Old-to-new maps: operator M_r A M_c^-1; metric M_r^-T G M_c^-1;
        coefficients M_r C M_c.T. Pairwise data has only permutations.

        Network representation is cleared after arbitrary basis changes.
        These numerical operations introduce no new formal bilinear forms.
        """
        a = self._array()
        if self.kind == MatrixKind.PAIRWISE:
            raise ValueError("pairwise data supports reordering, not linear basis change")
        maps = []
        for old, new, values in ((self.rows, rows, row_transition),
                                  (self.columns, columns, column_transition)):
            if not isinstance(new, MatrixAxis) or new.space is not old.space:
                raise ValueError("basis change must stay in the same Space")
            new.labels
            raw = np.asarray(values)
            if raw.dtype.kind not in "iuf":
                raise TypeError("transition must be real numeric")
            m = np.array(raw, dtype=float)
            if m.shape != (len(old.labels), len(old.labels)) or not np.isfinite(m).all():
                raise ValueError("transition must be finite and square in its space")
            s = np.linalg.svd(m, compute_uv=False)
            if np.count_nonzero(s > self.tolerance.threshold(num.scale(s))) != len(m):
                raise ValueError("transition must be invertible")
            maps.append(m)
        mr, mc = maps
        if self.kind == MatrixKind.FORM_COEFFICIENTS:
            result = mr @ a @ mc.T
        else:
            result = np.linalg.solve(mc.T, a.T).T
            result = mr @ result if self.kind == MatrixKind.OPERATOR else np.linalg.solve(mr.T, result)
        return self._new(result, rows=rows, columns=columns)
