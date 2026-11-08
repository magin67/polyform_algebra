"""Immutable contexts and data, with explicit arithmetic and ordered axes."""
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from math import isfinite
from numbers import Rational, Real
from uuid import uuid4

from .._core.scalars import exact


class PMGError(ValueError):
    """Base error with a stable machine-readable code."""
    code = 'PMG_ERROR'

    def __init__(self, message, *, parameters=(), diagnostics=None):
        super().__init__(message)
        self.parameters = tuple(parameters)
        self.diagnostics = diagnostics


class InvalidInput(PMGError):
    code = 'INVALID_INPUT'


class ContextMismatch(PMGError):
    code = 'CONTEXT_MISMATCH'


class ArithmeticMismatch(PMGError):
    code = 'ARITHMETIC_MISMATCH'


class NumericalUnresolved(PMGError):
    code = 'NUMERICAL_UNRESOLVED'


def _ordered(values):
    if not isinstance(values, Sequence) or isinstance(values, (str, bytes)):
        raise InvalidInput('expected an ordered sequence (list or tuple)')
    return tuple(values)


def _label(value):
    return type(value) in (str, int) or (type(value) is tuple and all(_label(x) for x in value))


def _labels(values):
    values = _ordered(values)
    if not all(_label(x) for x in values) or len(set(values)) != len(values):
        raise InvalidInput('labels must be unique immutable strings, integers or tuples; bool is forbidden')
    return values


def _number(value, arithmetic):
    if isinstance(value, bool):
        raise InvalidInput('bool is not a coefficient')
    if arithmetic == 'exact':
        if not isinstance(value, Rational):
            raise InvalidInput('exact arithmetic requires rational numbers, not floats or symbols')
        return exact(value)
    if not isinstance(value, Real):
        raise InvalidInput('numeric arithmetic requires finite real numbers')
    try:
        result = float(value)
    except (OverflowError, ValueError):
        raise NumericalUnresolved('number is outside float64 range') from None
    if not isfinite(result):
        raise InvalidInput('nonfinite coefficient')
    if result == 0 and value != 0:
        raise NumericalUnresolved('nonzero coefficient would underflow')
    return result


@dataclass(frozen=True)
class NumericPolicy:
    rtol: float = 1e-10
    conductance_atol: float = 0.0
    resistance_atol: float = 0.0
    coordinate_atol: float = 0.0

    def __post_init__(self):
        for name in ('rtol', 'conductance_atol', 'resistance_atol', 'coordinate_atol'):
            value = _number(getattr(self, name), 'numeric')
            if value < 0:
                raise InvalidInput('tolerances must be nonnegative')
            object.__setattr__(self, name, value)


@dataclass(frozen=True)
class Context:
    labels: tuple
    arithmetic: str = 'exact'
    epsilon: tuple = None
    identity: str = field(default_factory=lambda: uuid4().hex, init=False)

    def __post_init__(self):
        if self.arithmetic not in ('exact', 'numeric'):
            raise InvalidInput('arithmetic must be exact or numeric')
        labels = _labels(self.labels)
        weights = (1,) * len(labels) if self.epsilon is None else _ordered(self.epsilon)
        if len(weights) != len(labels):
            raise InvalidInput('epsilon size must match labels')
        object.__setattr__(self, 'labels', labels)
        object.__setattr__(self, 'epsilon', tuple(_number(x, self.arithmetic) for x in weights))

    def _derived(self, labels=None, arithmetic=None):
        labels = self.labels if labels is None else labels
        weights = dict(zip(self.labels, self.epsilon))
        result = Context(labels, arithmetic or self.arithmetic, tuple(weights[x] for x in labels))
        object.__setattr__(result, 'identity', self.identity)
        return result


def _context(value):
    if not isinstance(value, Context):
        raise InvalidInput('expected Context')


def _compatible(left, right):
    if left.arithmetic != right.arithmetic:
        raise ArithmeticMismatch('explicit to_numeric conversion is required')
    if left != right:
        raise ContextMismatch('space identity, order and epsilon must agree')


@dataclass(frozen=True)
class Coordinates:
    context: Context
    values: tuple

    def __post_init__(self):
        _context(self.context)
        values = self.values
        if isinstance(values, Mapping):
            if not all(_label(k) for k in values) or set(values) != set(self.context.labels):
                raise InvalidInput('coordinate mapping must contain exactly the context labels')
            values = tuple(values[x] for x in self.context.labels)
        else:
            values = _ordered(values)
        if len(values) != len(self.context.labels):
            raise InvalidInput('coordinate size mismatch')
        object.__setattr__(self, 'values', tuple(_number(x, self.context.arithmetic) for x in values))

    def __add__(self, other):
        if not isinstance(other, Coordinates):
            return NotImplemented
        _compatible(self.context, other.context)
        return Coordinates(self.context, tuple(a+b for a, b in zip(self.values, other.values)))

    def __neg__(self):
        return Coordinates(self.context, tuple(-x for x in self.values))

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, scalar):
        scalar = _number(scalar, self.context.arithmetic)
        return Coordinates(self.context, tuple(x * scalar for x in self.values))

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = _number(scalar, self.context.arithmetic)
        if scalar == 0:
            raise InvalidInput('division by zero')
        return Coordinates(self.context, tuple(x / scalar for x in self.values))


coordinates = Coordinates


@dataclass(frozen=True)
class MatrixData:
    values: tuple
    rows: Context
    columns: Context
    representation: str

    def __post_init__(self):
        _context(self.rows)
        _context(self.columns)
        if self.rows.arithmetic != self.columns.arithmetic:
            raise ArithmeticMismatch('matrix axes must use the same arithmetic')
        if self.representation not in ('laplacian', 'greenian', 'resistance', 'gram'):
            raise InvalidInput('unknown matrix representation')
        if self.representation != 'gram':
            _compatible(self.rows, self.columns)
        values = tuple(_ordered(row) for row in _ordered(self.values))
        if len(values) != len(self.rows.labels) or any(len(row) != len(self.columns.labels) for row in values):
            raise InvalidInput('matrix dimensions must match axes')
        object.__setattr__(self, 'values', tuple(tuple(_number(x, self.rows.arithmetic) for x in row) for row in values))

    def to_numpy(self):
        """Return an independent array; exact coefficients retain object dtype."""
        import numpy as np
        dtype = object if self.rows.arithmetic == 'exact' else float
        return np.array(self.values, dtype=dtype).reshape(len(self.rows.labels), len(self.columns.labels))


@dataclass(frozen=True)
class Network:
    context: Context
    edges: tuple
    policy: NumericPolicy = None
    components: tuple = field(init=False)

    def __post_init__(self):
        _context(self.context)
        if any(x != 1 for x in self.context.epsilon):
            raise InvalidInput('network vertices require epsilon=1')
        if self.policy is not None and not isinstance(self.policy, NumericPolicy):
            raise InvalidInput('expected NumericPolicy')
        if self.context.arithmetic == 'exact' and self.policy is not None:
            raise InvalidInput('numeric policy does not apply to exact networks')
        labels = self.context.labels
        index = {x: i for i, x in enumerate(labels)}
        prepared = {}
        for edge in _ordered(self.edges):
            edge = _ordered(edge)
            if len(edge) != 3:
                raise InvalidInput('edges must be triples (u, v, conductance)')
            u, v, c = edge
            if not _label(u) or not _label(v) or u not in index or v not in index or u == v:
                raise InvalidInput('edge endpoints must be distinct context labels')
            c = _number(c, self.context.arithmetic)
            if c < 0:
                raise InvalidInput('negative conductance')
            pair = tuple(sorted((index[u], index[v])))
            if pair in prepared:
                raise InvalidInput('Network data requires aggregated pairs; use network()')
            prepared[pair] = c
        edges = tuple((labels[i], labels[j], c) for (i, j), c in sorted(prepared.items()) if c != 0)
        neighbors = {x: set() for x in labels}
        for u, v, _ in edges:
            neighbors[u].add(v)
            neighbors[v].add(u)
        unseen, groups = set(labels), []
        for start in labels:
            if start not in unseen:
                continue
            pending, found = [start], set()
            unseen.remove(start)
            while pending:
                u = pending.pop()
                found.add(u)
                for v in neighbors[u] & unseen:
                    unseen.remove(v)
                    pending.append(v)
            groups.append(tuple(x for x in labels if x in found))
        object.__setattr__(self, 'edges', edges)
        object.__setattr__(self, 'components', tuple(groups))
        if self.context.arithmetic == 'numeric' and self.policy is None:
            object.__setattr__(self, 'policy', NumericPolicy())


def network(vertices, edges, *, arithmetic='exact', parallel='sum', policy=None):
    """Aggregate undirected conductances, preserving every declared vertex."""
    context = Context(vertices, arithmetic)
    if parallel not in ('sum', 'reject'):
        raise InvalidInput('parallel must be sum or reject')
    prepared = {}
    for edge in _ordered(edges):
        # Reuse direct-container validation before aggregation (including zero edges).
        validated = Network(context, (edge,), policy)
        u, v, c = edge
        key = frozenset((u, v))
        if key in prepared and parallel == 'reject':
            raise InvalidInput('parallel edge rejected')
        previous = prepared.get(key, (u, v, 0))
        weight = validated.edges[0][2] if validated.edges else _number(0, arithmetic)
        prepared[key] = (previous[0], previous[1], _number(previous[2] + weight, arithmetic))
    return Network(context, tuple(prepared.values()), policy)


@dataclass(frozen=True)
class Diagnostics:
    operation: str
    arithmetic: str
    policy: NumericPolicy = None
    contract_version: str = '1'
    checks: tuple = ()

    def __post_init__(self):
        if self.arithmetic not in ('exact', 'numeric'):
            raise InvalidInput('invalid diagnostic arithmetic')
        if self.policy is not None and not isinstance(self.policy, NumericPolicy):
            raise InvalidInput('invalid diagnostic policy')
        checks = _ordered(self.checks)
        if not all(isinstance(x, str) for x in checks):
            raise InvalidInput('diagnostic checks must be strings')
        object.__setattr__(self, 'checks', checks)

    @property
    def status(self):
        return 'exact' if self.arithmetic == 'exact' else 'approximate'


@dataclass(frozen=True)
class Calculation:
    value: object
    diagnostics: Diagnostics


def _calculation(value, operation, arithmetic, policy=None):
    return Calculation(value, Diagnostics(operation, arithmetic, policy, checks=('validated inputs',)))


def reorder(obj, vertices):
    """Permute data and axes together; never identify independent spaces."""
    labels = _labels(vertices)
    context = obj.rows if isinstance(obj, MatrixData) else getattr(obj, 'context', None)
    _context(context)
    if set(labels) != set(context.labels):
        raise ContextMismatch('reorder requires a full permutation')
    derived = context._derived(labels)
    positions = [context.labels.index(x) for x in labels]
    if isinstance(obj, Coordinates):
        result = Coordinates(derived, tuple(obj.values[i] for i in positions))
    elif isinstance(obj, Network):
        result = Network(derived, obj.edges, obj.policy)
    elif isinstance(obj, MatrixData):
        _compatible(obj.rows, obj.columns)
        result = MatrixData(tuple(tuple(obj.values[i][j] for j in positions) for i in positions), derived, derived, obj.representation)
    else:
        raise InvalidInput('unsupported reorder object')
    return _calculation(result, 'reorder', context.arithmetic, getattr(obj, 'policy', None))


def to_numeric(obj, *, policy=None):
    """Explicit finite conversion; reject overflow and nonzero underflow."""
    policy = NumericPolicy() if policy is None else policy
    if not isinstance(policy, NumericPolicy):
        raise InvalidInput('expected NumericPolicy')
    if isinstance(obj, Coordinates):
        result = Coordinates(obj.context._derived(arithmetic='numeric'), obj.values)
    elif isinstance(obj, Network):
        result = Network(obj.context._derived(arithmetic='numeric'), obj.edges, policy)
    elif isinstance(obj, MatrixData):
        result = MatrixData(obj.values, obj.rows._derived(arithmetic='numeric'), obj.columns._derived(arithmetic='numeric'), obj.representation)
    else:
        raise InvalidInput('to_numeric requires network, coordinates or matrix data')
    return _calculation(result, 'to_numeric', 'numeric', policy)
