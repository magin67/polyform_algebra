"""Finite undirected network geometry, exact rational and scaled numeric paths.

PMG v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0:
Introduction / Laplacian, Green Matrix, and Effective Resistance Matrix;
Inner Product of Graph Vectors; From Lengths to Areas and Volumes.
"""
from dataclasses import dataclass
from fractions import Fraction
from math import exp, factorial, isfinite, log

import numpy as np
import sympy as sp

from .data import (
    Context, Coordinates, MatrixData, Network, NumericPolicy, Calculation,
    Diagnostics, InvalidInput, DomainError, NumericalUnresolved, UnsupportedOperation, INFINITY,
    _compatible, _label, _labels, _ordered, network,
    _context, _matrix_number,
)
from .._core.scalars import exact


@dataclass(frozen=True)
class SolveCheck:
    labels: tuple
    scale: float
    rank: int
    condition: float
    residual: float


@dataclass(frozen=True)
class PairValues:
    context: Context
    pairs: tuple
    values: tuple

    def __post_init__(self):
        _context(self.context)
        pairs = tuple(_ordered(pair) for pair in _ordered(self.pairs))
        if any(len(pair) != 2 or any(not _label(x) or x not in self.context.labels for x in pair) for pair in pairs):
            raise InvalidInput('invalid pair labels')
        if len(set(pairs)) != len(pairs):
            raise InvalidInput('duplicate ordered pair')
        values = _ordered(self.values)
        if len(values) != len(pairs):
            raise InvalidInput('pair/value size mismatch')
        object.__setattr__(self, 'pairs', pairs)
        object.__setattr__(self, 'values', tuple(_matrix_number(x,self.context.arithmetic,'resistance') for x in values))


@dataclass(frozen=True)
class CircumsphereResult:
    """Squared radius and affine coordinates in the ordered vertex context."""

    resistance_radius: Fraction
    center_coordinates: Coordinates


def _require(net):
    if not isinstance(net, Network):
        raise InvalidInput('expected Network')


def _result(net, operation, value, numerical=(), corrections=()):
    return Calculation(value, Diagnostics(operation, net.context.arithmetic, net.policy,
                       checks=('validated domain',), numerical=tuple(numerical), corrections=tuple(corrections)))


def _matrix(values, context, representation):
    return MatrixData(tuple(tuple(row) for row in values), context, context, representation)


def _finite(array):
    if not np.isfinite(array).all():
        raise NumericalUnresolved('result is outside the finite float64 range')
    return array


def _groups(net):
    index = {label: i for i, label in enumerate(net.context.labels)}
    return [tuple(index[x] for x in group) for group in net.components]


def laplacian(net):
    """M01: sum conductances on the diagonal; retain isolates and vertex order."""
    _require(net)
    n = len(net.context.labels)
    index = {x: i for i, x in enumerate(net.context.labels)}
    values = [[0 for _ in range(n)] for _ in range(n)]
    for u, v, c in net.edges:
        i, j = index[u], index[v]
        values[i][j] -= c
        values[j][i] -= c
        values[i][i] += c
        values[j][j] += c
    if net.context.arithmetic == 'numeric':
        _finite(np.array(values, dtype=float))
    return _result(net, 'laplacian', _matrix(values, net.context, 'laplacian'))


def _inverse_centered(block, labels, policy):
    """Invert only the n-1 positive modes, checking rank and forward-error budget.

    The zero mode is known from network topology, never inferred by deleting
    weak edges. Scaling avoids using conductance tolerances on inverse units.
    """
    n = len(block)
    if n == 1:
        return np.zeros((1, 1)), SolveCheck(labels, 0.0, 0, 1.0, 0.0)
    unit = float(np.max(np.abs(block)))
    if not unit or not isfinite(unit):
        raise NumericalUnresolved('component has no resolvable positive scale')
    a = block / unit
    try:
        eigenvalues, vectors = np.linalg.eigh(a)
    except np.linalg.LinAlgError as error:
        raise NumericalUnresolved('centered eigensolver did not converge') from error
    cutoff = max(policy.rtol, 100*n*np.finfo(float).eps) * max(abs(eigenvalues))
    if eigenvalues[0] < -cutoff:
        raise InvalidInput('centered matrix is not positive semidefinite')
    if abs(eigenvalues[0]) > cutoff or np.any(eigenvalues[1:] <= cutoff):
        raise NumericalUnresolved('connected component rank is unresolved')
    condition = float(eigenvalues[-1] / eigenvalues[1])
    budget = max(policy.rtol, 100*n*np.finfo(float).eps)
    if 10*n*np.finfo(float).eps*condition > budget:
        raise NumericalUnresolved('condition estimate exceeds the accuracy budget')
    inverse = (vectors[:, 1:] / eigenvalues[1:]) @ vectors[:, 1:].T
    j = np.eye(n) - np.ones((n, n))/n
    # Explicit projection restores the centroid gauge after roundoff.
    inverse = j @ inverse @ j
    residual = float(np.linalg.norm(a @ inverse-j, ord=np.inf))
    if not isfinite(residual) or residual > budget:
        raise NumericalUnresolved('centered inverse residual exceeds tolerance')
    with np.errstate(over='ignore', invalid='ignore'):
        result = _finite(inverse / unit)
    return result, SolveCheck(labels, unit, n-1, condition, residual)


def green(net):
    """M02: centered inverse per original connected component, zero across blocks."""
    _require(net)
    l = laplacian(net).value.values
    n = len(l)
    result = [[0 for _ in range(n)] for _ in range(n)]
    checks = []
    for group, labels in zip(_groups(net), net.components):
        m = len(group)
        block = [[l[i][j] for j in group] for i in group]
        if net.context.arithmetic == 'exact':
            gauge = sp.ones(m)/m
            inverse = (sp.Matrix(block)+gauge).inv()-gauge
            values = [[exact(x) for x in row] for row in inverse.tolist()]
        else:
            values, check = _inverse_centered(np.array(block), labels, net.policy)
            checks.append(check)
        for i, a in enumerate(group):
            for j, b in enumerate(group):
                result[a][b] = values[i][j]
    return _result(net, 'green', _matrix(result, net.context, 'greenian'), checks)


def resistance(net, *, pairs=None):
    """M03: resistance is a squared difference length; disconnected pairs are infinite."""
    _require(net)
    labels = net.context.labels
    index = {x: i for i, x in enumerate(labels)}
    selected = None
    if pairs is not None:
        selected = tuple(_ordered(pair) for pair in _ordered(pairs))
        for pair in selected:
            if len(pair) != 2 or any(not _label(x) or x not in index for x in pair):
                raise InvalidInput('each requested pair must contain two known labels')
        if len(set(selected)) != len(selected):
            raise InvalidInput('duplicate ordered pair')
    calc = green(net)
    g = calc.value.values
    component = {x: k for k, group in enumerate(net.components) for x in group}
    infinity = INFINITY if net.context.arithmetic == 'exact' else float('inf')

    def pair_value(u, v):
        if component[u] != component[v]:
            return infinity
        if u == v:
            return 0
        i, j = index[u], index[v]
        value = (g[i][i] + g[j][j]) - (g[i][j] + g[j][i])
        if net.context.arithmetic == 'numeric' and (not isfinite(value) or value <= 0):
            raise NumericalUnresolved('connected pair resistance is not resolvably positive')
        return value

    if selected is not None:
        value = PairValues(net.context, selected, tuple(pair_value(*p) for p in selected))
    else:
        value = _matrix([[pair_value(a, b) for b in labels] for a in labels], net.context, 'resistance')
    return _result(net, 'resistance', value, calc.diagnostics.numerical)


def circumsphere(net):
    """Return the exact circumcenter and squared radius of a resistance simplex.

    The network must be exact, nonempty and connected. The sphere is taken in
    the affine hull of its vertices. Center coordinates may be negative.
    """
    _require(net)
    if net.context.arithmetic != 'exact':
        raise UnsupportedOperation('core circumsphere requires exact arithmetic')
    if len(net.components) != 1:
        raise DomainError('circumsphere requires a nonempty connected network')

    resistances = resistance(net).value.values
    labels = net.context.labels
    index = {label: i for i, label in enumerate(labels)}
    center = [Fraction(1) for _ in labels]
    for a, b, conductance in net.edges:
        i, j = index[a], index[b]
        # Each edge contributes its spanning-tree inclusion probability to
        # both expected degrees. The center is one minus half that degree.
        half_probability = conductance * resistances[i][j] / 2
        center[i] -= half_probability
        center[j] -= half_probability

    # Equal-radius equations give R @ center = 2 * radius * ones.
    radius = sum(r * alpha for r, alpha in zip(resistances[0], center)) / 2
    result = CircumsphereResult(radius, Coordinates(net.context, center))
    return _result(net, 'circumsphere', result)


def spanning_tree_weight(net, *, componentwise=False):
    """M04: Kirchhoff cofactor, with explicit empty/disconnected conventions."""
    _require(net)
    if type(componentwise) is not bool:
        raise InvalidInput('componentwise must be bool')
    if not componentwise and len(net.components) != 1:
        return _result(net, 'spanning_tree_weight', 0)
    l = laplacian(net).value.values
    values, checks = [], []
    for group, labels in zip(_groups(net), net.components):
        if len(group) == 1:
            value = 1
        elif net.context.arithmetic == 'exact':
            value = exact(sp.Matrix([[l[i][j] for j in group[:-1]] for i in group[:-1]]).det())
        else:
            block = np.array([[l[i][j] for j in group] for i in group])
            _, check = _inverse_centered(block, labels, net.policy)
            checks.append(check)
            sign, determinant_log = np.linalg.slogdet(block[:-1, :-1]/check.scale)
            determinant_log += (len(group)-1)*log(check.scale)
            if sign <= 0 or determinant_log > log(np.finfo(float).max) or determinant_log < log(np.nextafter(0., 1.)):
                raise NumericalUnresolved('spanning-tree weight exceeds float64 range')
            value = exp(determinant_log)
        values.append((labels, value))
    return _result(net, 'spanning_tree_weight', values if componentwise else values[0][1], checks)


def network_from_laplacian(values, vertices, *, arithmetic='exact', policy=None):
    """C02: validate signs and sums before building a network; never prune an edge."""
    context = Context(vertices, arithmetic)
    data = MatrixData(values, context, context, 'laplacian')
    # Validate policy even on empty input using the regular constructor.
    empty = Network(context, (), policy)
    a = data.values
    n = len(a)
    edges, corrections = [], []
    if arithmetic == 'exact':
        if any(sum(row) != 0 for row in a) or any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
            raise InvalidInput('laplacian must be symmetric with zero row sums')
    else:
        for i in range(n):
            scale = max(map(abs, a[i]), default=0)
            if scale and abs(sum(x/scale for x in a[i])) > empty.policy.conductance_atol/scale + empty.policy.rtol*sum(abs(x/scale) for x in a[i]):
                raise InvalidInput('laplacian row sum exceeds tolerance')
        for i in range(n):
            for j in range(i):
                scale = max(abs(a[i][j]), abs(a[j][i]))
                if abs(a[i][j]-a[j][i]) > empty.policy.conductance_atol + empty.policy.rtol*scale:
                    raise InvalidInput('laplacian symmetry exceeds tolerance')
                if a[i][j] != a[j][i]:
                    corrections.append(('symmetrized', i, j))
    if any(a[i][j] > 0 for i in range(n) for j in range(n) if i != j) or any(a[i][i] < 0 for i in range(n)):
        raise InvalidInput('laplacian signs do not describe nonnegative conductances')
    for i in range(n):
        for j in range(i+1, n):
            c = -a[i][j] if arithmetic == 'exact' else -(a[i][j]/2+a[j][i]/2)
            if c:
                edges.append((context.labels[i], context.labels[j], c))
    net = Network(context, tuple(edges), policy)
    rebuilt = laplacian(net).value.values
    if arithmetic == 'numeric':
        for i in range(n):
            if rebuilt[i][i] != a[i][i]:
                corrections.append(('rebuilt diagonal', i))
    return _result(net, 'network_from_laplacian', net, corrections=corrections)


def network_from_resistance(values, vertices, *, arithmetic='exact', policy=None):
    """C03: double-center finite blocks, then validate inverse network signs.

    Numeric inverse reconstruction refuses uncertain near-zero conductances:
    choosing a sparse topology from roundoff requires additional information.
    """
    context = Context(vertices, arithmetic)
    data = MatrixData(values, context, context, 'resistance')
    empty = Network(context, (), policy)
    a = data.values
    n = len(a)
    infinite = lambda x: x is INFINITY or (arithmetic == 'numeric' and x == float('inf'))
    for i in range(n):
        if a[i][i] != 0:
            raise InvalidInput('resistance diagonal must be zero')
        for j in range(n):
            if infinite(a[i][j]) != infinite(a[j][i]):
                raise InvalidInput('infinity pattern must be symmetric')
            if not infinite(a[i][j]) and (a[i][j] < 0 or a[i][j] != a[j][i]):
                raise InvalidInput('resistance must be symmetric and nonnegative')
    unseen, groups = set(range(n)), []
    for i in range(n):
        if i not in unseen:
            continue
        group = tuple(j for j in range(n) if not infinite(a[i][j]))
        if not set(group) <= unseen or any(infinite(a[j][k]) for j in group for k in group):
            raise InvalidInput('finite resistance entries must form complete component blocks')
        unseen.difference_update(group)
        groups.append(group)
    edges, checks = [], []
    for group in groups:
        m = len(group)
        if m == 1:
            continue
        labels = tuple(context.labels[i] for i in group)
        if arithmetic == 'exact':
            r = sp.Matrix([[a[i][j] for j in group] for i in group])
            j = sp.eye(m)-sp.ones(m)/m
            g = -j*r*j/2
            if g[:-1, :-1].is_positive_definite is not True:
                raise InvalidInput('resistance block has no positive definite metric on the zero-sum space')
            gauge = sp.ones(m)/m
            l = (g+gauge).inv()-gauge
            weights = [[-exact(l[i, k]) for k in range(m)] for i in range(m)]
        else:
            r = np.array([[a[i][j] for j in group] for i in group])
            unit = float(np.max(r))
            if not unit:
                raise InvalidInput('zero resistance block cannot contain multiple vertices')
            j = np.eye(m)-np.ones((m,m))/m
            g = -j @ (r/unit) @ j/2
            l, check = _inverse_centered(g, labels, empty.policy)
            checks.append(SolveCheck(labels, check.scale*unit, check.rank, check.condition, check.residual))
            weights = _finite(-l/unit)
            uncertainty = max(empty.policy.rtol, 100*m*np.finfo(float).eps)*float(np.max(np.abs(weights)))
        for i in range(m):
            for k in range(i+1, m):
                c = weights[i][k]
                if arithmetic == 'numeric' and abs(c) <= uncertainty:
                    raise NumericalUnresolved('inverse resistance has an uncertain zero conductance; use exact data')
                if c < 0:
                    raise InvalidInput('inverse resistance has negative conductance')
                if c:
                    edges.append((labels[i], labels[k], c))
    net = Network(context, tuple(edges), policy)
    if arithmetic == 'numeric':
        restored = resistance(net)
        for group in groups:
            original = np.array([[a[i][j] for j in group] for i in group])
            recovered = np.array([[restored.value.values[i][j] for j in group] for i in group])
            unit = float(np.max(original,initial=0))
            if unit and np.max(np.abs(original/unit-recovered/unit)) > empty.policy.rtol + empty.policy.resistance_atol/unit:
                raise NumericalUnresolved('reconstructed resistance exceeds tolerance')
        checks.extend(restored.diagnostics.numerical)
    return _result(net, 'network_from_resistance', net, checks)


def _vectors(net, vectors):
    vectors = _ordered(vectors)
    if not net.context.labels:
        raise DomainError('metric vectors require a nonempty network')
    for vector in vectors:
        if not isinstance(vector, Coordinates):
            raise InvalidInput('metric arguments must be Coordinates')
        _compatible(net.context, vector.context)
        for group in _groups(net):
            entries = [vector.values[i] for i in group]
            if net.context.arithmetic == 'exact':
                valid = sum(entries) == 0
            else:
                unit = max(map(abs, entries), default=0)
                valid = not unit or abs(sum(x/unit for x in entries)) <= net.policy.coordinate_atol/unit + net.policy.rtol*sum(abs(x/unit) for x in entries)
            if not valid:
                raise DomainError('metric vectors must have zero sum on every component')
    return vectors


def gram(net, U, V=None):
    """G02: U.T Gr V with axes labeled by vector positions, not graph vertices."""
    _require(net)
    U = _vectors(net, U)
    V = U if V is None else _vectors(net, V)
    calc = green(net)
    n = len(net.context.labels)
    if net.context.arithmetic == 'exact':
        u = sp.Matrix(n, len(U), lambda i,j: U[j].values[i])
        v = sp.Matrix(n, len(V), lambda i,j: V[j].values[i])
        answer = u.T*sp.Matrix(calc.value.values)*v
        values = [[exact(answer[i,j]) for j in range(len(V))] for i in range(len(U))]
    else:
        u = np.array([x.values for x in U], dtype=float).reshape(len(U), n)
        v = np.array([x.values for x in V], dtype=float).reshape(len(V), n)
        uscale = np.max(np.abs(u),axis=1,initial=0)
        vscale = np.max(np.abs(v),axis=1,initial=0)
        unorm = u/np.where(uscale==0,1,uscale)[:,None]
        vnorm = v/np.where(vscale==0,1,vscale)[:,None]
        with np.errstate(over='ignore', invalid='ignore'):
            normalized = _finite(unorm @ calc.value.to_numpy() @ vnorm.T)
        values = []
        for i in range(len(U)):
            row = []
            for j in range(len(V)):
                coefficient = float(normalized[i,j])
                if U[i] == V[j] and uscale[i] and coefficient <= 0:
                    raise NumericalUnresolved('nonzero vector has an unresolved self-norm')
                if not coefficient or not uscale[i] or not vscale[j]:
                    value = 0.0
                else:
                    value = _scaled_value(coefficient, log(uscale[i])+log(vscale[j]))
                row.append(value)
            values.append(row)
    rows = Context(tuple(range(len(U))), net.context.arithmetic)
    columns = rows if V is U else Context(tuple(range(len(V))), net.context.arithmetic)
    return _result(net, 'gram', MatrixData(values, rows, columns, 'gram'), calc.diagnostics.numerical)


def inner_product(net, u, v):
    """G01: finite resistance inner product on the componentwise zero-sum space."""
    calc = gram(net, (u,), (v,))
    value = calc.value.values[0][0]
    corrections = []
    if net.context.arithmetic == 'numeric' and u == v and value < 0:
        # A negative self-norm is never turned positive with abs().
        raise NumericalUnresolved('negative numeric self-norm')
    return _result(net, 'inner_product', value, calc.diagnostics.numerical, corrections)


def blade_inner_product(net, U, V):
    """G03: determinant of the cross Gram matrix, including orientation signs."""
    _require(net)
    U, V = _ordered(U), _ordered(V)
    if len(U) != len(V):
        raise DomainError('exterior grades must match')
    if not U:
        return _result(net, 'blade_inner_product', 1)
    calc = gram(net, U, V)
    if net.context.arithmetic == 'exact':
        value = exact(sp.Matrix(calc.value.values).det())
    else:
        # An exactly dependent list of supplied binary floats has zero blade.
        # This algebraic check must not discard a merely small positive volume.
        dependent = any(sp.Matrix([[sp.Rational(float(x)) for x in vector.values] for vector in vectors]).rank() < len(vectors)
                        for vectors in (U,V))
        if dependent:
            value = 0.0
        else:
            a = calc.value.to_numpy()
            unit = float(np.max(np.abs(a)))
            if not unit:
                value = 0.0
            else:
                sign, determinant_log = np.linalg.slogdet(a/unit)
                if U == V and sign <= 0:
                    raise NumericalUnresolved('positive blade norm cannot be resolved')
                if sign:
                    condition = float(np.linalg.cond(a/unit))
                    if 10*len(U)*np.finfo(float).eps*condition > max(net.policy.rtol,100*len(U)*np.finfo(float).eps):
                        raise NumericalUnresolved('Gram determinant condition exceeds accuracy budget')
                value = _scaled_value(float(sign), float(determinant_log)+len(U)*log(unit)) if sign else 0.0
    return _result(net, 'blade_inner_product', value, calc.diagnostics.numerical)


def simplex_volume_squared(net, vertices):
    """G04: determinant of the edge Gram matrix divided by (k!) squared."""
    _require(net)
    vertices = _labels(vertices)
    if not vertices or any(x not in net.context.labels for x in vertices):
        raise DomainError('simplex requires distinct known vertices')
    if not any(set(vertices) <= set(group) for group in net.components):
        raise DomainError('simplex vertices must lie in one component')
    vectors = []
    for vertex in vertices[1:]:
        vectors.append(Coordinates(net.context, tuple(int(x == vertex)-int(x == vertices[0]) for x in net.context.labels)))
    calc = blade_inner_product(net, vectors, vectors)
    divisor = factorial(len(vectors))**2
    if net.context.arithmetic == 'exact':
        value = Fraction(calc.value, divisor)
    else:
        value = _scaled_value(calc.value,-log(divisor)) if calc.value else 0.0
    if value < 0:
        raise NumericalUnresolved('negative squared volume')
    return _result(net, 'simplex_volume_squared', value, calc.diagnostics.numerical)


def _scaled_value(coefficient, log_scale):
    """Restore a nonzero product without intermediate underflow or overflow."""
    exponent = log(abs(coefficient))+log_scale
    if not isfinite(exponent) or exponent > log(np.finfo(float).max) or exponent < log(np.nextafter(0.,1.)):
        raise NumericalUnresolved('nonzero metric value is outside float64 range')
    value = exp(exponent)
    if not isfinite(value) or value == 0:
        raise NumericalUnresolved('nonzero metric value is not representable')
    return value if coefficient > 0 else -value
