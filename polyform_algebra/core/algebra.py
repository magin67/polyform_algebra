"""Exact restricted quadratic algebra (A01-A12), with bounded expansions.

PMG v5 c909ddf8a590c69ed0a82fe48ca496e88fb707e0: Core Specification,
sections Exterior objects, Forms, Spanning-tree coefficient and Metric identity.
The public grammar contains only squares of products of atoms and boundaries.
Ordered tensor coefficient pairs below are PRIVATE comparison/evaluation data,
not a public constructor for general bilinear forms. See CORE_ALGEBRA_API.md.
"""
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from numbers import Rational, Real

import sympy as sp

from .data import (
    Context, Calculation, Diagnostics, InvalidInput, DomainError,
    UnsupportedOperation, ResourceLimit, _context, _compatible, _number, _ordered,
)
from .network_math import _require, green


@dataclass(frozen=True)
class AlgebraLimits:
    """Bound each sparse intermediate and each Cartesian expansion before allocation."""
    max_terms: int = 100000

    def __post_init__(self):
        if type(self.max_terms) is not int or self.max_terms < 1:
            raise InvalidInput('max_terms must be a positive integer')


def _limits(limits=None, *objects):
    if limits is None:
        limits = AlgebraLimits(min((x.limits.max_terms for x in objects), default=100000))
    if not isinstance(limits, AlgebraLimits):
        raise InvalidInput('expected AlgebraLimits')
    return limits


def _check(size, limits):
    if size > limits.max_terms:
        raise ResourceLimit('formal expansion exceeds max_terms',
                            parameters=(('max_terms', limits.max_terms), ('requested', size)))


def _exact_context(context):
    _context(context)
    if context.arithmetic != 'exact':
        raise UnsupportedOperation('formal algebra is available only in exact arithmetic')


def _rational(value):
    if not isinstance(value, (Rational, Real)):
        raise UnsupportedOperation('formal coefficients must be rational; symbolic coefficients are unsupported')
    return _number(value, 'exact')


def _accumulate(data, key, coefficient, limits):
    value = data.get(key, 0) + coefficient
    if value:
        # Check before insertion, including terms that could cancel later.
        _check(len(data) + int(key not in data), limits)
        data[key] = value
    else:
        data.pop(key, None)


def _indices(context, indices):
    indices = _ordered(indices)
    if any(type(i) is not int or not 0 <= i < len(context.labels) for i in indices):
        raise InvalidInput('atom indices must be integer positions in the context')
    return indices


def _orientation(indices):
    if len(set(indices)) != len(indices):
        return (), 0
    inversions = sum(a > b for i, a in enumerate(indices) for b in indices[i+1:])
    return tuple(sorted(indices)), (-1 if inversions % 2 else 1)


def _items(values):
    values = values.items() if isinstance(values, Mapping) else _ordered(values)
    for entry in values:
        entry = _ordered(entry)
        if len(entry) != 2:
            raise InvalidInput('expected pairs of term and coefficient')
        yield entry


def _wedge(left, right, limits):
    _check(len(left) * len(right), limits)
    result = {}
    for a, x in left.items():
        for b, y in right.items():
            key, sign = _orientation(a + b)
            if sign:
                _accumulate(result, key, sign * x * y, limits)
    return result


def _boundary(values, context, limits):
    _check(sum(len(blade) for blade in values), limits)
    result = {}
    for blade, coefficient in values.items():
        for r, i in enumerate(blade):
            _accumulate(result, blade[:r] + blade[r+1:],
                        (-1)**r * context.epsilon[i] * coefficient, limits)
    return result


class _Linear:
    """Exact arithmetic returns data objects; named operations return Calculation."""
    def __add__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        _compatible(self.context, other.context)
        limits = _limits(None, self, other)
        return type(self)(self.context, self._entries + other._entries, limits=limits)

    def __neg__(self):
        return self * -1

    def __sub__(self, other):
        if type(other) is not type(self):
            return NotImplemented
        return self + (-other)

    def __mul__(self, scalar):
        scalar = _rational(scalar)
        return type(self)(self.context, tuple((key, scalar*c) for key, c in self._entries), limits=self.limits)

    __rmul__ = __mul__

    def __truediv__(self, scalar):
        scalar = _rational(scalar)
        if not scalar:
            raise InvalidInput('division by zero')
        return self * (1 / scalar)


@dataclass(frozen=True)
class Exterior(_Linear):
    context: Context
    coefficients: tuple
    limits: AlgebraLimits = field(default_factory=AlgebraLimits, compare=False, repr=False)

    def __post_init__(self):
        _exact_context(self.context)
        limits = _limits(self.limits)
        result = {}
        for blade, coefficient in _items(self.coefficients):
            blade = _indices(self.context, blade)
            coefficient = _rational(coefficient)
            key, sign = _orientation(blade)
            if sign:
                _accumulate(result, key, sign * coefficient, limits)
        object.__setattr__(self, 'coefficients', tuple(sorted(result.items())))

    @property
    def _entries(self):
        return self.coefficients


exterior = Exterior


def _factor_grade(factor):
    kind, indices = factor
    return len(indices) - (kind == 'boundary')


@dataclass(frozen=True)
class SimpleTerm:
    context: Context
    factors: tuple = ()
    coefficient: Fraction = Fraction(1)
    limits: AlgebraLimits = field(default_factory=AlgebraLimits, compare=False, repr=False)

    def __post_init__(self):
        _exact_context(self.context)
        limits = _limits(self.limits)
        coefficient = _rational(self.coefficient)
        factors = []
        for factor in _ordered(self.factors):
            if isinstance(factor, SimpleTerm):
                _compatible(self.context, factor.context)
                coefficient *= factor.coefficient
                factors.extend(factor.factors)
                _check(len(factors), limits)
                continue
            factor = _ordered(factor)
            if len(factor) != 2 or factor[0] not in ('atom', 'boundary'):
                raise UnsupportedOperation('SimpleTerm factors must be atoms, boundaries or SimpleTerms')
            kind, indices = factor
            indices = _indices(self.context, indices)
            if kind == 'atom' and len(indices) != 1:
                raise InvalidInput('an atom requires exactly one index')
            if kind == 'boundary':
                if not indices or len(set(indices)) != len(indices):
                    raise InvalidInput('compact boundary requires distinct nonempty indices')
                if any(self.context.epsilon[i] != 1 for i in indices):
                    raise UnsupportedOperation('compact boundary requires epsilon=1; use boundary(Exterior)')
            indices, sign = _orientation(indices)
            coefficient *= sign
            if kind == 'boundary' and len(indices) == 1:
                continue  # The boundary of a unit-weight atom is the scalar 1.
            factors.append((kind, indices))
            _check(len(factors), limits)
        # Sorting factors carries the graded exterior sign, not a plain sign
        # per swap. Squaring later removes this sign from a quadratic form.
        for i, a in enumerate(factors):
            for b in factors[i+1:]:
                if a > b:
                    coefficient *= (-1)**(_factor_grade(a) * _factor_grade(b))
        object.__setattr__(self, 'factors', tuple(sorted(factors)))
        object.__setattr__(self, 'coefficient', coefficient)

    @property
    def grade(self):
        return sum(_factor_grade(f) for f in self.factors)

    def __neg__(self):
        return self * -1

    def __mul__(self, scalar):
        return SimpleTerm(self.context, self.factors, self.coefficient * _rational(scalar), self.limits)

    __rmul__ = __mul__

    def as_exterior(self, *, limits=None):
        """Expand the supported grammar; this does not make arbitrary sums squareable."""
        limits = _limits(limits, self)
        return Exterior(self.context, tuple(_expand(self, limits).items()), limits)


simple_term = SimpleTerm


def Atom(context, index, *, limits=None):
    """A basis atom by its zero-based position, not by its display label."""
    return SimpleTerm(context, (('atom', (index,)),), limits=_limits(limits))


def Boundary(context, indices, *, limits=None):
    """A compact simplex boundary; orientation follows the supplied index order."""
    return SimpleTerm(context, (('boundary', indices),), limits=_limits(limits))


def Unit(context, *, limits=None):
    """The scalar exterior unit; square(Unit(context)) gives the form unit e."""
    return SimpleTerm(context, limits=_limits(limits))


def _expand(term, limits):
    if not term.coefficient or term.grade > len(term.context.labels):
        return {}
    result = {(): term.coefficient}
    for kind, indices in term.factors:
        factor = {indices: Fraction(1)}
        if kind == 'boundary':
            factor = _boundary(factor, term.context, limits)
        result = _wedge(result, factor, limits)
        if not result:
            break
    return result


@dataclass(frozen=True)
class Form(_Linear):
    context: Context
    terms: tuple
    limits: AlgebraLimits = field(default_factory=AlgebraLimits, compare=False, repr=False)

    def __post_init__(self):
        _exact_context(self.context)
        limits = _limits(self.limits)
        result = {}
        for term, coefficient in _items(self.terms):
            if not isinstance(term, SimpleTerm):
                raise UnsupportedOperation('form terms must be squares of SimpleTerm, not arbitrary Exterior objects')
            _compatible(self.context, term.context)
            coefficient = _rational(coefficient) * term.coefficient**2
            if coefficient and _expand(term, limits):
                key = SimpleTerm(self.context, term.factors, limits=limits)
                _accumulate(result, key, coefficient, limits)
        object.__setattr__(self, 'terms', tuple(sorted(result.items(), key=lambda item: item[0].factors)))

    @property
    def _entries(self):
        return self.terms


form = Form


def _answer(value, operation, limits):
    return Calculation(value, Diagnostics(operation, 'exact',
        checks=('validated exact context and grammar', 'bounded formal expansion',
                'max_terms={}'.format(limits.max_terms)),
        numerical=(), corrections=()))


def _external(obj, limits):
    if isinstance(obj, SimpleTerm):
        return obj.as_exterior(limits=limits)
    if not isinstance(obj, Exterior):
        raise InvalidInput('expected Exterior or SimpleTerm, not a form')
    _check(len(obj.coefficients), limits)
    return obj


def wedge(x, y, *, limits=None):
    """A01: bilinear oriented exterior product; repeated indices annihilate a blade."""
    limits = _limits(limits, *(o for o in (x, y) if isinstance(o, (Exterior, SimpleTerm))))
    x, y = _external(x, limits), _external(y, limits)
    limits = _limits(limits, x, y)
    _compatible(x.context, y.context)
    values = _wedge(dict(x.coefficients), dict(y.coefficients), limits)
    return _answer(Exterior(x.context, tuple(values.items()), limits), 'wedge', limits)


def boundary(x, *, limits=None):
    """A02: contraction by epsilon, including non-unit rational affine weights."""
    limits = _limits(limits, *((x,) if isinstance(x, (Exterior, SimpleTerm)) else ()))
    x = _external(x, limits)
    values = _boundary(dict(x.coefficients), x.context, limits)
    return _answer(Exterior(x.context, tuple(values.items()), limits), 'boundary', limits)


def grade_part(x, k, *, limits=None):
    """A03: select by exterior degree, never by number of vertices in the support."""
    if not isinstance(x, (Exterior, Form)):
        raise InvalidInput('grade_part requires Exterior or Form')
    if type(k) is not int or k < 0:
        raise DomainError('grade must be a nonnegative integer')
    limits = _limits(limits, x)
    entries = tuple((key, c) for key, c in x._entries
                    if (len(key) if isinstance(x, Exterior) else key.grade) == k)
    return _answer(type(x)(x.context, entries, limits), 'grade_part', limits)


def square(term, *, limits=None):
    """A04: a formal square of a supported simple structure, not its self-wedge."""
    if not isinstance(term, SimpleTerm):
        raise UnsupportedOperation('square accepts only SimpleTerm; arbitrary sums are unsupported')
    limits = _limits(limits, term)
    return _answer(Form(term.context, ((term, 1),), limits), 'square', limits)


def _forms(f, g=None):
    if not isinstance(f, Form) or (g is not None and not isinstance(g, Form)):
        raise InvalidInput('expected Form')
    if g is not None:
        _compatible(f.context, g.context)


def _product(f, g, limits):
    _forms(f, g)
    _check(len(f.terms) * len(g.terms), limits)
    result = {}
    for x, a in f.terms:
        for y, b in g.terms:
            term = SimpleTerm(f.context, (x, y), limits=limits)
            if _expand(term, limits):
                # Both tensor factors get the same orientation, hence sign^2.
                coefficient = a * b * term.coefficient**2
                term = SimpleTerm(f.context, term.factors, limits=limits)
                _accumulate(result, term, coefficient, limits)
    return Form(f.context, tuple(result.items()), limits)


def form_product(f, g, *, limits=None):
    """A05: multiply supported squares and extend bilinearly; forms commute."""
    _forms(f, g)
    limits = _limits(limits, f, g)
    return _answer(_product(f, g, limits), 'form_product', limits)


def _pairs(f, limits):
    """Private full ordered tensor coefficients; off-diagonal entries occur twice."""
    _forms(f)
    _check(len(f.terms), limits)
    result = {}
    for term, c in f.terms:
        expansion = _expand(term, limits)
        _check(len(expansion)**2, limits)
        for a, x in expansion.items():
            for b, y in expansion.items():
                _accumulate(result, (a, b), c*x*y, limits)
    return result


def equivalent(f, g, *, limits=None):
    """A06: exact coordinate-pair comparison, including scalars and top grade."""
    _forms(f, g)
    limits = _limits(limits, f, g)
    return _answer(_pairs(f, limits) == _pairs(g, limits), 'equivalent', limits)


def polar_form(u, v, *, limits=None):
    """A07: {u,v}, without 1/2, in atom and pair-boundary squares only."""
    limits = _limits(limits, *(o for o in (u, v) if isinstance(o, (Exterior, SimpleTerm))))
    u, v = _external(u, limits), _external(v, limits)
    _compatible(u.context, v.context)
    if any(len(blade) != 1 for x in (u, v) for blade, _ in x.coefficients):
        raise UnsupportedOperation('polarization supports only grade one (or zero)')
    if any(e != 1 for e in u.context.epsilon):
        raise UnsupportedOperation('polarization requires a vertex basis with epsilon=1')
    a, b = dict(u.coefficients), dict(v.coefficients)
    sa, sb = sum(a.values()), sum(b.values())
    terms = {}
    for i in range(len(u.context.labels)):
        ai, bi = a.get((i,), 0), b.get((i,), 0)
        _accumulate(terms, Atom(u.context, i, limits=limits), sa*bi + sb*ai, limits)
        for j in range(i+1, len(u.context.labels)):
            c = -(ai*b.get((j,), 0) + bi*a.get((j,), 0))
            if c:
                _accumulate(terms, Boundary(u.context, (i, j), limits=limits), c, limits)
    return _answer(Form(u.context, tuple(terms.items()), limits), 'polar_form', limits)


def _closed_pairs(f, limits, *, vertex_space=False):
    if vertex_space and (not f.context.labels or any(e != 1 for e in f.context.epsilon)):
        raise DomainError('requires a nonempty vertex space with epsilon=1')
    pairs = _pairs(f, limits)
    # Test contraction in each tensor slot AFTER summing all terms. This
    # accepts cancellation of non-boundary terms, not just syntactic origins.
    for side in (0, 1):
        contracted = {}
        _check(sum(len(pair[side]) for pair in pairs), limits)
        for pair, c in pairs.items():
            blade = pair[side]
            for r, i in enumerate(blade):
                reduced = blade[:r] + blade[r+1:]
                key = (reduced, pair[1]) if side == 0 else (pair[0], reduced)
                _accumulate(contracted, key, (-1)**r*f.context.epsilon[i]*c, limits)
        if contracted:
            raise DomainError('form does not belong to the boundary subalgebra')
    return pairs


def form_exp(p, *, limits=None):
    """A08: exact finite exponential via commuting decomposable square-zero terms.

    Every supported positive-grade SimpleTerm is decomposable, so its square
    as a form is nilpotent of index two. Thus exp(sum c*q)=product(e+c*q),
    which includes the factorials of the power series without floating division.
    """
    _forms(p)
    limits = _limits(limits, p)
    _check(len(p.terms), limits)
    if any(term.grade == 0 for term, _ in p.terms):
        raise UnsupportedOperation('nonzero scalar exponential is outside the rational formal domain')
    unit = Form(p.context, ((Unit(p.context, limits=limits), 1),), limits)
    result = unit
    for term, coefficient in p.terms:
        factor = Form(p.context, unit.terms + ((term, coefficient),), limits)
        result = _product(result, factor, limits)
    return _answer(result, 'form_exp', limits)


def laplacian_form(net, *, limits=None):
    """A09: sum c_ij [(ij)]^2, retaining every vertex, including isolates."""
    _require(net)
    _exact_context(net.context)
    limits = _limits(limits)
    _check(len(net.edges), limits)
    index = {x: i for i, x in enumerate(net.context.labels)}
    terms = tuple((Boundary(net.context, (index[u], index[v]), limits=limits), c) for u, v, c in net.edges)
    return _answer(Form(net.context, terms, limits), 'laplacian_form', limits)


def spanning_coefficient(f, *, limits=None):
    """A10: extract fixed T_n, not the highest nonzero term of this expression."""
    _forms(f)
    limits = _limits(limits, f)
    pairs = _closed_pairs(f, limits, vertex_space=True)
    # In partial(a_0...a_n-1), the blade omitting a_0 has coefficient +1.
    # Its diagonal tensor entry therefore reads the unique top coefficient.
    top = tuple(range(1, len(f.context.labels)))
    return _answer(pairs.get((top, top), Fraction(0)), 'spanning_coefficient', limits)


def potential(m, f, *, limits=None):
    """A11: tau(M F); well-defined even when tau(M)=0, without normalization."""
    _forms(m, f)
    limits = _limits(limits, m, f)
    _closed_pairs(m, limits, vertex_space=True)
    _closed_pairs(f, limits, vertex_space=True)
    value = spanning_coefficient(_product(m, f, limits), limits=limits).value
    return _answer(value, 'potential', limits)


def evaluate_form(net, f, *, limits=None):
    """A12: contract boundary tensor coefficients with exact Greenian minors.

    The centered Greenian is an extension of the metric from W to V. Since
    the complete form has been checked to lie on W in both slots, its minors
    give the intrinsic exterior metric even if individual terms are not on W.
    No exponential is built, and evaluations of factors are not multiplied.
    """
    _require(net)
    _forms(f)
    _compatible(net.context, f.context)
    limits = _limits(limits, f)
    if not net.context.labels or len(net.components) != 1:
        raise DomainError('metric evaluation requires a nonempty connected network')
    pairs = _closed_pairs(f, limits, vertex_space=True)
    g = green(net).value.values
    value = Fraction(0)
    for (a, b), c in pairs.items():
        determinant = sp.Matrix([[g[i][j] for j in b] for i in a]).det() if a else 1
        value += c * _rational(determinant)
    return _answer(value, 'evaluate_form', limits)
