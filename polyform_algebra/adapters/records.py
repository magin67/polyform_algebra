"""Versioned JSON-safe computation records; no model SDK or executable input.

Transport contract: docs/PMG_minimal_core_technical_spec_2026-11-08.md, section 13.
The dispatcher deliberately exposes a finite set of network operations only.
"""
from dataclasses import fields, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
from math import isfinite
from numbers import Integral, Rational, Real

from .. import get_version
from .. import core

DEFINITION_REVISION = '4bea8e6e52c1299510dab05e5d96dbe6ea79d796'
DEFINITION_URL = ('https://github.com/magin67/pmg-notes/blob/' + DEFINITION_REVISION
                  + '/pmg-library/en/Getting%20Started%20and%20Foundations/Foundations/Core%20specification.md')

# Explicit keywords prevent accidental exposure of new implementation parameters.
_OPERATIONS = {
    'laplacian': (core.laplacian, ()),
    'green': (core.green, ()),
    'resistance': (core.resistance, ('pairs',)),
    'spanning_tree_weight': (core.spanning_tree_weight, ()),
    'dirichlet': (core.dirichlet, ('boundary',)),
    'kron_reduce': (core.kron_reduce, ('keep',)),
    'contract': (core.contract, ('partition',)),
    'induced_subgraph': (core.induced_subgraph, ('keep',)),
    'simplex_volume_squared': (core.simplex_volume_squared, ('vertices',)),
}
OPERATIONS = tuple(_OPERATIONS)


def canonical_json(value):
    """Serialize already encoded data without NaN, Infinity or locale dependence."""
    return json.dumps(value, ensure_ascii=True, sort_keys=True, separators=(',', ':'), allow_nan=False)


def encode(value):
    """Encode supported values explicitly; never stringify arbitrary Python objects.

    Typed integers, strings and tuples preserve vertex-label distinctions.
    Fractions use decimal strings, avoiding JSON integer precision limits.
    Dataclass records are output data only, not instructions for reconstruction.
    """
    if value is core.INFINITY:
        return {'type': 'positive_infinity'}
    if value is None or type(value) is bool:
        return value
    if isinstance(value, Integral):
        return {'type': 'integer', 'value': str(int(value))}
    if isinstance(value, Rational):
        q = Fraction(value)
        return {'type': 'rational', 'numerator': str(q.numerator), 'denominator': str(q.denominator)}
    if isinstance(value, Real):
        x = float(value)
        if x == float('inf'):
            return {'type': 'positive_infinity'}
        if not isfinite(x):
            raise core.InvalidInput('NaN and negative infinity are not transport values')
        return {'type': 'float', 'value': x}
    if type(value) is str:
        return {'type': 'string', 'value': value}
    if type(value) in (tuple, list, range):
        return {'type': 'tuple' if type(value) is tuple else 'list', 'items': [encode(x) for x in value]}
    if type(value) is dict:
        items = [[encode(k), encode(v)] for k, v in value.items()]
        # Mapping order can define a contraction's output vertex order.
        return {'type': 'mapping', 'items': items}
    if is_dataclass(value) and type(value).__module__.startswith('polyform_algebra.core.'):
        return {'type': type(value).__name__, 'fields': {
            f.name: encode(getattr(value, f.name)) for f in fields(value) if not f.name.startswith('_')}}
    raise core.InvalidInput('unsupported transport value type')


def calculation_record(calculation, *, input_data, definition_revision=DEFINITION_REVISION):
    """Record a successful core calculation, including diagnostic precision policy.

    The digest covers the operation, input, arithmetic and policy, not runtime
    context UUIDs. Callers supply explicit input data, not opaque live objects.
    """
    if not isinstance(calculation, core.Calculation):
        raise core.InvalidInput('expected Calculation')
    d = calculation.diagnostics
    payload = encode({'operation': d.operation, 'input': input_data,
                      'arithmetic': d.arithmetic, 'policy': d.policy})
    return _record(d.operation, d.arithmetic, payload, calculation=calculation,
                   definition_revision=definition_revision)


def _record(operation, arithmetic, payload, *, calculation=None, error=None,
            definition_revision=DEFINITION_REVISION):
    return {
        'schema_version': '1', 'operation': operation, 'core_version': get_version(),
        'definition_revision': definition_revision,
        'definition_url': DEFINITION_URL if definition_revision == DEFINITION_REVISION else None,
        'arithmetic': arithmetic, 'input': payload,
        'input_digest': sha256(canonical_json(payload).encode('utf-8')).hexdigest() if payload is not None else None,
        'status': calculation.diagnostics.status if calculation is not None else 'error',
        'value': encode(calculation.value) if calculation is not None else None,
        'diagnostics': encode(calculation.diagnostics if calculation is not None else error.diagnostics),
        'error': None if error is None else {'code': error.code, 'message': str(error),
                                           'parameters': encode(error.parameters)},
    }


def execute_network(operation, vertices, edges, *, arithmetic='exact', policy=None, arguments=None):
    """Execute one registered operation on Python data and return a JSON-safe record.

    This is not a JSON parser or an agent server. Invalid non-encodable inputs
    have null input/digest rather than an invented representation. Resource and
    domain errors remain typed exceptions inside the core and become records here.
    """
    payload = None
    try:
        if type(operation) is not str or type(arithmetic) is not str:
            raise core.InvalidInput('operation and arithmetic must be strings')
        effective_policy = core.NumericPolicy() if arithmetic == 'numeric' and policy is None else policy
        arguments = {} if arguments is None else arguments
        payload = encode({'operation': operation, 'vertices': vertices, 'edges': edges,
                          'arithmetic': arithmetic, 'policy': effective_policy, 'arguments': arguments})
        if operation not in _OPERATIONS:
            raise core.UnsupportedOperation('operation is not registered')
        function, allowed = _OPERATIONS[operation]
        if type(arguments) is not dict or any(key not in allowed for key in arguments):
            raise core.InvalidInput('unexpected operation arguments')
        required = set(allowed) - ({'pairs'} if operation == 'resistance' else set())
        if not required.issubset(arguments):
            raise core.InvalidInput('missing operation arguments')
        net = core.network(vertices, edges, arithmetic=arithmetic, policy=effective_policy)
        result = function(net, **arguments)
        return _record(operation, arithmetic, payload, calculation=result)
    except core.PMGError as error:
        return _record(operation if type(operation) is str else None,
                       arithmetic if type(arithmetic) is str else None, payload, error=error)
