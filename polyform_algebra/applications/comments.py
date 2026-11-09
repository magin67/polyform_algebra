"""Version 1 contracts for complete employee/document/comment snapshots.

These functions validate exchange data, not the mathematical meaning of scores.
One comment record represents one event; aggregation belongs to a later scenario.
"""
from datetime import datetime
from collections.abc import Mapping
from fractions import Fraction
import hashlib
import json
import math
import re


class DataContractError(ValueError):
    """Invalid or unsupported application exchange payload."""


def canonical_json(value):
    """Stable UTF-8 JSON text; preserve array order and reject nonfinite numbers."""
    try:
        text = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)
        text.encode('utf-8')
        return text
    except (TypeError, ValueError, UnicodeError) as error:
        raise DataContractError('payload must contain finite JSON values') from error


def _copy(value):
    if isinstance(value, Mapping):
        value = dict(value)
    # Validate JSON object keys before serialization can silently coerce integers.
    def check(x):
        if type(x) is dict:
            if any(type(k) is not str for k in x):
                raise DataContractError('JSON object keys must be strings')
            for v in x.values():
                check(v)
        elif type(x) is list:
            for v in x:
                check(v)
        elif x is not None and type(x) not in (str, int, float, bool):
            raise DataContractError('unsupported JSON value')
    check(value)
    return json.loads(canonical_json(value))


def _fields(obj, required, optional=()):
    if type(obj) is not dict or not set(required) <= set(obj) or set(obj)-set(required)-set(optional):
        raise DataContractError('missing or unknown fields; expected ' + ', '.join(required))


def _text(value):
    if type(value) is not str or not value.strip():
        raise DataContractError('expected a nonempty string; identifiers are never numbers')
    return value


def _time(value):
    _text(value)
    if not re.fullmatch(r'\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(?:\.\d{1,6})?(?:Z|[+-]\d\d:\d\d)',value):
        raise DataContractError('expected ISO timestamp with seconds and UTC offset')
    try:
        return datetime.fromisoformat(value.replace('Z','+00:00'))
    except ValueError as error:
        raise DataContractError('invalid timestamp') from error


def _rows(value):
    if type(value) is not list:
        raise DataContractError('expected an array of records')
    return value


def validate_snapshot(payload):
    """Return a detached complete snapshot; reject duplicates and dangling references."""
    p = _copy(payload)
    _fields(p, ('schema_version','kind','source_system','snapshot_id','period','exported_at',
                'complete','employees','documents','comments'), ('filters',))
    if p['schema_version'] != '1' or p['kind'] != 'comment_snapshot' or p['complete'] is not True:
        raise DataContractError('only complete comment_snapshot schema 1 is supported')
    _text(p['source_system']); _text(p['snapshot_id'])
    _fields(p['period'], ('start','end'))
    start, end = _time(p['period']['start']), _time(p['period']['end'])
    if start >= end or _time(p['exported_at']) < end:
        raise DataContractError('period must be nonempty and end no later than export time')
    if 'filters' in p and type(p['filters']) is not dict:
        raise DataContractError('filters must be an object')
    ids = {}
    for table in ('employees','documents'):
        ids[table] = set()
        for row in _rows(p[table]):
            _fields(row, ('id','label'))
            identifier = _text(row['id']); _text(row['label'])
            if identifier in ids[table]:
                raise DataContractError('duplicate '+table+' identifier')
            ids[table].add(identifier)
    seen = set()
    for row in _rows(p['comments']):
        _fields(row, ('id','employee_id','document_id','occurred_at'))
        identifier = _text(row['id'])
        if identifier in seen:
            raise DataContractError('duplicate comment identifier')
        seen.add(identifier)
        for key, table in (('employee_id','employees'),('document_id','documents')):
            if _text(row[key]) not in ids[table]:
                raise DataContractError('unknown '+key)
        if not start <= _time(row['occurred_at']) < end:
            raise DataContractError('comment outside half-open period [start,end)')
    return p


def snapshot_digest(snapshot):
    """Hash the entire validated snapshot, including IDs, labels, filters and row order."""
    return hashlib.sha256(canonical_json(validate_snapshot(snapshot)).encode('utf-8')).hexdigest()


def _number(value):
    _fields(value, ('type',), ('numerator','denominator','value'))
    if value['type'] == 'rational':
        _fields(value, ('type','numerator','denominator'))
        n, d = value['numerator'], value['denominator']
        if type(n) is not str or type(d) is not str or not re.fullmatch(r'-?(0|[1-9]\d*)',n) or not re.fullmatch(r'[1-9]\d*',d):
            raise DataContractError('rational requires decimal integer strings and a positive denominator')
        try:
            return Fraction(int(n),int(d))
        except ValueError as error:
            raise DataContractError('rational exceeds integer conversion limits') from error
    if value['type'] == 'number':
        _fields(value, ('type','value'))
        x = value['value']
        if type(x) not in (int,float) or (type(x) is float and not math.isfinite(x)):
            raise DataContractError('expected finite numeric value')
        return x
    if value == {'type':'positive_infinity'}:
        return float('inf')
    raise DataContractError('unknown value type')


def validate_result(payload, snapshot):
    """Validate results against the exact source snapshot, preserving undefined values."""
    source = validate_snapshot(snapshot)
    p = _copy(payload)
    _fields(p, ('schema_version','kind','run_id','source_system','snapshot_id','snapshot_digest',
                'model','library_version','created_at','status','metrics','edges','messages'))
    if p['schema_version'] != '1' or p['kind'] != 'comment_results' or p['status'] != 'completed':
        raise DataContractError('only completed comment_results schema 1 is supported')
    for key in ('source_system','snapshot_id'):
        if p[key] != source[key]:
            raise DataContractError('result belongs to another source or snapshot')
    if p['snapshot_digest'] != snapshot_digest(source):
        raise DataContractError('snapshot digest mismatch')
    _text(p['run_id']); _text(p['library_version'])
    if _time(p['created_at']) < _time(source['exported_at']):
        raise DataContractError('result predates snapshot export')
    _fields(p['model'], ('id','version','parameters'))
    _text(p['model']['id']); _text(p['model']['version'])
    if type(p['model']['parameters']) is not dict:
        raise DataContractError('model parameters must be an object')
    ids = {row['id'] for row in source['employees']}
    seen = set()
    for row in _rows(p['metrics']):
        _fields(row, ('employee_id','metric','status','value'), ('component','normalization','reason'))
        identifier, metric = _text(row['employee_id']), _text(row['metric'])
        if identifier not in ids or (identifier,metric) in seen:
            raise DataContractError('unknown employee or duplicate metric')
        seen.add((identifier,metric))
        for key in ('component','normalization','reason'):
            if key in row:
                _text(row[key])
        if row['status'] == 'undefined':
            if row['value'] is not None or 'reason' not in row:
                raise DataContractError('undefined metric requires null value and a reason')
        elif row['status'] in ('exact','approximate'):
            _number(row['value'])
            if row['status'] == 'exact' and row['value']['type'] == 'number':
                raise DataContractError('exact metrics require rational values or typed infinity')
        else:
            raise DataContractError('unknown metric status')
    seen = set()
    for row in _rows(p['edges']):
        _fields(row, ('source_id','target_id','weight'))
        u, v = _text(row['source_id']), _text(row['target_id'])
        pair = frozenset((u,v))
        if u not in ids or v not in ids or u == v or pair in seen:
            raise DataContractError('unknown endpoint, loop or repeated undirected edge')
        seen.add(pair)
        weight = _number(row['weight'])
        if weight <= 0 or weight == float('inf'):
            raise DataContractError('edge conductance must be positive and finite')
    for message in _rows(p['messages']):
        _text(message)
    return p
