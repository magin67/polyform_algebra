"""Delivery boundaries and transport contract, technical specification section 13.

Mathematical fixtures use the PMG Introduction weighted chain (c=1,2)
and the pinned v5 sources in tests/fixtures/reference_corpus_v1.json.
"""
from fractions import Fraction as Q
import json
from pathlib import Path
import subprocess
import sys

import pytest

from polyform_algebra import core
from polyform_algebra.adapters.records import encode, execute_network, calculation_record, canonical_json


def test_core_import_without_optional_or_thematic_modules():
    """Core operations must work even when optional and legacy modules cannot import."""
    script = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in ('networkx','pandas','scipy') or fullname.startswith((
            'polyform_algebra.matrix','polyform_algebra.objects','polyform_algebra.combinations',
            'polyform_algebra.adapters','polyform_algebra._core.basis','polyform_algebra._core.space')):
            raise AssertionError('forbidden dependency: '+fullname)
sys.meta_path.insert(0, Block())
from polyform_algebra.core import network, resistance, laplacian_form, spanning_coefficient, form_exp
from fractions import Fraction
for mode in ('exact','numeric'):
    net=network((0,1,2),((0,1,1),(1,2,2)),arithmetic=mode)
    assert abs(float(resistance(net).value.values[0][2])-1.5)<1e-12
net=network((0,1,2),((0,1,1),(1,2,2)))
assert spanning_coefficient(form_exp(laplacian_form(net).value).value).value == 2
'''
    result = subprocess.run([sys.executable, '-c', script], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_lazy_exports_preserve_class_identity():
    """Historic imports and module-qualified classes must remain identical after lazy loading."""
    import polyform_algebra as p
    from polyform_algebra.objects.element import Point
    from polyform_algebra.combinations.polyform import Polyform
    assert p.Point is Point and p.Polyform is Polyform
    namespace = {}
    exec('from polyform_algebra import *', namespace)
    assert all(name in namespace for name in p.__all__)
    assert set(p.__all__).issubset(dir(p))


def test_transport_rationals_infinity_and_typed_labels():
    """Large rationals are decimal strings; label 1 and label '1' stay distinct."""
    q = Q(10**80+1, 10**71+3)
    assert encode(q) == {'type':'rational','numerator':str(q.numerator),'denominator':str(q.denominator)}
    assert encode(core.INFINITY) == encode(float('inf')) == {'type':'positive_infinity'}
    assert encode(1) != encode('1') != encode((1,))
    json.loads(canonical_json(encode([q, core.INFINITY, (1, '1')])))


@pytest.mark.parametrize('value', [float('nan'), -float('inf'), object(), {1, 2}])
def test_transport_rejects_unsupported_values(value):
    """No nonstandard JSON numbers or repr-based object serialization can enter a record."""
    with pytest.raises(core.InvalidInput):
        encode(value)


@pytest.mark.parametrize('mode', ['exact', 'numeric'])
def test_network_records_match_direct_core(mode):
    """Chain resistance is 1+1/2; all registered dispatch paths preserve core results."""
    vertices, edges = (0,1,2), ((0,1,1),(1,2,2))
    args = {'dirichlet': {'boundary': (0,2)}, 'kron_reduce': {'keep': (0,2)},
            'contract': {'partition': {'a': (0,1), 'b': (2,)}},
            'induced_subgraph': {'keep': (0,2)}, 'simplex_volume_squared': {'vertices': (0,2)}}
    from polyform_algebra.adapters.records import OPERATIONS
    net = core.network(vertices, edges, arithmetic=mode)
    for operation in OPERATIONS:
        record = execute_network(operation, vertices, edges, arithmetic=mode, arguments=args.get(operation))
        assert record['error'] is None, record
        direct = getattr(core, operation)(net, **args.get(operation, {}))
        # Context UUIDs describe runtime identities, not numerical differences.
        def strip_identity(value):
            if isinstance(value, dict):
                return {k:strip_identity(v) for k,v in value.items() if k != 'identity'}
            if isinstance(value, list):
                return list(map(strip_identity,value))
            return value
        assert strip_identity(record['value']) == strip_identity(encode(direct.value))
        assert record['status'] == ('exact' if mode == 'exact' else 'approximate')
        json.loads(canonical_json(record))


def test_digest_and_error_contract():
    """Digest includes mode, tolerance, input and semantically meaningful mapping order."""
    vertices, edges = (0,1), ((0,1,2),)
    a = execute_network('resistance',vertices,edges,arithmetic='numeric')
    b = execute_network('resistance',vertices,edges,arithmetic='numeric',policy=core.NumericPolicy())
    assert a['input_digest'] == b['input_digest']
    for changed in [execute_network('resistance',vertices,edges),
                    execute_network('resistance',vertices,((0,1,3),),arithmetic='numeric'),
                    execute_network('resistance',vertices,edges,arithmetic='numeric',policy=core.NumericPolicy(rtol=1e-8))]:
        assert changed['input_digest'] != a['input_digest']
    assert encode({'a': 1, 'b': 2}) != encode({'b': 2, 'a': 1})
    left = execute_network('contract', vertices, edges, arguments={'partition': {'a': (0,), 'b': (1,)}})
    right = execute_network('contract', vertices, edges, arguments={'partition': {'b': (1,), 'a': (0,)}})
    assert left['input_digest'] != right['input_digest']
    for record, code in [
        (execute_network('__import__("os")',vertices,edges),'UNSUPPORTED_OPERATION'),
        (execute_network('resistance',vertices,((0,1,-1),)),'INVALID_INPUT'),
        (execute_network('dirichlet',vertices,edges),'INVALID_INPUT'),
        (execute_network('resistance',vertices,edges,arguments={'unknown':1}),'INVALID_INPUT'),
        (execute_network('resistance',vertices,object()),'INVALID_INPUT')]:
        assert record['status'] == 'error' and record['value'] is None
        assert record['error']['code'] == code
        json.loads(canonical_json(record))


def test_disconnection_is_success_and_general_calculation_record():
    """An isolated pair returns typed infinity, not an error or a made-up finite value."""
    record = execute_network('resistance',(1,'1'),(),arguments={'pairs': ((1,'1'),)})
    assert record['status'] == 'exact'
    assert 'positive_infinity' in canonical_json(record['value'])
    result = core.spanning_tree_weight(core.network((0,1),((0,1,Q(1,3)),)))
    saved = calculation_record(result,input_data={'edges':((0,1,Q(1,3)),)})
    assert saved['value'] == encode(Q(1,3)) and saved['error'] is None


def test_shared_topology_preserves_legacy_components():
    """Shared traversal preserves input order, isolates and every nonzero bridge."""
    import numpy as np
    from polyform_algebra._core.topology import connected_components
    from polyform_algebra.matrix._numeric import components
    a = np.array([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
    assert components(a) == connected_components(range(4), ((0,1),(2,3))) == ((0,1),(2,3))
    a[1,2] = a[2,1] = 1
    assert components(a) == ((0,1,2,3),)
