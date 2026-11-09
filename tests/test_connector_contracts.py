"""Connector extension points and dependency boundaries; no external service is used."""
import ast
from pathlib import Path
import subprocess
import sys

import pytest
from polyform_algebra.connectors import DataSource, ResultSink


@pytest.mark.parametrize('contract',[DataSource,ResultSink])
def test_connector_requires_implementation(contract):
    """An abstract transport cannot pretend to read data or publish a result."""
    with pytest.raises(TypeError):
        contract()


def test_read_and_publish_are_separate_capabilities():
    """Applications can consume read-only sources without obtaining a write operation."""
    class Source(DataSource):
        def read_snapshot(self, query):
            return {'snapshot_id':'s1','objects': [{'entity_id':query['id']}]}
    source = Source()
    assert not isinstance(source,ResultSink) and not hasattr(source,'publish')
    assert source.read_snapshot({'id':'employee-1'})['objects'][0]['entity_id'] == 'employee-1'
    class Sink(ResultSink):
        def publish(self, result, *, idempotency_key):
            return {'run_id':result['run_id'],'key':idempotency_key}
    assert Sink().publish({'run_id':'r1'},idempotency_key='r1')['key'] == 'r1'
    with pytest.raises(TypeError):
        Sink().publish({'run_id':'r1'})


def test_import_boundaries_without_optional_dependencies():
    """Connectors load no mathematics; importing table/graph adapters loads no transport or optional SDK."""
    script = '''
import importlib.abc, sys
class Block(importlib.abc.MetaPathFinder):
    blocked = ('numpy','sympy','networkx','pandas','scipy','requests','httpx')
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in self.blocked:
            raise AssertionError('unexpected dependency: '+fullname)
finder=Block()
sys.meta_path.insert(0,finder)
import polyform_algebra.connectors
assert not any(x.startswith(('polyform_algebra.core','polyform_algebra.matrix','polyform_algebra.adapters')) for x in sys.modules)
finder.blocked = ('networkx','pandas','scipy','requests','httpx')
import polyform_algebra.adapters.pandas
import polyform_algebra.adapters.matrix_networkx
import polyform_algebra.adapters.networkx
'''
    result = subprocess.run([sys.executable,'-c',script],capture_output=True,text=True)
    assert result.returncode == 0, result.stderr


def test_math_does_not_depend_on_connectors_or_adapters():
    """Only the two explicitly retained matrix compatibility shims may import adapters."""
    import polyform_algebra
    root = Path(polyform_algebra.__file__).parent
    for area in ('core','matrix'):
        for path in (root/area).glob('*.py'):
            if area == 'matrix' and path.name in ('pandas.py','networkx.py'):
                continue
            tree = ast.parse(path.read_text(encoding='utf-8'))
            for node in ast.walk(tree):
                if isinstance(node,ast.ImportFrom):
                    modules = [node.module or ''] + [a.name for a in node.names]
                elif isinstance(node,ast.Import):
                    modules = [a.name for a in node.names]
                else:
                    continue
                assert not any(set(m.split('.')) & {'adapters','connectors','applications'} for m in modules), path
