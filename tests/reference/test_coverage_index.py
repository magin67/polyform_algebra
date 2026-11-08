"""Keep the acceptance index navigable; this is not a mathematical proof."""
import ast
import json
from pathlib import Path


def test_all_requirements_have_existing_tests():
    """Every T01-T42 selector must resolve to a real test, preventing stale coverage claims."""
    root = Path(__file__).resolve().parents[2]
    coverage = json.loads((root / 'tests/fixtures/reference_coverage.json').read_text())
    assert set(coverage) == {'T{:02}'.format(i) for i in range(1, 43)}
    for selectors in coverage.values():
        assert selectors
        for selector in selectors:
            path, _, name = selector.partition('::')
            tree = ast.parse((root / path).read_text(encoding='utf-8'))
            tests = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)
                     and n.name.startswith('test_')}
            assert tests
            if name:
                assert name in tests, selector
