"""Execute reviewed, pinned computational snippets from both public languages.

Source: https://magin67.github.io/pmg-notes/ru/pythonexamples/
v5 revision c909ddf8a590c69ed0a82fe48ca496e88fb707e0. The JSON fixture
contains verbatim Python blocks, excluding installers and subprocess wrappers.
The six module wrappers are already exercised by test_examples.py.
"""
import contextlib
import io
import json
from pathlib import Path

import numpy as np
import pytest


PAGES = json.loads((Path(__file__).parent/'fixtures/pythonexamples_v5.json').read_text(encoding='utf-8'))


@pytest.mark.parametrize('page', PAGES, ids=[p['path'] for p in PAGES])
def test_published_computational_blocks(page):
    """Public notebook-style blocks run unchanged and retain the note's numbers."""
    namespace = {'__name__': '__public_example__'}
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        for block in page['blocks']:
            exec(compile(block, page['source'], 'exec'), namespace)
    if 'tau_two' in namespace:
        assert namespace['tau_one'] == 5
        assert namespace['r12_one'] == pytest.approx(2/5)
        assert namespace['r13_one'] == pytest.approx(3/5)
        assert namespace['tau_two'] == 8
        assert namespace['r23_two'] == pytest.approx(1/2)
        assert namespace['mixed_coefficient'] == 1
    elif 'area_squared' in namespace:
        assert namespace['area_squared'] == pytest.approx(1/12)
        assert np.linalg.det(namespace['gram']) == pytest.approx(1/3)
    elif 'potential' in namespace:
        # A 2 A current drops 2 V across c=1 and 1 V across c=2.
        np.testing.assert_allclose(namespace['potential'], [3, 1, 0])
        assert namespace['potential'] @ namespace['matrix'] @ namespace['potential'] == pytest.approx(6)
    elif 'updated' in namespace:
        assert 'cycle tau = 1' in stream.getvalue()
        assert 'cycle R12 = 3' in stream.getvalue()
        assert namespace['updated'].potential(namespace['quForm'].one()) == 0
        with pytest.raises(ValueError):
            namespace['updated'].norm(namespace['f13'])
    elif 'name' in namespace:
        assert 'P5' in stream.getvalue() and 'C4' in stream.getvalue()
        assert namespace['metric'].potential(namespace['quForm'].one()) == 4
        assert namespace['metric'].norm(namespace['form']) == pytest.approx(3/4)
    else:
        assert namespace['tau'] == 3
        assert namespace['r12'] == pytest.approx(2/3)


def test_translated_computational_blocks_match():
    """Both published languages currently carry the same six computational programs."""
    en = sorted('\n'.join(p['blocks']) for p in PAGES if '/en/' in p['path'])
    ru = sorted('\n'.join(p['blocks']) for p in PAGES if '/ru/' in p['path'])
    assert len(en) == len(ru) == 6
    assert en == ru
