"""Exact series resistance as a traceable record, PMG Introduction weighted chain.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/en/Getting%20Started%20and%20Foundations/Introduction
Definition revision is pinned by the record adapter; this example needs no AI SDK.
"""
from polyform_algebra.adapters.records import execute_network, canonical_json


def main():
    result = execute_network('resistance', ('a','b','c'),
                             (('a','b',1),('b','c',2)),
                             arguments={'pairs': (('a','c'),)})
    assert result['status'] == 'exact'
    # Two resistors in series give 1/1 + 1/2 = 3/2.
    assert result['value']['fields']['values']['items'][0] == {
        'type':'rational','numerator':'3','denominator':'2'}
    print(canonical_json(result))


if __name__ == '__main__':
    main()
