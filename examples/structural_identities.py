"""Check the tetrahedron identity and the rank on four vertices.

PMG v5, Basis and Dimension of the Space of Partition Forms:
https://github.com/magin67/pmg-notes/blob/c1c017893e96e221061b1e33bfc034891bcd0518/pmg-library/en/Polyform%20Algebra%20and%20Combinatorics/Structural%20identities/Basis%20and%20Dimension%20of%20the%20Space%20of%20Partition%20Forms.md
"""

from sympy import Matrix

from polyform_algebra.core import Context, equivalent, form_coefficients, partition_form


def main():
    context = Context(('a', 'b', 'c', 'd'))
    partitions = (
        (('a', 'b'), ('c', 'd')), (('a', 'c'), ('b', 'd')), (('a', 'd'), ('b', 'c')),
        (('a', 'b', 'c'),), (('a', 'b', 'd'),), (('a', 'c', 'd'),), (('b', 'c', 'd'),),
    )
    forms = [partition_form(context, blocks).value for blocks in partitions]
    left = sum(forms[1:3], forms[0])
    right = sum(forms[4:], forms[3])
    assert equivalent(left, right).value
    assert not equivalent(left, right + forms[3]).value

    columns = [dict(form_coefficients(form).value) for form in forms]
    coordinates = sorted(set().union(*columns))
    matrix = Matrix([[column.get(pair, 0) for column in columns] for pair in coordinates])
    print('Tetrahedron identity:', equivalent(left, right).value)
    print('Rank:', matrix.rank(), '| Independent relations:', len(matrix.nullspace()))


if __name__ == '__main__':
    main()
