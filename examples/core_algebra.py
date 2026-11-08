"""Exact formal algebra, finite graph exponential and normalized potentials.

Source: PMG v5, weighted triangle and metric identity:
https://github.com/magin67/pmg-notes/blob/c909ddf8a590c69ed0a82fe48ca496e88fb707e0/pmg-library/en/Getting%20Started%20and%20Foundations/Introduction/Laplacian%20Exponential.md
"""
from fractions import Fraction

from polyform_algebra.core import (
    Context, Atom, Boundary, Unit, exterior, boundary, wedge, square,
    form_product, form_exp, laplacian_form, spanning_coefficient,
    potential, evaluate_form, equivalent, network,
)


def main():
    net = network(('a', 'b', 'c'), (('a', 'b', 2), ('b', 'c', 3), ('a', 'c', 4)))
    c = net.context
    l = laplacian_form(net).value
    m = form_exp(l).value
    f = square(Boundary(c, (0, 1))).value
    tau = spanning_coefficient(m).value
    print('Spanning coefficient:', tau)  # 2*3 + 2*4 + 3*4 = 26.
    print('Potential of [(ab)]^2:', potential(m, f).value)
    print('Normalized potential:', potential(m, f).value/tau)
    print('Direct metric evaluation:', evaluate_form(net, f).value)
    inverse = form_exp(-l).value
    product = form_product(m, inverse).value
    print('exp(L) exp(-L) = e:', equivalent(product, square(Unit(c)).value).value)
    print('Small exact coefficient:', (Fraction(1, 10**12)*f).terms[0][1])

    weighted = Context(('a', 'b', 'c'), epsilon=(1, 2, 1))
    ab = boundary(exterior(weighted, {(0, 1): 1})).value
    bc = boundary(exterior(weighted, {(1, 2): 1})).value
    abc = boundary(exterior(weighted, {(0, 1, 2): 1})).value
    print('Weighted boundary identity:', wedge(ab, bc).value == 2*abc)
    print('Orientation:', wedge(Atom(c, 0), Atom(c, 1)).value.coefficients)


if __name__ == '__main__':
    main()
