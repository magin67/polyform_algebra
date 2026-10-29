"""Products of boundaries and their quadratic forms in the PMG introduction.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from itertools import combinations
from polyform_algebra import Point, Simplex, Polysimplex, Polyform, quForm


def main():
    """Show merging, vanishing products and the four-vertex structural identity."""
    a, b, c, d = Point.create_list("abcd")
    ab, bc = Simplex([(a,b)]), Simplex([(b,c)])
    print("Boundary product:", ab * bc)
    print("Polysimplex:", Polysimplex({ab: 2, bc: -1}))
    f_ab, f_bc, f_cd = quForm(ab), quForm(bc), quForm([(c,d)])
    f_abc = quForm([(a,b,c)])
    f_ab_cd = quForm([(a,b),(c,d)])
    print("Merged square:", f_ab*f_bc)
    print("Shared edge gives zero:", f_ab*f_abc)
    print("Disjoint boundaries:", f_ab*f_cd)
    print("Joined components:", f_ab_cd*f_bc)
    print("Cycle gives zero:", f_ab_cd*quForm([(b,c),(d,a)]))
    p, q = 2*f_ab-f_bc, f_cd.as_polyform()
    print("P:", p)
    print("Q:", q)
    print("P*Q:", p*q)
    print("3P:", 3*p)
    print("P squared:", p**2)
    # Use ordered tuples for boundaries; sets are not valid simplex components.
    triangles = Polyform.from_pairs([([t], 1) for t in combinations([a,b,c,d], 3)])
    pairs = Polyform.from_pairs([
        ([(a,b),(c,d)],1), ([(a,c),(b,d)],1), ([(a,d),(b,c)],1)
    ])
    equivalent = triangles.is_equivalent(pairs)
    assert equivalent
    print("Four-vertex identity:", equivalent)


if __name__ == "__main__":
    main()
