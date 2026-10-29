"""Path norms and a scalar product on a cycle in PMG.

Source: https://github.com/magin67/pmg-notes/tree/v5/pmg-library/ru/Введение
"""
from itertools import combinations
import networkx as nx
from polyform_algebra import Point, Polyform, quForm


def main():
    """Evaluate quadratic forms after taking the network Laplacian exponential."""
    a,b,c,d = Point.create_list("abcd")
    triangles = Polyform.from_pairs([([t],1) for t in combinations([a,b,c,d],3)])
    pairs = Polyform.from_pairs([
        ([(a,b),(c,d)],1), ([(a,c),(b,d)],1), ([(a,d),(b,c)],1)
    ])
    print("Four-vertex identity:", triangles.is_equivalent(pairs))
    path = quForm([(a,b)]) + quForm([(b,c)]) + quForm([(c,d)])
    print("Path endpoint norm:", path.exp().norm(quForm([(a,d)])))
    laplacian, vertices = Polyform.from_networkx(nx.cycle_graph(8), return_points="list")
    scalar_form = (vertices[4]-vertices[0]) @ (vertices[5]-vertices[3])
    print("Cycle scalar product:", laplacian.exp().norm(scalar_form))


if __name__ == "__main__":
    main()
