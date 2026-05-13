import sys, os
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from objects.element import Point, Vector
from objects.quform import quForm

from combinations.polysimplex import Polysimplex
from combinations.polyform import Polyform
from _core.basis import Basis

import numpy as np
import networkx as nx

a, b, c, d = Point.create_list(['a', 'b', 'c', 'd'])

# Формы
# print("Формы:")

f_ab = quForm([(a, b)])
f_bc = quForm([(b, c)])
f_ac = quForm([(a, c)])
f_cd = quForm([(c, d)])
f_ad = quForm([(a, d)])

f_abc = quForm([(a, b, c)])

f_ab_cd = quForm([(a, b), (c, d)])
f_bc_da = quForm([(b, c), (d, a)])

# Полиформы
print("Полиформы:")

# P1 = Polyform({f_ab: 2, f_bc: -1})
P = Polyform(2*f_ab - f_bc)
Q = Polyform(f_cd)

print("P: ", P)           # 2*quForm([[a, b]]) + -1*quForm([[b, c]])
print("Q: ", Q)           # quForm([[c, d]])
print("P*Q: ", P * Q)       # 2*quForm([[a, b, c, d]]) + -1*quForm([[b, c, d]])
print("P*3: ", P * 3)       # 6*quForm([[a, b]]) + -3*quForm([[b, c]])
print("P^2: ", P ** 2)      # (2*quForm([[a, b]]) + -1*quForm([[b, c]]))^2

print("P*f_bc: ", P*f_bc)

# Эквивалентность. Проверяем тождество тетраэдра

L1 = Polyform.from_pairs([
    ([(a, b, c)], 1),
    ([(a, b, d)], 1),
    ([(a, c, d)], 1),
    ([(b, c, d)], 1)
])
L2 = Polyform.from_pairs([
    ([(a, b), (c, d)], 1),
    ([(a, c), (b, d)], 1),
    ([(a, d), (b, c)], 1)
])

print("L1: ", L1)
print("L2: ", L2)
print("L1 == L2: ", L1.is_equivalent(L2))

# print("P.laplacian: ", P.to_laplacian())

Cycle3 = Polyform({f_ab: 1, f_bc: 1, f_ac: 1})
Path4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1})
Cycle4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1, f_ad: 1})

G_petersen = nx.petersen_graph()
L_petersen = Polyform.from_networkx(G_petersen)

PathEig, B_eigen = Path4.to_eigenbasis()
print("Path.eigen: ", PathEig)
print("Path.eigen.expand: ", PathEig.expand())

points = list(Path4.get_elements())

print("Basis of v: ", B_eigen)

B_points = B_eigen.create_basis_from_frames(0)

B_eigen.link_to(B_points)

print("Basis of p[0]: ", B_points[0][1])

v = Vector('v', frame = (d - a))
f_v = v[0]@v[0]

v.add_frame(frame=v.in_basis(0, B_points, B_eigen), basis=B_eigen)

f_veig = v[1].to_quadratic_form(B_eigen)

print("v = (d - a) in eig: ", v[1])
print("<v>^2 in eig: ", f_veig)

# Метрика:

W = Path4# Cycle3
Weig = PathEig # W.to_eigenbasis()

Pmetric = W.exp()
Peigmetric = Weig.exp()

# print("W sizes: ", Pmetric.sizes())
# print("Eig sizes: ", Peigmetric.sizes())

print("f_v.norm: ", Pmetric.norm(f_v))
print("f_veig.norm: ", Peigmetric.norm(f_veig))
# print("PathEig.Exp: ", PathEig.exp())
