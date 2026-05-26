from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.quform import quForm

from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.combinations.polyform import Polyform
from polyform_algebra._core.basis import Basis

import numpy as np
import networkx as nx

a, b, c, d = Point.create_list(['a', 'b', 'c', 'd'])

# Формы
f_ab = quForm([(a, b)])
f_bc = quForm([(b, c)])
f_ac = quForm([(a, c)])
f_cd = quForm([(c, d)])
f_ad = quForm([(a, d)])

f_abc = quForm([(a, b, c)])

f_ab_cd = quForm([(a, b), (c, d)])
f_bc_da = quForm([(b, c), (d, a)])

# Полиформы

Cycle3 = Polyform({f_ab: 1, f_bc: 1, f_ac: 1})
Path4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1})
Cycle4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1, f_ad: 1})

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
print("[v]^2 in eig: ", f_veig)

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
