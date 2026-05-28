from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.quform import quForm

from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.combinations.polyform import Polyform
from polyform_algebra._core.basis import Basis

import numpy as np
import networkx as nx

a, b, c, d, e, f, g, h = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'])

# Формы
f_ab = quForm([(a, b)])
f_bc = quForm([(b, c)])
f_cd = quForm([(c, d)])
f_de = quForm([(d, e)])

f_ac = quForm([(a, c)])
f_ad = quForm([(a, d)])
f_ae = quForm([(a, e)])

f_abc = quForm([(a, b, c)])
f_abd = quForm([(a, b, d)])
f_acd = quForm([(a, c, d)])
f_bcd = quForm([(b, c, d)])

f_ab_cd = quForm([(a, b), (c, d)])
f_bc_da = quForm([(b, c), (d, a)])
f_ad_bc = quForm([(a, d), (b, c)])

f_abcd = quForm([(a, b, c, d)])

# Полиформы
print("Полиформы:")

# P1 = Polyform({f_ab: 2, f_bc: -1})
P = Polyform(2*f_ab - f_bc)
Q = Polyform(f_cd)

# print("P: ", P)           # 2*quForm([[a, b]]) + -1*quForm([[b, c]])
# print("Q: ", Q)           # quForm([[c, d]])
# print("P*Q: ", P * Q)       # 2*quForm([[a, b, c, d]]) + -1*quForm([[b, c, d]])
# print("P*3: ", P * 3)       # 6*quForm([[a, b]]) + -3*quForm([[b, c]])
# print("P^2: ", P ** 2)      # (2*quForm([[a, b]]) + -1*quForm([[b, c]]))^2

# print("P*f_bc: ", P*f_bc)

# print("qu1*qu1: ", quForm.one()*quForm.one())
# print("Poly1*Poly1: ", Polyform.one() * Polyform.one())

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

# print("L1: ", L1)
# print("L2: ", L2)
# print("L1*ab: ", L1*f_ab)
# print("L2*ab: ", L2*f_ab)
# print("[(a, b, c)]*[(a, b)]: ", quForm([(a, b, c)])*quForm([(a, b)]))

# print("L1 == L2: ", L1.is_equivalent(L2))

Path3 = f_ab + f_bc
Path4 = Path3 + f_cd
Path5 = Path4 + f_de
Cycle3 = Path3 + f_ac
Cycle4 = Path4 + f_ad
Cycle5 = Path5 + f_ae

plMetricPath4 = Path4.exp()
plMetric = plMetricPath4
print("Metric: ", plMetric)
print("norm(ac): ", plMetric.norm(f_ac))

# --- Скалярное произведение векторов ---
plMetricCycle = Cycle5.exp()
plMetric = plMetricCycle #plMetricPath4

sprod = (b - a)@(c - e)
print("potential{(v1)@(v2)}: ", plMetric.potential(sprod))

# --- Скалярное произведение треугольников ---

f_acd_bcd = f_acd + f_bcd - f_ab_cd
f_abc_bcd = f_abc + f_bcd - f_ad_bc

plMetric = Path5.exp() #plMetricCycle #plMetricPath4
print("norm{acd, bcd}: ", plMetric.potential(f_acd_bcd))
print("norm{abc, bcd}: ", plMetric.potential(f_abc_bcd))