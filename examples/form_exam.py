from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.simplex import Simplex
from polyform_algebra.objects.quform import quForm

from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.combinations.polyform import Polyform

a, b, c, d = Point.create_list(['a', 'b', 'c', 'd'])

# Симплексы
print("Симплексы:")

s_ab = Simplex([(a, b)])
s_bc = Simplex([(b, c)])
print(s_ab, "*", s_bc, "=", s_ab * s_bc)  # [(a, b, c)]

s_abc = Simplex([(a, b, c)])
print("s_abc: ", s_abc)  # должно вывести [(a, b, c)]

# Полисимплексы
print("Полисимплексы:")

pS = Polysimplex({s_ab: 2, s_bc: -1})
pQ = Polysimplex(s_ab)

print(pS)
print(pQ)           # quForm([[c, d]])

# Формы
print("Формы:")

f_ab = quForm([(a, b)])
f_bc = quForm([(b, c)])
f_cd = quForm([(c, d)])

f_abc = quForm(s_abc)          # квадратичная форма от симплекса

f_ab_cd = quForm([(a, b), (c, d)])
f_bc_da = quForm([(b, c), (d, a)])

print('f_abc: ', f_abc)

# Умножение форм

print(f_abc, "*", f_bc, "=", f_abc * f_bc)
print(f_ab, "*", f_bc, "=", f_ab * f_bc)  # [(a, b, c)]
print(f_ab, "*", f_abc, "=", f_ab * f_abc)  # [] (нулевая)
print(f_ab, "*", f_cd, "=", f_ab * f_cd)  # [{a, b}, {c, d}]
print(f_ab_cd, "*", f_bc, "=", f_ab_cd * f_bc)  # [{a, b, c, d}]
print(f_ab_cd, "*", f_bc_da, "=", f_ab_cd * f_bc_da)  # [] (нулевая)

# Полиформы
print("Полиформы:")

P = Polyform({f_ab: 2, f_bc: -1})
Q = Polyform(f_cd)

print("P: ", P)           # 2*quForm([[a, b]]) + -1*quForm([[b, c]])
print("Q: ", Q)           # quForm([[c, d]])
print("P*Q: ", P * Q)       # 2*quForm([[a, b, c, d]]) + -1*quForm([[b, c, d]])
print("P*3: ", P * 3)       # 6*quForm([[a, b]]) + -3*quForm([[b, c]])
print("P^2: ", P ** 2)      # (2*quForm([[a, b]]) + -1*quForm([[b, c]]))^2


# И эквивалентность. Тут проверяем тождество тетраэдра
L1 = Polyform([([{a, b, c}], 1), ([{a, b, d}], 1), ([{a, c, d}], 1), ([(b, c, d)], 1)])
L2 = Polyform([([{a, b}, {c, d}], 1), ([{a, c}, {b, d}], 1), ([{a, d}, {b, c}], 1)])
print(L1.is_equivalent(L2))