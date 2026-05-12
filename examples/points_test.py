import sys, os
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

# from _core.lincomb import LinearCombination

from src.objects.element import Point, Vector
from src.objects.simplex import Simplex
from src.combinations.polyform import Polyform

# a, b = Point('a'), Point('b')
a, b, c = Point.create_list(['a', 'b', 'c'])

print("1: ", a*b, a[0])

lincmb = (a + b) / 2
# print(lc)                     # 0.5*Point(a) + 0.5*Point(b)
# print(lc.multiplicity())      # 1.0

# Создание точки из линейной комбинации
d = Point('d', frame=lincmb)    # сохраняет lincmb
print("2: ", d[0]*a[0])         # та же комбинация

v = Vector('v', frame=b - a)  # multiplicity = 0

s1 = Simplex([a, b])
print("3: ", s1.boundary())   # 1*[b] + (-1)*[a]  -> [b] - [a]
s2 = Simplex([a, v])
print("4: ", s2.boundary())   # v имеет mult=0 -> только слагаемое для a: 1*[v] -> [v]

s3 = Simplex([(a, b), a])
print("s3: ", s3)

s4 = Simplex([(a, b), c])
print("s4: ", s4.to_polysimplex())

s5 = Simplex([(a, b), c])
print("m5: ", s2*s4)

u = Vector('u', frame=2*c - b - a)

print("v form: ", v[0]@v[0])
print("u form: ", u[0]@u[0])

print("Вектор на вектор v@u form: ", v[0]@u[0])
print("Точка на точку d@d form: ", d[0]@d[0])
print("Вектор на точку v@d form: ", v[0]@d[0])