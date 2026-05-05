import sys, os
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

# from _core.lincomb import LinearCombination

from objects.element import Point, Vector
from objects.simplex import Simplex

# a, b = Point('a'), Point('b')
a, b, c = Point.create_list(['a', 'b', 'c'])

print(a*b, a[0])

lc = (a + b) / 2
# print(lc)                     # 0.5*Point(a) + 0.5*Point(b)
# print(lc.multiplicity())      # 1.0

# Создание точки из линейной комбинации
d = Point('d', frame=lc)       # сохраняет lc
print(d[0]*a[0])                 # та же комбинация

v = Vector('v', frame=b - a)  # multiplicity = 0

s1 = Simplex([a, b])
print(s1.boundary())   # 1*[b] + (-1)*[a]  -> [b] - [a]
s2 = Simplex([a, v])
print(s2.boundary())   # v имеет mult=0 -> только слагаемое для a: 1*[v] -> [v]