import sys, os
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from objects.element import Point #, Vector
from _core.lincomb import LinearCombination

a, b, c = Point.create_list(['a', 'b', 'c'])

print(a*b)

lc = (a + b) / 2
print(lc)                     # 0.5*Point(a) + 0.5*Point(b)
print(lc.multiplicity())      # 1.0

# Создание точки из линейной комбинации
c = Point('c', data=lc)       # сохраняет lc как c.data
print(c.data)                 # та же комбинация