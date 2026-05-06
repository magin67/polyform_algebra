import sys, os
# sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'src'))

from objects.element import Point, Vector
from objects.simplex import Simplex
from objects.quform import quForm

a, b, c = Point.create_list(['a', 'b', 'c'])


s1 = Simplex([(a, b)])
q1 = quForm(s1)          # квадратичная форма от симплекса

print('q1: ', q1)

s2 = Simplex([(b, c)])
q2 = quForm([b,c])

# Умножение квадратичных форм
s3 = s1*s2
q3 = q1 * q2
print('q3: ', q3)

print('s3: ', s3)