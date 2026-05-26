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

print("L1*ab: ", L1*f_ab)
print("L2*ab: ", L2*f_ab)

print("[(a, b, c)]*[(a, b)]: ", quForm([(a, b, c)])*quForm([(a, b)]))

print("L1 == L2: ", L1.is_equivalent(L2))

# Cycle3 = Polyform({f_ab: 1, f_bc: 1, f_ac: 1})
# Path4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1})
# Cycle4 = Polyform({f_ab: 1, f_bc: 1, f_cd: 1, f_ad: 1})

# print("Path4.exp: ", Path4.exp())

# print("P.laplacian: ", P.to_laplacian())

# После последних изменений отвалилось умножение полиформ на форму. Самого кода умножения в классе полиформ нет, скорее всего используется код из линейных комбинаций. Он там такой:
#     def __mul__(self, other):
#         if isinstance(other, (int, float)): return self.__rmul__(other)
#         if isinstance(other, LinearCombination):
#             if self.is_zero() or other.is_zero(): return self.zero()
#             if self.is_one(): return other
#             if other.is_one(): return self
#             if self.term_type != other.term_type:
#                 raise TypeError("Cannot multiply combinations with different term types")
#             return self._mul_terms(other)
#         if isinstance(other, Monoid): # Проверка совместимости типа
#             if self.term_type is not None and not isinstance(other, self.term_type):
#                 raise TypeError(f"Cannot multiply {self.__class__.__name__} of {self.term_type} by {type(other)}")
#             return self * self.__class__({other: 1})
#         raise TypeError(f"Cannot multiply {self.__class__.__name__} by {type(other)}")

#     def __rmul__(self, scalar):
#         if isinstance(scalar, (int, float)):
#             new = self.__class__()
#             new.terms = Counter({k: v * scalar for k, v in self.terms.items()})
#             new.term_type = self.term_type   # добавлено
#             new._finalize()
#             new._update_cached_grades()
#             return new
#         raise TypeError(f"Cannot multiply {self.__class__.__name__} by {type(scalar)}")