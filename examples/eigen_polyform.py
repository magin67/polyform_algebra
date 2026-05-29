import networkx as nx

from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.quform import quForm

from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.combinations.polyform import Polyform
from polyform_algebra._core.basis import Basis

# ---------- Переменные ----------


# 1. Создаем полиформу и ее собственное представление

# Cycle4 = Polyform([([{a, b}], 1), ([(b, c)], 1), ([(c, d)], 1), ([(d, a)], 1)]) # 4-цикл
# Cycle5 = Polyform([([{a, b}], 1), ([(b, c)], 1), ([(c, d)], 1), ([(d, e)], 1), ([(e, a)], 1)]) # 5-цикл
G_cycle = nx.cycle_graph(5)

poly = Polyform.from_networkx(G_cycle)
eig_poly = PolyBForm.from_eigenvectors(poly, True) # И преобразуем ее в полиформу на собственных векторах
print(eig_poly)
print(eig_poly*eig_poly/2)
# print(eig_poly.basis, eig_poly.basis.size)

# 2. Берем две экспоненты - от исходной полиформы и от собственной

W = poly.exp()
eW = eig_poly.exp()
print(W.sizes(), eW.sizes()) # размеры разные

# 3. Создаем аффинный вектор
av = AffineVector({a: -1, b: 2, d: -1})

quv = av.to_polyform() # полиформа вектора в исходном базисе
euv = av.to_polyform_in_basis(eig_poly.basis) # полиформа вектора в собственном базисе

print("Потенциалы: ", W.potential(quv), eW.potential(euv)) # потенциалы разные - зависят от базиса
print("Нормы: ", W.norm(quv), eW.norm(euv)) # а нормы одинаковые - не зависят от базиса
