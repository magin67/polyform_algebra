import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import networkx as nx
import sympy as sp

# from src import AffineVector, Polyform, PolyBForm

from src.combinations.polyform import Polyform
from src.combinations.poly_bf import PolyBForm
from src.combinations.affine_vector import AffineVector

# ---------- Переменные ----------
a, b, c, d, e, f, g, h, r, x, y, z = sp.symbols('a b c d e f g h r x y z')


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
