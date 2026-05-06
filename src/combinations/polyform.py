# @title Полиформа - линейная комбинация форм

from collections import Counter
import numpy as np

from .polysimplex import Polysimplex
from objects.quform import quForm
from objects.element import Point, Vector   # новый класс вектора

class PolyForm(Polysimplex):
    @classmethod
    def one(cls):
        return cls({quForm.one(): 1})

    @classmethod
    def from_networkx(cls, G):
        nodes = sorted(G.nodes())
        # Создаём точки с именами, равными строковому представлению узлов
        points = {node: Point(str(node)) for node in nodes}
        terms = Counter()
        for u, v, data in G.edges(data=True):
            weight = data.get('weight', 1)
            if weight != 0:
                a = points[u]
                b = points[v]
                # Используем кортеж (a,b) — граница ориентированного ребра
                form = quForm([(a, b)])
                terms[form] += weight
        return cls(terms)

    @classmethod
    def from_laplacian(cls, laplacian_matrix, vertex_labels, tolerance=1e-10):
        L = np.asarray(laplacian_matrix)
        n = L.shape[0]
        if L.shape[1] != n:
            raise ValueError("Matrix must be square")
        if len(vertex_labels) != n:
            raise ValueError("Number of vertex labels must match matrix size")
        if not np.allclose(L, L.T, atol=tolerance):
            raise ValueError("Matrix must be symmetric")
        terms = Counter()
        for i in range(n):
            for j in range(i+1, n):
                weight = -L[i, j]
                if abs(weight) > tolerance:
                    a, b = sorted([vertex_labels[i], vertex_labels[j]], key=str)
                    form = quForm([(a, b)])
                    if not form.is_zero():
                        terms[form] += weight
        return cls(terms)

    @classmethod
    def from_components(cls, components):
        result = cls.one()
        for comp in components:
            if isinstance(comp, Vector):
                # Преобразуем вектор в quForm: для вектора v = b - a -> quForm([(a,b)])?
                # Старое: AffineVector.to_polyform() возвращал Polyform с термом quForm.
                # Реализуем через представление вектора:
                # Если у вектора есть frame (линейная комбинация) 
                # и он имеет multiplicity=0, то frame — это b - a.
                # Но для упрощения пока считаем, что вектор уже есть как элемент.
                # Временно: vector -> quForm([tuple(вектор?)]) - неверно.
                # Пропустим, пока не определимся.
                raise NotImplementedError("Conversion from Vector to PolyForm not implemented")
            elif isinstance(comp, quForm):
                result = result * cls({comp: 1})
            else:
                raise TypeError(f"Unsupported component type: {type(comp)}")
        return result

    def __init__(self, data=None):
        if data is None:
            super().__init__(term_type=quForm)
            return
        prepared = {}
        if isinstance(data, dict):
            items = data.items()
        elif isinstance(data, (list, tuple)):
            items = data
        elif isinstance(data, quForm):
            items = [(data, 1)]
        elif isinstance(data, Polysimplex):
            # Если передан Polysimplex, преобразуем его термы (Simplex) в quForm
            for term, coeff in data.terms.items():
                qterm = quForm(term)
                self.add_term(prepared, qterm, coeff)
            super().__init__(prepared, term_type=quForm)
            return
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        for term, coeff in items:
            if coeff == 0:
                continue
            if not isinstance(term, quForm):
                term = quForm(term)
            self.add_term(prepared, term, coeff)
        super().__init__(prepared, term_type=quForm)

    @staticmethod
    def add_term(prepared, term, coeff):
        cterm = term.canonical()
        key = quForm(cterm.components, sign=1)
        prepared[key] = prepared.get(key, 0) + coeff * cterm.sign

    def __str__(self):
        if self.is_zero():
            return "0"
        parts = []
        for term, coeff in self.terms.items():
            term_str = str(term)   # quForm выводит f[...]
            if coeff == 1:
                parts.append(term_str)
            elif coeff == -1:
                parts.append(f"-{term_str}")
            else:
                parts.append(f"{coeff}*{term_str}")
        result = parts[0]
        for part in parts[1:]:
            if part[0] == '-':
                result += f" - {part[1:]}"
            else:
                result += f" + {part}"
        return result

    def __repr__(self):
        return f"PolyForm({self.__str__()})"

    # --- Методы, специфичные для квадратичных форм ---
    def potential(self, obj_poly):
        if self.is_zero():
            return 0
        g_metric = self.grade()
        g_obj = obj_poly.grade()
        if self.is_homogeneous():
            metric_part = self
        else:
            if g_obj < 0 or g_metric < g_obj:
                return 0
            metric_part = self.extract_grade(g_metric - g_obj)
        if metric_part.is_zero():
            return 0
        prod = metric_part * obj_poly
        if prod.is_zero():
            return 0
        lead = prod.leading_coefficient()
        if isinstance(lead, list):
            return sum(lead)
        return lead

    def norm(self, obj_poly):
        if self.is_zero():
            raise ValueError("Cannot compute norm with zero metric")
        pot = self.potential(obj_poly)
        metric_lead = self.leading_coefficient()
        if isinstance(metric_lead, list):
            metric_lead = sum(metric_lead)
        if metric_lead == 0:
            raise ValueError("Metric leading coefficient is zero")
        return pot / metric_lead

    def is_equivalent(self, other, tolerance=1e-10):
        if not isinstance(other, PolyForm):
            return False
        diff = self - other
        if diff.is_zero(tolerance):
            return True
        vertices = sorted(diff.get_elements(), key=str)
        # Генерируем базисные векторы a_i - a_j как новые векторы
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                # Создаём вектор v = b - a (как линейную комбинацию)
                # Представление вектора: b - a — это quForm([(a,b)])? Нет, вектор — это элемент, а не форма.
                # Для проверки эквивалентности нужно умножить diff на вектор (как quForm?).
                # В старом коде использовался AffineVector.to_polyform() -> Polyform с термом quForm.
                # Создадим quForm из разности двух базисных точек.
                # Но проще: создать симплекс [(a,b)] и взять его границу? Нет.
                # Временно: предполагаем, что вектор задаётся как quForm([(a,b)]), что соответствует (b-a).
                q = quForm([(vertices[i], vertices[j])])  # это b - a? Да, quForm([(a,b)]) соответствует b-a.
                prod = diff * PolyForm({q: 1})
                if not prod.is_zero(tolerance):
                    return False
        return True

    def to_laplacian(self):
        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")
        vertices = sorted(self.get_elements(), key=str)
        n = len(vertices)
        idx = {v: i for i, v in enumerate(vertices)}
        L = np.zeros((n, n))
        for term, coeff in self.terms.items():
            # term — quForm (симплекс), его компоненты — кортежи или атомы
            if len(term.components) != 1:
                raise ValueError("Форма должна состоять из одной компоненты")
            comp = term.components[0]
            if isinstance(comp, tuple) and len(comp) == 2:
                a, b = comp
                if a not in idx or b not in idx:
                    continue
                i, j = idx[a], idx[b]
                L[i, j] -= coeff
                L[j, i] -= coeff
                L[i, i] += coeff
                L[j, j] += coeff
            else:
                # возможен атомарный случай: term = quForm([a,b]) — это два атома?
                # Но для формы 1-го порядка должны быть только кортежи из 2 элементов
                raise ValueError("Неподдерживаемый формат формы для лапласиана")
        return L, idx

    # def eigenvectors(self, normalized=False, tolerance=1e-10):
    #     if not self.is_homogeneous() or self.grade() != 1:
    #         raise ValueError("Полиформа должна быть однородной 1-го порядка")
    #     L, idx = self.to_laplacian()
    #     eig_vals, eig_vecs = np.linalg.eigh(L)
    #     eig_vals, eig_vecs = np.real(eig_vals), np.real(eig_vecs)
    #     result = []
    #     for i in range(len(eig_vals)):
    #         ev = eig_vals[i]
    #         if abs(ev) < tolerance:
    #             continue
    #         vec = eig_vecs[:, i]
    #         if not normalized:
    #             vec = vec * np.sqrt(ev)
    #         # Собираем аффинный вектор как линейную комбинацию точек (коэффициенты при вершинах)
    #         coeffs = {}
    #         for vertex, j in idx.items():
    #             val = self._round_value(vec[j], tolerance)
    #             if val != 0.0:
    #                 coeffs[vertex] = val
    #         # Создаём новый Vector (с multiplicity=0)
    #         # Вектор задаётся как линейная комбинация точек, например, sum coeff * Point(vertex)
    #         # Но нам нужен объект Vector, который можно хранить в представлении.
    #         # Можно создать Vector без явного имени, передав ему PolySimplex?
    #         # Вектор — это элемент с multiplicity=0, и его frame может быть суммой точек с коэффициентами.
    #         # Создадим просто экземпляр Vector, передав ему rep = PolySimplex({...}).
    #         # Пока упростим: создадим Vector с именем, но без frame.
    #         avec = Vector(name=f"eigenvec_{i}", multiplicity=0)
    #         # Установим frame для вектора как линейную комбинацию.
    #         # frame — это Polysimplex, где термы — singletons (точки).
    #         polysimp = Polysimplex()
    #         for vertex, val in coeffs.items():
    #             # создаём симплекс из одной точки
    #             p = Simplex([vertex], sign=1)
    #             polysimp += Polysimplex({p: val})
    #         avec.add_representation(polysimp)
    #         result.append((ev, avec))
    #     return result