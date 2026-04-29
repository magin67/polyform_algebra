# @title Полиформа - линейная комбинация форм

from collections import Counter
import sympy as sp
import numpy as np

from .._core.lincomb import LinearCombination
# from .._core.basis import Basis

from ..objects.vector import Vector
from .affine_vector import AffineVector
from ..objects.quform import quForm
from ..objects.face import Face

class Polyform(LinearCombination):
    @classmethod
    def from_networkx(cls, G): # Создаёт полиформу 1-го порядка из (неориентированного) графа networkx
        nodes = sorted(G.nodes())
        symbols = {node: sp.symbols(chr(ord('a')+i)) for i, node in enumerate(nodes)}
        terms = Counter()
        for u, v, data in G.edges(data=True):
            weight = data.get('weight', 1)
            if weight != 0:
                a, b = symbols[u], symbols[v]
                form = quForm([{a, b}])
                terms[form] += weight
        return cls(terms)

    @classmethod
    def from_laplacian(cls, laplacian_matrix, vertex_labels, tolerance=1e-10):
        """
        Создаёт полиформу 1-го порядка из матрицы Лапласиана.
        Параметры:
            laplacian_matrix: numpy.ndarray или list списков, квадратная матрица.
            vertex_labels: список меток вершин (символы или строки), порядок соответствует строкам/столбцам.
            tolerance: порог для отбрасывания малых значений.

        Возвращает:
            Polyform: полиформа, соответствующая лапласиану.
        """
        L = np.asarray(laplacian_matrix)
        n = L.shape[0]
        if L.shape[1] != n:
            raise ValueError("Matrix must be square")
        if len(vertex_labels) != n:
            raise ValueError("Number of vertex labels must match matrix size")

        # Проверка, что матрица симметрична (с допуском)
        if not np.allclose(L, L.T, atol=tolerance):
            raise ValueError("Matrix must be symmetric")

        terms = Counter()
        for i in range(n):
            for j in range(i+1, n):
                weight = -L[i, j]  # так как L_ij = -w_ij для i≠j
                if abs(weight) > tolerance: # Приводим к каноническому виду (сортировка вершин)
                    a, b = sorted([vertex_labels[i], vertex_labels[j]], key=str)
                    form = quForm([{a, b}])
                    if not form.is_zero(): terms[form] += weight
        return cls(terms)

    @classmethod
    def from_components(cls, components):
        """
        Создаёт полиформу как произведение компонент.
        Каждая компонента может быть AffineVector, Face или quForm.
        """
        result = cls(quForm.one())
        for comp in components:
            if isinstance(comp, AffineVector):
                poly = comp.to_polyform()
            elif isinstance(comp, Face):
                poly = cls(comp.to_quform())
            elif isinstance(comp, quForm):
                poly = cls(comp)
            else:
                raise TypeError(f"Unsupported component type: {type(comp)}")
            result = result * poly
        return result

    def __init__(self, data=None, term_type=quForm):
        if data is None:
            super().__init__(term_type=term_type)
            return
        prepared = {}
        if isinstance(data, dict):
            for key, coeff in data.items():
                if coeff == 0: continue
                if not isinstance(key, term_type): key = term_type(key)
                prepared[key] = prepared.get(key, 0) + coeff
        elif isinstance(data, (list, tuple)):
            for item in data:
                if len(item) != 2:
                    raise ValueError("Each item must be (form_data, coeff)")
                key, coeff = item
                if coeff == 0: continue
                if not isinstance(key, term_type): key = term_type(key)
                prepared[key] = prepared.get(key, 0) + coeff
        elif isinstance(data, term_type):
            prepared = {data: 1}
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        super().__init__(prepared, term_type=term_type)

    def extract_max_coeff(self, grade: int = None, comp_count: int = None, by_abs: bool = True) -> 'Polyform':
        """
        Возвращает полиформу из слагаемых, у которых коэффициент максимален,
        с возможной фильтрацией по грейду и/или количеству компонент.

        grade: если задан, учитываются только формы этого грейда.
        comp_count: если задан, учитываются только формы с таким числом компонент.
        by_abs: если True, сравниваются абсолютные значения коэффициентов.
        """
        filtered = []
        for form, coeff in self.terms.items():
            if grade is not None and form.grade() != grade: continue
            if comp_count is not None and form.size() != comp_count: continue
            filtered.append((form, coeff))
        if not filtered: return Polyform()
        if by_abs:
            max_val = max(abs(c) for _, c in filtered)
            selected = [(form, coeff) for form, coeff in filtered if abs(coeff) == max_val]
        else:
            max_val = max(c for _, c in filtered)
            selected = [(form, coeff) for form, coeff in filtered if coeff == max_val]
        return Polyform(dict(selected))

    def extract_max_coeff_all_grades(self, comp_count: int = None, by_abs: bool = True) -> 'Polyform':
        """
        Для каждого грейда оставляет слагаемые с максимальным коэффициентом,
        с возможной фильтрацией по количеству компонент. Объединяет в одну полиформу.
        """
        max_grade = self.max_grade()
        result_terms = {}
        for g in range(max_grade + 1):
            sub = self.extract_max_coeff(grade=g, comp_count=comp_count, by_abs=by_abs)
            for form, coeff in sub.terms.items():
                result_terms[form] = result_terms.get(form, 0) + coeff
        return Polyform(result_terms)

    def extract_diagonal_components(self) -> 'Polyform':
        """
        Возвращает полиформу, содержащую для каждого грейда k только те формы,
        у которых количество компонент = max_grade - k + 1.
        """
        mg = self.max_grade()
        if mg < 0: return Polyform()
        result_terms = {}
        for form, coeff in self.terms.items():
            k = form.grade()
            comp_cnt = form.size()
            if comp_cnt == mg - k + 1:
                result_terms[form] = result_terms.get(form, 0) + coeff
        return Polyform(result_terms)

    # Метрика
    def exp(self):
        """
        Быстрое вычисление экспоненты для полиформы, диагональной в ортогональном базисе.
        Предполагается, что self = Σ coeff_i * b_i, где b_i — bForm из одного вектора.
        Возвращает Product(1 + coeff_i * b_i).
        """
        if not self.is_homogeneous() or self.grade() != 1: return super().exp() # Если не диагональная, используем стандартный (медленный) метод из родителя
        self._skip_basis_update = True
        result = Polyform({quForm.one(): 1})   # единица
        for frm, coeff in self.terms.items():
            # Создаём терм (1 + coeff * bf) как полиформу
            term_terms = {quForm.one(): 1}
            if coeff != 0:
                term_terms[frm] = term_terms.get(frm, 0) + coeff
            term = Polyform(term_terms)
            result = result * term
        self._skip_basis_update = False
        return result

    # Потенциал и норма
    def potential(self, obj_poly: 'Polyform'):
        if self.is_zero(): return 0
        g_metric = self.grade()
        g_obj = obj_poly.grade()
        if self.is_homogeneous(): # если полиформы однорордная, то считаем, что уже выбран нужный слой
          metric_part = self
        else: # иначе считаем, что задана общая полиформа, из которой надо предварительно извлечь нужный слой
          if g_obj < 0 or g_metric < g_obj: return 0
          metric_part = self.extract_grade(g_metric - g_obj)
        if metric_part.is_zero(): return 0
        prod = metric_part * obj_poly
        if prod.is_zero(): return 0
        lead = prod.leading_coefficient()
        if isinstance(lead, list): return sum(lead)
        return lead

    def norm(self, obj_poly: 'Polyform'):
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
        # Проверяет эквивалентность двух полиформ с точностью до нулевых комбинаций
        if not isinstance(other, Polyform): return False
        diff = self - other
        # Если diff уже нулевая, то эквивалентны
        if diff.is_zero(): return True
        # Получаем все вершины из разности
        vertices = sorted(diff.get_vertices(), key=str)
        # Генерируем все базисные векторы a_i - a_j
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                v = AffineVector({vertices[i]: 1, vertices[j]: -1})
                prod = diff * v.to_polyform()
                if not prod.is_zero(): return False
        return True

    # Спектры
    def to_laplacian(self):
        """Преобразует однородную полиформу 1-го порядка в матрицу лапласиана.
        Возвращает (matrix, vertex_index), где matrix — np.ndarray,
        vertex_index — dict {символ: индекс}.
        """
        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")

        # Собираем все вершины из всех форм
        vertices = sorted(self.get_elements(), key=str)
        n = len(vertices)
        idx = {v: i for i, v in enumerate(vertices)}

        L = np.zeros((n, n))
        for form, coeff in self.terms.items():
            # Форма 1-го порядка — это один кортеж из двух вершин
            if len(form.components) != 1:
                raise ValueError("Форма должна состоять из одной компоненты")
            a, b = form.components[0]  # компонента — кортеж (вершина1, вершина2)
            if a not in idx or b not in idx: continue  # на всякий случай
            i, j = idx[a], idx[b]
            L[i, j] -= coeff
            L[j, i] -= coeff
            L[i, i] += coeff
            L[j, j] += coeff
        return L, idx

    def eigenvectors(self, normalized=False, tolerance=1e-10):
        """
        Возвращает список (собственное_значение, аффинный_вектор) для полиформы 1-го порядка.
        Вектор масштабирован так, что его норма (сумма квадратов) равна собственному значению.
        Параметр tolerance задаёт точность округления (значения меньше tolerance считаются нулём,
        а близкие к целому округляются до целого).
        """
        if not self.is_homogeneous() or self.max_grade() != 1:
            raise ValueError("Полиформа должна быть однородной 1-го порядка")

        L, idx = self.to_laplacian()
        # Приводим idx к формату {символ: индекс} (если нужно)
        if isinstance(next(iter(idx.keys())), int):
            idx = {v: k for k, v in idx.items()}
        eig_vals, eig_vecs = np.linalg.eigh(L)
        eig_vals, eig_vecs = np.real(eig_vals), np.real(eig_vecs)

        result = []
        for i in range(len(eig_vals)):
            ev = eig_vals[i]
            if abs(ev) < tolerance: continue # пропускаем векторы с нулевым собственным значением
            vec = eig_vecs[:, i]
            if normalized: vec_scaled = vec
            else: vec_scaled = vec * np.sqrt(ev)
            coeffs = {}
            for vertex, j in idx.items():
                val = self._round_value(vec_scaled[j], tolerance)
                if val != 0.0: coeffs[vertex] = val
            avec = Vector(coeffs)
            result.append((ev, avec))
        return result