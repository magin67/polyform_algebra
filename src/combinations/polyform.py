# @title Полиформа - линейная комбинация форм

from collections import Counter
import numpy as np

from combinations.polysimplex import Polysimplex
from objects.quform import quForm
from objects.element import Element, Point   # новый класс вектора

class Polyform(Polysimplex):
    @classmethod
    def zero(cls):
        return cls()
    
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
    def from_pairs(cls, pairs): # pairs: список пар (key, coeff), где key может быть списком или кортежем
        new_data = {}
        for key, value in pairs:
            new_data[quForm(key)] = value
        return cls(new_data)

    def __init__(self, data=None):
        super().__init__(data, term_type=quForm)

    @staticmethod
    def add_term(prepared, term, coeff):
        cterm = term.canonical()
        key = quForm(cterm.components, sign=1)
        prepared[key] = prepared.get(key, 0) + coeff * cterm.sign

    def __str__(self):
        if self.is_zero(): return "<>"
        if self.is_one(): return "<1>"
        parts = []
        for term, coeff in self.terms.items():
            term_str = str(term)  # quForm выводит с префиксом f
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
        return f"Polyform({self.__str__()})"

    def is_equivalent(self, other, tolerance=1e-10):
        # Проверяет эквивалентность двух полиформ с точностью до нулевых комбинаций
        if not isinstance(other, Polyform): return False
        diff = self - other
        # Если diff уже нулевая, то эквивалентны
        if diff.is_zero(): return True
        # Получаем все вершины из разности
        vertices = sorted(diff.get_elements(kind='points'), key=str)
        # Генерируем все базисные векторы a_i - a_j
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                plf = quForm([(vertices[i], vertices[j])])
                prod = diff * plf
                if not prod.is_zero(): return False
        return True

    # --- Методы метрики ---
    def potential(self, obj_polyform):
        # obj_polyform - quForm or Polyform
        if self.is_zero(): return 0
        g_metric = self.grade()
        g_obj = obj_polyform.grade()
        if self.is_homogeneous():
            metric_part = self
        else:
            if g_obj < 0 or g_metric < g_obj: return 0
            metric_part = self.extract_grade(g_metric - g_obj)
        if metric_part.is_zero(): return 0
        prod = metric_part * obj_polyform
        if prod.is_zero(): return 0
        lead = prod.leading_coefficient()
        if isinstance(lead, list): return sum(lead)
        return lead

    def norm(self, obj_poly):
        pot = self.potential(obj_poly)
        if pot == 0: return 0
        metric_lead = self.leading_coefficient()
        if isinstance(metric_lead, list):
            metric_lead = sum(metric_lead)
        if metric_lead == 0:
            raise ValueError("Metric leading coefficient is zero")
        return pot / metric_lead

    @staticmethod
    def scalar_product(frame1, frame2, tolerance=1e-10):
        # Проверка, что фреймы состоят только из точек
        for frame in (frame1, frame2):
            if frame.get_elements('all') != frame.get_elements('points'):
                raise ValueError("Фрейм содержит векторы или границы")

        mult1 = frame1.multiplicity()
        mult2 = frame2.multiplicity()

        points = sorted(frame1.get_elements('points') | frame2.get_elements('points'), key=str)

        def get_coeffs(frame):
            coeffs = {pt: 0.0 for pt in points}
            for term, coeff in frame.terms.items():
                if len(term.components) == 1 and not isinstance(term.components[0], tuple):
                    pt = term.components[0]
                    coeffs[pt] += coeff
            return coeffs

        coeff1 = get_coeffs(frame1)
        coeff2 = get_coeffs(frame2)

        terms_counter = Counter()

        # Векторная часть (всегда)
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                pi, pj = points[i], points[j]
                vi = coeff1.get(pi, 0)
                vj = coeff1.get(pj, 0)
                wi = coeff2.get(pi, 0)
                wj = coeff2.get(pj, 0)
                coeff = - (vi * wj + vj * wi) / 2
                coeff = Polyform._round_value(coeff, tolerance)
                if abs(coeff) > tolerance:
                    form = quForm([(pi, pj)])
                    terms_counter[form] += coeff

        # Точечная часть
        if abs(mult1) > tolerance or abs(mult2) > tolerance:
            # Вычисляем f1 и f2 как линейные комбинации quForm([point])
            def point_part(coeffs, mult):
                part = Counter()
                for pt, c in coeffs.items():
                    if abs(c) > tolerance:
                        form = quForm([pt])
                        part[form] += c * mult
                return part
            f1 = point_part(coeff1, mult2)
            f2 = point_part(coeff2, mult1)
            total_point = Counter()
            for form, c in f1.items():
                total_point[form] += c
            for form, c in f2.items():
                total_point[form] += c
            for form, c in total_point.items():
                coeff = Polyform._round_value(c / 2, tolerance)
                if abs(coeff) > tolerance:
                    terms_counter[form] += coeff

        return Polyform(terms_counter)

    def to_matrix(self):
        """
        Преобразует полиформу 1-го порядка в матрицу.
        Для термов quForm([(a,b)]) добавляет внедиагональные элементы,
        для термов quForm([a]) — диагональные.
        Возвращает (matrix, index), где index — словарь {элемент: индекс}.
        """
        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")

        # Собираем все элементы (точки или векторы) из полиформы
        elements = sorted(self.get_elements('all'), key=str)
        n = len(elements)
        idx = {elem: i for i, elem in enumerate(elements)}
        M = np.zeros((n, n))

        for term, coeff in self.terms.items():
            if len(term.components) != 1:
                raise ValueError("Форма должна состоять из одной компоненты")
            comp = term.components[0]
            if isinstance(comp, tuple) and len(comp) == 2: # Граничная форма (b - a)
                a, b = comp
                if a not in idx or b not in idx: continue
                i, j = idx[a], idx[b]
                M[i, j] -= coeff
                M[j, i] -= coeff
                M[i, i] += coeff
                M[j, j] += coeff
            else: # Атомарная форма (точка или вектор)
                if comp not in idx: continue
                i = idx[comp]
                M[i, i] += coeff

        return M, idx

    # --- Диагонализация - приведение к базису из собственных векторов ---

    def to_eigenbasis(self, tolerance=1e-10):
        """
        Преобразует однородную полиформу 1-го порядка в диагональную форму в базисе собственных векторов.
        Возвращает Polyform: Σ λ_i * quForm([v_i]), где v_i — собственный вектор (как объект Vector),
        λ_i — собственное число.
        """
        # import numpy as np
        # from objects.element import Vector
        # from .polysimplex import Polysimplex
        from _core.basis import Basis

        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")

        # Получаем матрицу (лапласиан) и отображение элементов->индексы
        M, idx = self.to_matrix()   # M — матрица, idx — {element: index}
        n = M.shape[0]
        if n == 0: return Polyform.zero(), Basis([])

        # Вычисляем собственные значения и векторы (отсортированы по возрастанию)
        eig_vals, eig_vecs = np.linalg.eigh(M)

        # Словарь для новой полиформы и список для базиса
        new_terms = {}
        basis_elements = []

        for k in range(n):
            ev = eig_vals[k]
            vec = eig_vecs[:, k]

            # Нормируем вектор (если уже нормирован, то норма = 1, но для нулевого собственного числа он может быть не нормирован)
            norm = np.linalg.norm(vec)
            if norm > tolerance: vec = vec / norm
            else: continue # Если вектор нулевой (practically), пропускаем? Не должно быть.

            # Строим фрейм (Polysimplex) для собственного вектора в исходном базисе
            frame_terms = {}
            for elem, i in idx.items():
                coeff = vec[i]
                if abs(coeff) < tolerance: continue
                frame_terms[elem] = coeff
            frame = Polysimplex(frame_terms)

            # Создаём элемент с этим фреймом и добавляем в базис
            mult = frame.multiplicity()
            elem = Element(name='v' + str(k), frame=frame, multiplicity=mult)
            basis_elements.append(elem)

            # Добавляем в полиформу (если собственное число не нулевое), коэффициент = собственное число.
            if abs(ev) > tolerance:
                qform = quForm([elem])
                new_terms[qform] = new_terms.get(qform, 0.0) + ev

        # Создаём ортонормированный базис
        basis = Basis(basis_elements, is_orthonormal = True)
        # Возвращаем полиформу и базис
        return Polyform(new_terms), basis

    def expand(self, index=0):
        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Expand применим только к однородным полиформам 1-го порядка")
        result = Polyform.zero()
        for term, coeff in self.terms.items():
            # term — quForm([el])
            el = term.components[0]   # это Element
            if not isinstance(el, Element):
                raise ValueError("Терм должен содержать элемент")
            frame = el[index]   # Polysimplex
            bilinear = Polyform.scalar_product(frame, frame)  # получаем полиформу нормы элемента
            result += bilinear * coeff
        return result

    def exp(self, tolerance=1e-10):
        # Для полиформ: экспонента = Product (1 + term) для всех термов
        result = Polyform.one()
        for term, coeff in self.terms.items():
            # term — quForm, coeff — коэффициент
            factor = Polyform({quForm.one(): 1, term: coeff})
            result = result * factor
        return result