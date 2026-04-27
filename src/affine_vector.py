# @title Аффинный вектор - линейная комбинация точек с нулевой суммой коэффициентов.

from collections import Counter
import sympy as sp

from .lincomb import LinearCombination
from .simplex import Simplex
from .polysimplex import PolySimplex

class AffineVector(PolySimplex):
    def __init__(self, data=None):
        if data is None:
            super().__init__()
            self._check_sum_zero()
            return

        if isinstance(data, dict):
            prepared = {}
            for key, coeff in data.items():
                self._add_to_prepared(prepared, key, coeff)
            super().__init__(prepared)
        elif isinstance(data, (list, tuple)):
            prepared = {}
            for key, coeff in data:
                self._add_to_prepared(prepared, key, coeff)
            super().__init__(prepared)
        elif isinstance(data, AffineVector):
            super().__init__(data.terms)
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        self._check_sum_zero()
        self.multiplicity = 0

    def __repr__(self):
        if self.is_zero(): return "Vector()"
        # Используем use_repr=False для чистого вывода (без лишних кавычек и имён классов)
        inner = self._format_terms_from_list(list(self.terms.items()), use_repr=False)
        return f"Vector({inner})"

    def __str__(self):
        if self.is_zero(): return "0"
        parts = []
        for term, coeff in self.terms.items():
            coeff_str = self._format_coeff(coeff)
            term_str = str(term)
            if coeff_str == "1":
                parts.append(term_str)
            elif coeff_str == "-1":
                parts.append(f"-{term_str}")
            else:
                parts.append(f"{coeff_str}*{term_str}")
        return " + ".join(parts).replace("+ -", "- ")

    def _add_to_prepared(self, prepared, key, coeff):
        if coeff == 0:
            return
        if isinstance(key, tuple) and len(key) == 2:
            a, b = key
            prepared[Simplex([b])] = prepared.get(Simplex([b]), 0) + coeff
            prepared[Simplex([a])] = prepared.get(Simplex([a]), 0) - coeff
        elif isinstance(key, sp.Basic):
            prepared[Simplex([key])] = prepared.get(Simplex([key]), 0) + coeff
        else:
            raise ValueError(f"Unsupported key: {key}")

    def canonical(self): # Каноническое представление вектора (сортировка вершин, нормализация)
        # Сортируем вершины по строковому представлению
        items = sorted(self.terms.items(), key=lambda x: str(x[0].components[0]))
        # Приводим коэффициенты к рациональным числам (можно и float, но лучше рациональные)
        # Для простоты используем кортеж (вершина, коэффициент)
        return tuple((v, c) for v, c in items)

    def __hash__(self):
        return hash(self.canonical())

    def __eq__(self, other):
        if not isinstance(other, AffineVector):
            return False
        return self.canonical() == other.canonical()

    def _check_sum_zero(self):
        total = sum(coeff for _, coeff in self.terms.items())
        if abs(total) > 1e-10:
            raise ValueError("Sum of coefficients must be zero for an affine vector")

    def dot(self, other): # Точечное произведение: сумма коэффициентов при одинаковых вершинах
        if not isinstance(other, AffineVector):
            raise TypeError("Can only compute dot product with another AffineVector")
        total = 0
        for term, coeff in self.terms.items():
            total += coeff * other.terms.get(term, 0)
        return total

    def norm(self): # Квадрат евклидовой нормы (сумму квадратов коэффициентов)
        return self.dot(self)

    def bilinear_form(self, other): # Билинейная форма двух векторов - это форма их скалярного произведения
        from .quform import quForm
        from .polyform import Polyform
        points = set()
        for term in self.terms:
            points.update(term.components)
        for term in other.terms:
            points.update(term.components)
        points = sorted(points, key=str)
        coeffs_self = {p: self.terms.get(Simplex([p]), 0) for p in points}
        coeffs_other = {p: other.terms.get(Simplex([p]), 0) for p in points}
        terms = Counter()
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                pi, pj = points[i], points[j]
                vi = coeffs_self.get(pi, 0)
                vj = coeffs_self.get(pj, 0)
                ui = coeffs_other.get(pi, 0)
                uj = coeffs_other.get(pj, 0)
                coeff = - (vi * uj + vj * ui) / 2
                # Округление
                coeff = self._round_value(coeff)   # используем статический метод из LinearCombination
                if coeff == 0: continue
                a, b = sorted([pi, pj], key=str)
                form = quForm([{a, b}])
                if not form.is_zero(): terms[form] += coeff
        # Дополнительное округление накопленных коэффициентов
        for form in list(terms.keys()):
            val = terms[form]
            new_val = self._round_value(val)
            if new_val == 0: del terms[form]
            else: terms[form] = new_val
        return Polyform(terms)

    def __matmul__(self, other): # Оператор @ для удобного вызова bilinear_form
        if not isinstance(other, AffineVector):
            raise TypeError("Can only compute bilinear form with another AffineVector")
        return self.bilinear_form(other)

    def to_polyform(self): # Полиформа квадратной формы вектора
        return self @ self

    def in_basis(self, basis): # Возвращает LinearCombination, представляющую данный вектор в другом базисе.
        from .lincomb import LinearCombination
        return LinearCombination(basis.coordinates_of(self))

    def to_polyform_in_basis(self, basis): # Полиформа вектора в базисе ортогональных векторов
        from .bform import bForm
        from .poly_bf import PolyBForm
        if not basis.is_orthogonal: raise NotImplementedError("Only orthogonal bases are supported")
        lc = self.in_basis(basis)
        terms = {}
        for b, coeff in lc.terms.items():
            if coeff != 0:
                bf = bForm([b])
                terms[bf] = terms.get(bf, 0) + coeff**2
        return PolyBForm(terms)
