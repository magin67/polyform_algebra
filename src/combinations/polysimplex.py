# @title Полиcимплекс – линейная комбинация симплексов

from collections import Counter

from src._core.lincomb import LinearCombination
from src.objects.simplex import Simplex

class Polysimplex(LinearCombination):
    @classmethod
    def zero(cls):
        return cls()

    @classmethod
    def one(cls, term_type=None):
        if term_type is None:
            term_type = Simplex
        return cls({term_type.one(): 1})

    def __init__(self, data=None, term_type=Simplex):
        self.term_type = term_type
        if data is None:
            super().__init__(term_type=term_type)
            return

        def normalize(term):
            if not isinstance(term, term_type):
                term = term_type(term)
            if term.is_zero():
                return None, None
            key, sign_coeff = term, 1

            if term_type.__name__ == 'Simplex':
                cterm = term.canonical()
                key = term_type(cterm.components, sign=1)
                sign_coeff = cterm.sign

            return key, sign_coeff

        prepared = {}
        if isinstance(data, dict):
            items = data.items()
        elif isinstance(data, (list, tuple)):
            items = data
        elif isinstance(data, term_type):
            items = [(data, 1)]
        elif isinstance(data, Polysimplex):
            for term, coeff in data.terms.items():
                key, sign = normalize(term)
                if key is None:
                    continue
                prepared[key] = prepared.get(key, 0) + coeff * sign
            super().__init__(prepared, term_type=term_type)
            return
        else:
            raise TypeError(f"Unsupported type: {type(data)}")

        for term, coeff in items:
            if coeff == 0: continue
            key, sign = normalize(term)
            if key is None: continue
            prepared[key] = prepared.get(key, 0) + coeff * sign

        super().__init__(prepared, term_type=term_type)

    def kind(self):
        if self.is_zero(): return 'atomic'   # или 'zero'
        atomic = boundary = False
        for term in self.terms:
            k = term.kind()
            if k == 'atomic': atomic = True
            elif k == 'boundary': boundary = True
            else: return 'mixed'
        if atomic and boundary: return 'mixed'
        if atomic: return 'atomic'
        if boundary: return 'boundary'
        return 'atomic'  # пустое (хотя пустой не бывает, т.к. is_zero уже обработан)

    def __eq__(self, other):
        if not isinstance(other, Polysimplex): return False
        return (self - other).is_zero()

    def __str__(self):
        if self.is_zero(): return "[0]"
        parts = []
        for term, coeff in self.terms.items():
            # Нормализуем знак для вывода
            if hasattr(term, 'sign') and term.sign == -1:
                coeff = -coeff
                term_str = str(term).lstrip('-')
            else:
                term_str = str(term)
            if coeff == 1: parts.append(term_str)
            elif coeff == -1: parts.append(f"-{term_str}")
            else: parts.append(f"{coeff}*{term_str}")
        result = parts[0]
        for part in parts[1:]:
            if part[0] == '-':
                result += f" - {part[1:]}"
            else:
                result += f" + {part}"
        return result

    def __repr__(self):
        if self.is_zero():
            return "Polysimplex()"
        return f"Polysimplex({self.__str__()})"

    def __matmul__(self, other):
        from .polyform import Polyform
        return Polyform.scalar_product(self, other)

    def get_elements(self, kind='all'):
        """
        Возвращает множество элементов (точек или векторов), входящих в полисимплекс.
        kind: 'all' - все элементы,
            'points' - только точки (multiplicity == 1),
            'vectors' - только векторы (multiplicity == 0).
        """
        elements = set()

        def add_if_matches(elem):
            if kind == 'all':
                elements.add(elem)
            elif kind == 'points' and hasattr(elem, 'multiplicity') and elem.multiplicity == 1:
                elements.add(elem)
            elif kind == 'vectors' and hasattr(elem, 'multiplicity') and elem.multiplicity == 0:
                elements.add(elem)

        for term in self.terms:
            for comp in term.components:
                if isinstance(comp, tuple):
                    for elem in comp:
                        add_if_matches(elem)
                else:
                    add_if_matches(comp)
        return elements

    def terms_lists(self, as_atomic=True):
        """
        Возвращает (список_термов, список_коэффициентов).
        Если as_atomic=True и терм является симплексом из одного атома (не кортежа),
        то вместо симплекса возвращается сам атом (точка или вектор).
        """
        terms = []
        coeffs = []
        for term, coeff in self.terms.items():
            if as_atomic and len(term.components) == 1 and not isinstance(term.components[0], tuple):
                terms.append(term.components[0])
            else:
                terms.append(term)
            coeffs.append(coeff)
        return terms, coeffs

    def boundary(self): # Возвращает границу полисимплекса (применяет boundary к каждому симплексу).
        if self.is_zero(): return Polysimplex()
        result_terms = Counter()
        for simplex, coeff in self.terms.items():
            bnd = simplex.boundary()
            if not bnd.is_zero():
                for sub_simplex, sub_coeff in bnd.terms.items():
                    result_terms[sub_simplex] += coeff * sub_coeff
        return Polysimplex(result_terms)

    def is_cycle(self, tolerance=1e-10): # Проверяет, является ли полисимплекс циклом (граница равна нулю).
        return self.boundary().is_zero()

    def extract_max_coeff(self, grade=None, comp_count=None, by_abs=True):
        filtered = [(term, coeff) for term, coeff in self.terms.items()
                    if (grade is None or term.grade() == grade)
                    and (comp_count is None or term.size() == comp_count)]
        if not filtered:
            return self.__class__()
        if by_abs:
            max_val = max(abs(c) for _, c in filtered)
            selected = [(term, coeff) for term, coeff in filtered if abs(coeff) == max_val]
        else:
            max_val = max(c for _, c in filtered)
            selected = [(term, coeff) for term, coeff in filtered if coeff == max_val]
        return self.__class__(dict(selected), term_type=self.term_type)

    def extract_max_coeff_all_grades(self, comp_count=None, by_abs=True):
        max_grade = self.max_grade()
        result_terms = {}
        for g in range(max_grade + 1):
            sub = self.extract_max_coeff(grade=g, comp_count=comp_count, by_abs=by_abs)
            for term, coeff in sub.terms.items():
                result_terms[term] = result_terms.get(term, 0) + coeff
        return self.__class__(result_terms, term_type=self.term_type)

    def extract_diagonal_components(self):
        mg = self.max_grade()
        if mg < 0:
            return self.__class__()
        result_terms = {}
        for term, coeff in self.terms.items():
            k = term.grade()
            comp_cnt = term.size()
            if comp_cnt == mg - k + 1:
                result_terms[term] = result_terms.get(term, 0) + coeff
        return self.__class__(result_terms, term_type=self.term_type)

    def transform(self, from_basis, to_basis, tolerance=1e-10):
        from combinations.basis import Basis
        import numpy as np
        """
        Преобразует полисимплекс из базиса from_basis в базис to_basis.
        Требуется, чтобы базисы были связаны (from_basis.link_to(to_basis)).
        """
        if to_basis not in from_basis.transitions:
            raise ValueError("Базисы не связаны. Вызовите from_basis.link_to(to_basis).")
        M, _ = from_basis.transitions[to_basis]   # M: from -> to
        n = len(from_basis)
        coeffs = np.zeros(n)
        # Собираем коэффициенты при элементах from_basis
        for term, coeff in self.terms.items():
            # Ожидаем, что term — симплекс из одного атомарного элемента
            if len(term.components) != 1 or isinstance(term.components[0], tuple):
                raise TypeError("Терм должен быть атомарным симплексом [element]")
            elem = term.components[0]
            idx = from_basis.index(elem)
            if idx == -1:
                raise ValueError(f"Элемент {elem} не принадлежит базису from_basis")
            coeffs[idx] += coeff
        # Преобразуем коэффициенты
        new_coeffs = M @ coeffs
        # Строим новый полисимплекс
        new_terms = {}
        for j, c in enumerate(new_coeffs):
            if abs(c) > tolerance:
                elem = to_basis[j]
                s = Simplex([elem], sign=1)
                new_terms[s] = c
        return self.__class__(new_terms, term_type=self.term_type)
    
    def to_quadratic_form(self, basis, tolerance=1e-10):
        """
        Преобразует полисимплекс (линейную комбинацию векторов) в полиформу,
        предполагая, что базис ортонормирован. 
        Для каждого терма Simplex([v]) с коэффициентом c добавляет слагаемое c^2 * quForm([v]).
        """
        from .polyform import Polyform
        from objects.quform import quForm

        if not basis.is_orthonormal:
            raise ValueError("Метод требует ортонормированного базиса")
        # Для единицы (пустой полисимплекс) возвращаем единичную полиформу? Пока ноль.
        if self.is_zero(): return Polyform.zero()
        new_terms = {}
        for term, coeff in self.terms.items():
            # Ожидаем, что term — Simplex из одного элемента (вектора)
            if len(term.components) != 1 or isinstance(term.components[0], tuple):
                raise TypeError("Терм должен быть атомарным симплексом (один вектор)")
            vec = term.components[0]
            # Коэффициент в квадрате
            c2 = coeff * coeff
            if abs(c2) < tolerance: continue
            qf = quForm([vec])
            new_terms[qf] = new_terms.get(qf, 0) + c2
        # Округляем коэффициенты
        for qf in list(new_terms.keys()):
            val = new_terms[qf]
            rounded = self._round_value(val, tolerance)
            if rounded == 0:
                del new_terms[qf]
            else:
                new_terms[qf] = rounded
        return Polyform(new_terms)