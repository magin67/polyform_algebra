# @title Полиcимплекс – линейная комбинация симплексов

from collections import Counter

from .._core.lincomb import LinearCombination
from ..objects.simplex import Simplex

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
        """Normalize all supported inputs through one coefficient-accumulation path."""
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
                if cterm.is_zero(): return None, None
                key = term_type.one() if cterm.is_one() else term_type(cterm.components, sign=1)
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
            items = data.terms.items()
        else:
            raise TypeError(f"Unsupported type: {type(data)}")

        for term, coeff in items:
            if coeff == 0: continue
            key, sign = normalize(term)
            if key is None: continue
            prepared[key] = prepared.get(key, 0) + coeff * sign

        super().__init__(prepared, term_type=term_type)

    def _from_terms(self, terms):
        """Construct a result with the current term type; subclasses preserve context."""
        return self.__class__(terms, term_type=self.term_type)

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
        if self.is_zero(): return "Polysimplex()"
        return f"Polysimplex({self.__str__()})"

    def __matmul__(self, other):
        from .polyform import Polyform
        return Polyform.scalar_product(self, other)

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

    def as_coefficient_dict(self, basis_elements):
        """
        Возвращает словарь {элемент: коэффициент} для данного полисимплекса
        относительно заданного списка базисных элементов.
        Предполагается, что термы полисимплекса — симплексы из одного атомарного элемента.
        """
        coeffs = {elem: 0.0 for elem in basis_elements}
        for term, coeff in self.terms.items():
            if len(term.components) == 1 and not isinstance(term.components[0], tuple):
                elem = term.components[0]
                if elem in coeffs:
                    coeffs[elem] += coeff
        return coeffs

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
        """Select all tied maxima after optional grade and component-count filters."""
        filtered = [(term, coeff) for term, coeff in self.terms.items()
                    if (grade is None or term.grade() == grade)
                    and (comp_count is None or term.size() == comp_count)]
        if not filtered:
            return self._from_terms({})
        if by_abs:
            max_val = max(abs(c) for _, c in filtered)
            selected = [(term, coeff) for term, coeff in filtered if abs(coeff) == max_val]
        else:
            max_val = max(c for _, c in filtered)
            selected = [(term, coeff) for term, coeff in filtered if coeff == max_val]
        return self._from_terms(dict(selected))

    def extract_max_coeff_all_grades(self, comp_count=None, by_abs=True):
        """Combine coefficient maxima selected independently within each grade."""
        max_grade = self.max_grade()
        result_terms = {}
        for g in range(max_grade + 1):
            sub = self.extract_max_coeff(grade=g, comp_count=comp_count, by_abs=by_abs)
            for term, coeff in sub.terms.items():
                result_terms[term] = result_terms.get(term, 0) + coeff
        return self._from_terms(result_terms)

    def extract_diagonal_components(self):
        """Retain terms whose size equals max_grade - term_grade + 1."""
        mg = self.max_grade()
        if mg < 0:
            return self._from_terms({})
        result_terms = {}
        for term, coeff in self.terms.items():
            k = term.grade()
            comp_cnt = term.size()
            if comp_cnt == mg - k + 1:
                result_terms[term] = result_terms.get(term, 0) + coeff
        return self._from_terms(result_terms)

    def transform(self, from_basis, to_basis, tolerance=1e-10):
        from .._core.basis import Basis
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
        return self._from_terms(new_terms)
    
    def to_diagonal_quform(self, basis, tolerance=1e-10):
        """Return the diagonal part sum(c_i**2 * [e_i]**2) for x=sum(c_i*e_i).

        Accept zero or atomic grade-one coordinates in the supplied Basis,
        declared orthonormal by the caller. No coordinate transformation or
        frame expansion is performed. The declaration is trusted: this method
        has no metric with which to verify orthogonality or normalization.

        This is not polarization or the full formal square x @ x: mixed terms
        are omitted. The result can replace the square for a metric evaluation
        only when that metric annihilates the mixed terms in this basis. It is
        not a substitute for the square in arbitrary algebraic operations.

        tolerance controls rounding of squared coefficients; the usual Polyform
        constructor rounding also applies. Negative coordinates are permitted.
        """
        from math import isfinite
        from numbers import Real
        from .._core.basis import Basis
        from ..objects.element import Element
        from .polyform import Polyform
        from ..objects.quform import quForm

        if isinstance(self, Polyform) or issubclass(self.term_type, quForm):
            raise TypeError("Expected external coordinates, not quadratic forms")
        if not isinstance(basis, Basis):
            raise TypeError("basis must be a Basis")
        if not basis.is_orthonormal:
            raise ValueError("Basis must be declared orthonormal")
        if any(not isinstance(elem, Element) for elem in basis.elements):
            raise TypeError("Basis must contain atomic Elements")
        if len(set(basis.elements)) != len(basis.elements):
            raise ValueError("Basis elements must be distinct")
        if (isinstance(tolerance, bool) or not isinstance(tolerance, Real) or
                not isfinite(tolerance) or tolerance < 0):
            raise ValueError("Tolerance must be finite and non-negative")
        new_terms = {}
        for term, coeff in self.terms.items():
            if (isinstance(term, quForm) or term.grade() != 1 or
                    len(term.components) != 1 or
                    not isinstance(term.components[0], Element)):
                raise TypeError("Expected atomic grade-one coordinates; expand boundaries first")
            elem = term.components[0]
            if elem not in basis:
                raise ValueError(f"Element {elem} does not belong to the supplied basis")
            if not isinstance(coeff, Real):
                raise TypeError("Coordinates must be real numbers")
            squared = coeff * coeff
            if not isfinite(squared):
                raise ValueError("Squared coordinates must be finite")
            rounded = self._round_value(squared, tolerance)
            if rounded:
                new_terms[quForm([elem])] = rounded
        return Polyform(new_terms)

    def to_quadratic_form(self, basis, tolerance=1e-10):
        """Compatibility alias for to_diagonal_quform; this is not a full square.

        Retains the old signature and delegates all validation and calculation
        to the explicitly named method. Use x @ x for supported polarization.
        """
        return self.to_diagonal_quform(basis, tolerance=tolerance)
