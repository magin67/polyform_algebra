# @title Полиcимплекс – линейная комбинация симплексов

from collections import Counter

from _core.lincomb import LinearCombination
from objects.simplex import Simplex

class Polysimplex(LinearCombination):
    @classmethod
    def one(cls):
        return cls({Simplex.one(): 1})

    @staticmethod
    def add_term(prepared, term, coeff):
        cterm = term.canonical() # Приводим симплекс к каноническому виду
        key = Simplex(cterm.components, sign=1) # Создаём симплекс без знака (sign=1)
        prepared[key] = prepared.get(key, 0) + coeff * cterm.sign # Умножаем коэффициент на знак канонического симплекса

    def __init__(self, data=None):
        if data is None:
            super().__init__(term_type=Simplex)
            return
        prepared = {}
        if isinstance(data, dict):
            for term, coeff in data.items():
                if not isinstance(term, Simplex):
                    if isinstance(term, (list, tuple)): term = Simplex(term)
                    else: raise TypeError(f"Cannot convert {term} to Simplex")
                self.add_term(prepared, term, coeff)
        elif isinstance(data, (list, tuple)):
            for term, coeff in data:
                if not isinstance(term, Simplex):
                    if isinstance(term, (list, tuple)): term = Simplex(term)
                    else: raise TypeError(f"Cannot convert {term} to Simplex")
                self.add_term(prepared, term, coeff)
        elif isinstance(data, Simplex):
            self.add_term(prepared, data, coeff)
        elif isinstance(data, Polysimplex):
            super().__init__(data.terms, term_type=Simplex)
            return
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        super().__init__(prepared, term_type=Simplex)

    def kind(self):
        if self.is_zero(): return 'atomic'   # или 'zero'
        atomic = boundary = False
        for term in self.terms:
            k = term.kind()
            if k == 'atomic': atomic = True
            elif k == 'boundary': boundary = True
            else:  # mixed
                return 'mixed'
        if atomic and boundary:
            return 'mixed'
        if atomic: return 'atomic'
        if boundary: return 'boundary'
        return 'atomic'  # пустое (хотя пустой не бывает, т.к. is_zero уже обработан)

    def __eq__(self, other):
        if not isinstance(other, Polysimplex): return False
        return (self - other).is_zero()

    def __str__(self):
        if self.is_zero(): return "0"
        parts = []
        for term, coeff in self.terms.items():
            # Нормализуем знак для вывода
            if hasattr(term, 'sign') and term.sign == -1:
                coeff = -coeff
                term_str = str(term).lstrip('-')
            else:
                term_str = str(term)
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
        if self.is_zero():
            return "Polysimplex()"
        return f"Polysimplex({self.__str__()})"

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
        return self.__class__(dict(selected))

    def extract_max_coeff_all_grades(self, comp_count=None, by_abs=True):
        max_grade = self.max_grade()
        result_terms = {}
        for g in range(max_grade + 1):
            sub = self.extract_max_coeff(grade=g, comp_count=comp_count, by_abs=by_abs)
            for term, coeff in sub.terms.items():
                result_terms[term] = result_terms.get(term, 0) + coeff
        return self.__class__(result_terms)

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
        return self.__class__(result_terms)