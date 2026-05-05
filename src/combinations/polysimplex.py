# @title Полиcимплекс – линейная комбинация симплексов

from collections import Counter

from _core.lincomb import LinearCombination
from objects.simplex import Simplex

class PolySimplex(LinearCombination):
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
        elif isinstance(data, PolySimplex):
            super().__init__(data.terms, term_type=Simplex)
            return
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        super().__init__(prepared, term_type=Simplex)

    def __eq__(self, other):
        if not isinstance(other, PolySimplex): return False
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
            return "PolySimplex()"
        return f"PolySimplex({self.__str__()})"

    def boundary(self): # Возвращает границу полисимплекса (применяет boundary к каждому симплексу).
        if self.is_zero(): return PolySimplex()
        result_terms = Counter()
        for simplex, coeff in self.terms.items():
            bnd = simplex.boundary()
            if not bnd.is_zero():
                for sub_simplex, sub_coeff in bnd.terms.items():
                    result_terms[sub_simplex] += coeff * sub_coeff
        return PolySimplex(result_terms)

    def is_cycle(self, tolerance=1e-10): # Проверяет, является ли полисимплекс циклом (граница равна нулю).
        return self.boundary().is_zero()
