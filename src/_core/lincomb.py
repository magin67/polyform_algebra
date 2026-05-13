# @title Линейная комбинация

# from typing import Dict, Union, List, Set, Tuple, Any
from abc import ABC, abstractmethod
from collections import Counter
import sympy as sp

from _core.monoid import Monoid

class LinearCombination(ABC):
    @staticmethod
    def _round_value(x, tol=1e-10): # Округляет число x с заданным допуском: обнуляет если |x|<tol, иначе округляет до целого, если близко
        if abs(x) < tol: return 0.0
        r = round(x)
        if abs(r - x) < tol: return float(r)
        return x

    def __init__(self, data=None, term_type=None):
        self.terms = Counter()
        self.term_type = term_type
        if data is None:
            self._update_cached_grades()
            return
        # Инициализация из словаря, списка пар, или другого объекта того же класса
        if isinstance(data, dict):
            for term, coeff in data.items(): self._add_term(term, coeff)
        elif isinstance(data, (list, tuple)):
            for term, coeff in data: self._add_term(term, coeff)
        elif isinstance(data, LinearCombination):
            self.terms = data.terms.copy()
            self.term_type = data.term_type
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        self._finalize()

    def _finalize(self, tolerance=1e-10): # Очищает нулевые термы, округляет коэффициенты и обновляет кэшированные грейды.
        # Удаляем термы с нулевым коэффициентом
        to_del = [t for t, c in self.terms.items() if c == 0]
        for t in to_del: del self.terms[t]
        # Округляем оставшиеся коэффициенты
        for t, c in list(self.terms.items()):
            rounded = self._round_value(c, tolerance)
            if rounded == 0: del self.terms[t]
            else: self.terms[t] = rounded
        self._update_cached_grades() # Обновляем кэшированные грейды

    def _create_empty(self): return self.__class__()

    def one(self): return self.__class__({self.term_type.one(): 1})

    def is_one(self): return len(self.terms) == 1 and list(self.terms.values())[0] == 1 and list(self.terms.keys())[0].is_one()

    def is_zero(self): return len(self.terms) == 0

    def is_zero_with_tolerance(self, tolerance=1e-10):
        if not self.terms: return True
        for coeff in self.terms.values():
            if abs(coeff) > tolerance: return False
        return True

    def _add_term(self, term, coeff):
        if coeff == 0: return
        if self.term_type is None:
            self.term_type = type(term)
        elif not isinstance(term, self.term_type):
            raise TypeError(f"Term {term} is not of type {self.term_type}")
        self.terms[term] += coeff

    def _update_cached_grades(self):
        if not self.terms:
            self._min_grade = -1
            self._max_grade = -1
            self._kind = 'empty'
            self._size = 0
        else:
            grades = [t.grade() for t in self.terms]
            self._min_grade = min(grades)
            self._max_grade = max(grades)
            self._kind = 'homogeneous' if self._min_grade == self._max_grade else 'mixed'
            self._size = len(self.terms)

    def size(self): return self._size

    def min_grade(self): return self._min_grade
    def max_grade(self): return self._max_grade
    def grade(self): return self.max_grade()   # для совместимости со старым кодом

    def is_homogeneous(self): return self._kind == 'homogeneous'

    def sizes(self): # Список размеров по грейдам от 0 до max_grade.
        if self.is_zero(): return []
        mg = self.max_grade()
        return [len([t for t in self.terms if t.grade() == g]) for g in range(mg+1)]

    def __add__(self, other):
        # Если other — число, преобразуем его в комбинацию с единичным термом
        if isinstance(other, (int, float)):
            other = self.__class__({self.term_type.one(): other})
        # Если other — одиночный терм, превращаем в комбинацию
        elif isinstance(other, self.term_type):
            other = self.__class__({other: 1})
        # Если other — другой тип LinearCombination, но с другим term_type — ошибка
        elif not isinstance(other, LinearCombination):
            raise TypeError(f"Cannot add {type(other)} to {self.__class__.__name__}")
        elif self.term_type != other.term_type:
            raise TypeError("Cannot add combinations with different term types")
        # Теперь other — LinearCombination
        new_terms = self.terms.copy()
        for term, coeff in other.terms.items():
            new_terms[term] = new_terms.get(term, 0) + coeff
        new = self._create_empty()
        new.terms = new_terms
        new.term_type = self.term_type
        new._finalize()
        return new

    def __radd__(self, other):
        return self + other

    def __sub__(self, other): return self + (-other)

    def __neg__(self):
        new = self.__class__()
        new.terms = Counter({k: -v for k, v in self.terms.items()})
        new.term_type = self.term_type
        new._update_cached_grades()
        return new

    def _mul_terms(self, other): # Умножение двух линейных комбинаций (общий алгоритм).
        if not isinstance(other, LinearCombination):
            raise TypeError(f"Cannot multiply {self.__class__.__name__} by {type(other)}")
        if self.term_type != other.term_type:
            raise TypeError("Cannot multiply combinations with different term types")
        result_terms = Counter()
        for term1, coeff1 in self.terms.items():
            for term2, coeff2 in other.terms.items():
                prod = term1 * term2
                if prod is not None and not (hasattr(prod, 'is_zero') and prod.is_zero()):
                    result_terms[prod] += coeff1 * coeff2
        result = self._create_empty()
        result.terms = result_terms
        result.term_type = self.term_type
        result._finalize()
        return result

    def __mul__(self, other):
        if isinstance(other, (int, float, sp.Expr)): return self.__rmul__(other)
        if isinstance(other, LinearCombination):
            if self.is_zero() or other.is_zero(): return self.zero()
            if self.is_one(): return other
            if other.is_one(): return self
            if self.term_type != other.term_type:
                raise TypeError("Cannot multiply combinations with different term types")
            return self._mul_terms(other)
        if isinstance(other, Monoid): # Проверка совместимости типа
            if self.term_type is not None and not isinstance(other, self.term_type):
                raise TypeError(f"Cannot multiply {self.__class__.__name__} of {self.term_type} by {type(other)}")
            return self * self.__class__({other: 1})
        raise TypeError(f"Cannot multiply {self.__class__.__name__} by {type(other)}")

    def __rmul__(self, scalar):
        if isinstance(scalar, (int, float, sp.Expr)):
            new = self.__class__()
            new.terms = Counter({k: v * scalar for k, v in self.terms.items()})
            new.term_type = self.term_type   # добавлено
            new._finalize()
            new._update_cached_grades()
            return new
        raise TypeError(f"Cannot multiply {self.__class__.__name__} by {type(scalar)}")

    def __truediv__(self, scalar):
        if isinstance(scalar, (int, float, sp.Expr)):
            if scalar == 0:
                raise ZeroDivisionError("Division by zero")
            # Используем умножение на обратное число
            return self * (1 / scalar)
        raise TypeError(f"Cannot divide {self.__class__.__name__} by {type(scalar)}")

    def _create_empty(self): # Создаёт пустой объект того же класса.
        return self.__class__()

    def __eq__(self, other):
        if not isinstance(other, self.__class__): return False
        return self.terms == other.terms

    def __iter__(self):
        return iter(self.terms.items())

    # Методы для работы с элементами (вершинами), общие для всех наследников
    def get_elements(self, tolerance=None):
        elements = set()
        for term, coeff in self.terms.items():
            if tolerance is not None and abs(coeff) <= tolerance: continue
            elements.update(term.get_elements())
        return elements

    def extract_grade(self, grade): # Возвращает новую линейную комбинацию того же типа, содержащую только слагаемые данного грейда.
        new_terms = {t: c for t, c in self.terms.items() if t.grade() == grade}
        return self.__class__(new_terms)

    def leading_coefficient(self): # коэффициент формы максимального грейда
        if self.is_zero(): return 0
        mg = self.max_grade()
        coeffs = [coeff for term, coeff in self.terms.items() if term.grade() == mg]
        if len(coeffs) == 1: return coeffs[0]
        return coeffs

    def leading_terms(self): # терм (или термы) максимального грейда
        if self.is_zero(): return []
        mg = self.max_grade()
        return [term for term in self.terms if term.grade() == mg]

    def multiplicity(self): # Кратность комбинации - сумма для всех термов: коэффициент * multiplicity(term)
        total = 0.0
        for term, coeff in self.terms.items():
            mult = getattr(term, 'multiplicity', 0)
            total += coeff * mult
        return total

    def sum_coeff(self, grade = None): # Сумма коэффициентов всех термов заданного грейда
        if grade == -1: return 0
        comb = self if grade is None else self.extract_grade(grade)
        return sum(comb.terms.values())

    def terms_lists(self): # Кортеж (список_термов, список_коэффициентов)
        terms = []
        coeffs = []
        for term, coeff in self.terms.items():
            terms.append(term)
            coeffs.append(coeff)
        return terms, coeffs

    # Вспомогательные методы для форматирования (могут быть переопределены)
    def _format_coeff(self, coeff):
        if isinstance(coeff, float):
            s = f"{coeff:.10f}".rstrip('0').rstrip('.')
            return s if s != '-0' else '0'
        return str(coeff)

    def _term_str(self, term, use_repr=False): # Строковое представление терма (без коэффициента)
        return repr(term) if use_repr else str(term)

    def _group_by_grade(self):
        groups = {}
        for term, coeff in self.terms.items():
            g = term.grade()
            groups.setdefault(g, []).append((term, coeff))
        return dict(sorted(groups.items()))

    def _format_terms_from_list(self, term_list, use_repr=False):
        parts = []
        for term, coeff in term_list:
            term_str = self._term_str(term, use_repr)
            if coeff == 1:
                parts.append(term_str)
            elif coeff == -1:
                parts.append(f"-{term_str}")
            else:
                coeff_str = str(coeff)
                if isinstance(coeff, float):
                    coeff_str = f"{coeff:.10f}".rstrip('0').rstrip('.')
                parts.append(f"{coeff_str}*{term_str}")
        if not parts: return "0"
        result = parts[0]
        for part in parts[1:]:
            if part[0] == '-': result += f" - {part[1:]}"
            else: result += f" + {part}"
        return result

    def __str__(self):
        if self.is_zero(): return "0"
        groups = self._group_by_grade()
        if len(groups) == 1:  # однородная
            return self._format_terms_from_list(list(self.terms.items()), use_repr=False)
        lines = []
        for g, terms in groups.items():
            line = f"Grade {g}: {self._format_terms_from_list(terms, use_repr=False)}"
            lines.append(line)
        return "\n".join(lines)

    def __repr__(self):
        if self.is_zero():
            return f"{self.__class__.__name__}()"
        return f"{self.__class__.__name__}({self._format_terms_from_list(list(self.terms.items()), use_repr=True)})"

    # Дополнительная математика
    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Only non-negative integer exponents are supported")
        if exponent == 0: return self.one()
        result = self.one()
        base = self
        while exponent:
            if exponent & 1: result = result * base
            base = base * base
            exponent >>= 1
        return result

    def exp(self, max_terms=None):
        result = self.zero()
        term = self.one()
        k = 0
        while True:
            term = term * self
            if term.is_zero(): break
            k += 1
            term = term * (1/k)
            result = result + term
            if max_terms is not None and k >= max_terms: break
        return result

