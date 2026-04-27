# @title Моноид - Это объекты, которые можно умножать друг на друга

from abc import ABC, abstractmethod
import sympy as sp

class Monoid(ABC):
    _zero_instances = {}
    _one_instances = {}

    def __init__(self, components, multiplicity=1, dual=False):
        self.components = tuple(components)
        self.multiplicity = multiplicity
        self.dual = dual

    @classmethod
    def zero(cls):
        if cls not in cls._zero_instances:
            obj = cls.__new__(cls)
            obj.components = ()
            obj.multiplicity = -1
            obj.dual = False
            cls._zero_instances[cls] = obj
        return cls._zero_instances[cls]

    @classmethod
    def one(cls):
        if cls not in cls._one_instances:
            obj = cls.__new__(cls)
            obj.components = ()
            obj.multiplicity = -2
            obj.dual = False
            cls._one_instances[cls] = obj
        return cls._one_instances[cls]

    def is_zero(self): return self is self.zero()

    def is_one(self): return self is self.one()

    def size(self): # Количество атомарных элементов (для нуля и единицы — 0)
        if self.is_zero() or self.is_one(): return 0
        return len(self.components)

    def grade(self): # Грейд (по умолчанию равен size, но наследники могут переопределить)
        return self.size()

    def get_elements(self): # Множество атомарных элементов, входящих в моноид
        if self.is_zero() or self.is_one(): return set()
        return set(self.components)

    def __eq__(self, other):
        if not isinstance(other, Monoid): return False
        if self.is_zero() or self.is_one(): return self is other
        return (self.components == other.components and
                self.multiplicity == other.multiplicity and
                self.dual == other.dual)

    def __hash__(self):
        if self.is_zero() or self.is_one(): return id(self)
        return hash((self.components, self.multiplicity, self.dual))

    def __mul__(self, other):
        if not isinstance(other, Monoid): return NotImplemented
        if other == 0: return self.zero()
        if other == 1: return self
        if not isinstance(other, Monoid): return NotImplemented
        if self.is_zero() or other.is_zero(): return self.zero()
        if self.is_one(): return other
        if other.is_one(): return self
        return NotImplemented

    # def __add__(self, other):
    #     if not isinstance(other, Monoid): return NotImplemented
    #     poly_class = self._poly_class()
    #     # Если other уже является линейной комбинацией, делегируем
    #     if isinstance(other, poly_class): return NotImplemented
    #     return poly_class({self: 1, other: 1})

    # def __sub__(self, other):
    #     return self + (-other)

    # def __neg__(self):
    #     poly_class = self._poly_class()
    #     return poly_class({self: -1})

    # def __rmul__(self, scalar):
    #     if isinstance(scalar, (int, float, sp.Expr)):
    #         poly_class = self._poly_class()
    #         return poly_class({self: scalar})
    #     return NotImplemented

    # def _poly_class(self): # Возвращает класс линейных комбинаций для данного моноида.
    #     # По умолчанию — общий LinearCombination. Наследники могут переопределить.
    #     from src.lincomb import LinearCombination
    #     return LinearCombination

    # --- Представление ---
    def __str__(self):
        if self.is_zero(): return "[0]"
        if self.is_one(): return "[1]"
        return self._str_components()

    def _str_components(self): # Строковое представление компонент (по умолчанию — список).
        return str(list(self.components))

    def __repr__(self):
        if self.is_zero(): return "[0]"
        if self.is_one(): return "[1]"
        return f"{list(self.components)}"