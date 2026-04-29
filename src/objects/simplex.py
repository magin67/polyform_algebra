# @title Симплекс – упорядоченный набор элементов. Каждый элемент списка уникален

import sympy as sp

from .._core.monoid import Monoid

class Simplex(Monoid):
    def __init__(self, elements, sign=1, dual=False):
        # Убираем дубликаты, сохраняя порядок первого вхождения
        seen = set()
        unique = []
        for e in elements:
            if e not in seen:
                seen.add(e)
                unique.append(e)
        super().__init__(unique)
        self.sign = sign

    def canonical(self):
        if self.is_zero():
            return self
        # Сортируем компоненты по строковому представлению
        indexed = list(enumerate(self.components))
        indexed.sort(key=lambda x: str(x[1]))
        perm = [i for i, _ in indexed]
        from sympy.combinatorics.permutations import Permutation
        sign_perm = -1 if Permutation(perm).is_odd else 1
        new_sign = self.sign * sign_perm
        sorted_components = tuple(elem for _, elem in indexed)
        return Simplex(sorted_components, sign=new_sign)

    def __eq__(self, other):
        if not isinstance(other, Simplex):
            return False
        c1 = self.canonical()
        c2 = other.canonical()
        return c1.components == c2.components and c1.sign == c2.sign

    def __hash__(self):
        c = self.canonical()
        return hash((c.components, c.sign))

    def __neg__(self):
        return Simplex(self.components, sign=-self.sign)

    def __mul__(self, other):
        result = super().__mul__(other)
        if result is not NotImplemented: return result
        if isinstance(other, Simplex):
            if set(self.components) & set(other.components): return self.zero()
            new_components = list(self.components) + list(other.components)
            new_sign = self.sign * other.sign
            return Simplex(new_components, new_sign).canonical()
        elif isinstance(other, sp.Basic):
            if other in self.components: return self.zero()
            new_components = list(self.components) + [other]
            return Simplex(new_components, self.sign).canonical()
        else:
            return NotImplemented

    def _poly_class(self):
        from ..combinations.polysimplex import PolySimplex
        return PolySimplex

    def permute(self, indices):
        if len(indices) != len(self.components):
            raise ValueError("Number of indices must match simplex size")
        if sorted(indices) != list(range(len(self.components))):
            raise ValueError("Indices must be a permutation of 0..size-1")
        new_components = [self.components[i] for i in indices]
        return Simplex(new_components, sign=self.sign)

    def boundary(self):
        from ..combinations.polysimplex import PolySimplex
        if self.is_zero(): return self.zero()
        n = len(self.components)
        if n == 0: return self.zero()
        terms = []
        for i in range(n):
            sub = [self.components[j] for j in range(n) if j != i]
            s = Simplex(sub, sign=1)
            terms.append((s, (-1) ** i))   # self.sign не используется
        return PolySimplex(terms)

    def __repr__(self):
        sign = "-" if self.sign == -1 else ""
        return f"{sign}{list(self.components)}"


    def __str__(self):
        if self.is_zero(): return "[]"
        if self.is_one(): return "[1]"
        sign = "-" if self.sign == -1 else ""
        inner = ", ".join(str(c) for c in self.components)   # str, а не repr
        return f"{sign}[{inner}]"

    # def __str__(self):
    #     sign_str = "-" if self.sign == -1 else ""
    #     return f"{sign_str}{list(self.components)}"
