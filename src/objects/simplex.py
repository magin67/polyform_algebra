# @title Симплекс – упорядоченный набор элементов. Каждый элемент списка уникален

from _core.monoid import Monoid

class Simplex(Monoid):
    def __init__(self, elements, sign=1, dual=False, multiplicity=None):
        # sign должен быть установлен до возможного вызова __hash__
        self.sign = sign
        if multiplicity is None:
            mult = 1
            for e in elements:
                mult *= getattr(e, 'multiplicity', 1)
        else:
            mult = multiplicity
        super().__init__(tuple(elements), multiplicity=mult, dual=dual)
        # Проверка дубликатов после инициализации компонентов
        if len(self.components) != len(set(self.components)):
            raise ValueError("Simplex contains duplicate elements")

    @staticmethod
    def _permutation_sign(perm): # +1 для чётной перестановки, -1 для нечётной
        sign = 1
        n = len(perm)
        for i in range(n):
            for j in range(i+1, n):
                if perm[i] > perm[j]: sign = -sign
        return sign

    def canonical(self):
        if self.is_zero(): return self
        # Сортируем компоненты по строковому представлению
        indexed = list(enumerate(self.components))
        indexed.sort(key=lambda x: str(x[1]))
        perm = [i for i, _ in indexed]   # перестановка исходных индексов
        sign_perm = self._permutation_sign(perm)
        new_sign = self.sign * sign_perm
        sorted_components = tuple(elem for _, elem in indexed)
        return Simplex(sorted_components, sign=new_sign)

    def __eq__(self, other):
        if not isinstance(other, Simplex): return False
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
        if result is not NotImplemented:
            return result
        if isinstance(other, Simplex):
            if set(self.components) & set(other.components):
                return self.zero()
            new_components = list(self.components) + list(other.components)
            new_sign = self.sign * other.sign
            return Simplex(new_components, new_sign).canonical()
        return NotImplemented

    def _poly_class(self):
        from combinations.polysimplex import PolySimplex
        return PolySimplex

    def permute(self, indices):
        if len(indices) != len(self.components):
            raise ValueError("Number of indices must match simplex size")
        if sorted(indices) != list(range(len(self.components))):
            raise ValueError("Indices must be a permutation of 0..size-1")
        new_components = [self.components[i] for i in indices]
        return Simplex(new_components, sign=self.sign)

    def boundary(self):
        from combinations.polysimplex import PolySimplex
        if self.is_zero(): return PolySimplex()
        n = len(self.components)
        if n == 0: return PolySimplex()
        result = PolySimplex()
        for i in range(n):
            mult_i = getattr(self.components[i], 'multiplicity', 0)
            if mult_i == 0: continue
            rest = [self.components[j] for j in range(n) if j != i]
            rest_simplex = Simplex(rest, sign=1)
            coeff = mult_i * ((-1) ** i)
            term = PolySimplex({rest_simplex: coeff})
            result += term
        return result

    def __repr__(self):
        sign = "-" if self.sign == -1 else ""
        return f"{sign}{list(self.components)}"

    def __str__(self):
        if self.is_zero():
            return "[]"
        if self.is_one():
            return "[1]"
        sign = "-" if self.sign == -1 else ""
        inner = ", ".join(str(c) for c in self.components)
        return f"{sign}[{inner}]"