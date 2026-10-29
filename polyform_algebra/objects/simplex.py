# @title Симплекс – упорядоченный набор элементов. Каждый элемент списка уникален

from .._core.monoid import Monoid

class Simplex(Monoid):
    @classmethod
    def zero(cls):
        obj = super().zero()
        obj.sign = 0
        return obj

    @classmethod
    def one(cls):
        obj = super().one()
        obj.sign = 1
        return obj

    def __init__(self, elements, sign=1, dual=False, multiplicity=None):
        self.sign = sign
        # Определяем, нужна ли полная нормализация
        need_normalize = False
        if len(elements) > 1:
            need_normalize = True

        elif len(elements) == 1:
            comp = elements[0]
            if isinstance(comp, tuple):
                if len(comp) > 1:
                    need_normalize = True
                elif len(comp) == 1:
                    elem = comp[0]
                    if hasattr(elem, 'multiplicity') and elem.multiplicity == 0:
                        elements = [elem]   # вектор выносим
                    else:
                        elements = []       # точка -> единица
                else:  # len(comp) == 0
                    # пустой кортеж эквивалентен единице
                    elements = []

        if need_normalize:
            normalized = self._full_normalize(elements)
            if normalized is None:
                # нулевой симплекс
                super().__init__([], multiplicity=-1, dual=dual)
                self.sign = 1
                return
            components = normalized
        else:
            components = list(elements)
        # Вычисление multiplicity
        if multiplicity is None:
            mult = 1
            for comp in components:
                if isinstance(comp, tuple):
                    mult *= 0
                else:
                    mult *= getattr(comp, 'multiplicity', 1)
        else:
            mult = multiplicity
        super().__init__(tuple(components), multiplicity=mult, dual=dual)

    def is_zero(self):
        return not self.components and self.multiplicity == -1

    def is_one(self):
        return not self.components and self.multiplicity != -1

    def kind(self):
        """Возвращает 'atomic', 'boundary' или 'mixed'."""
        if self.is_zero() or self.is_one():
            return 'atomic'   # единицу и ноль условно считаем атомными
        has_tuple = any(isinstance(comp, tuple) for comp in self.components)
        if not has_tuple: return 'atomic'
        # есть хотя бы один кортеж
        if all(isinstance(comp, tuple) for comp in self.components):
            return 'boundary'
        return 'mixed'

    def get_elements(self, kind='all'):
        """
        Return atomic elements, optionally filtered by affine weight.

        ``points`` selects weight one, ``vectors`` weight zero; neither
        selection asserts metric locality or a nonzero numerical norm.
        """
        elements = set()

        def add_if_matches(elem):
            if kind == 'all':
                elements.add(elem)
            elif kind == 'points' and elem.affine_weight == 1:
                elements.add(elem)
            elif kind == 'vectors' and elem.affine_weight == 0:
                elements.add(elem)

        for comp in self.components:
            if isinstance(comp, tuple):
                for elem in comp:
                    add_if_matches(elem)
            else:
                add_if_matches(comp)
        return elements

    # --------------------------------------------------------------------------
    # Вспомогательные методы для нормализации
    # --------------------------------------------------------------------------
    def _full_normalize(self, components):
        """Канонические факторы с градуированными знаками и слиянием по общей вершине."""
        factors = list(components)
        def degree(x):
            return len(x)-1 if isinstance(x, tuple) else 1
        while True:
            # Ориентация каждой границы и порядок факторов нормализуются отдельно.
            for i, factor in enumerate(factors):
                if isinstance(factor, tuple):
                    if len(set(factor)) != len(factor): return None
                    order = sorted(range(len(factor)), key=lambda k: str(factor[k]))
                    self.sign *= self._permutation_sign(order)
                    factors[i] = tuple(factor[k] for k in order)
            order = sorted(range(len(factors)),
                           key=lambda k: (isinstance(factors[k], tuple), str(factors[k])))
            for i, old_i in enumerate(order):
                for old_j in order[i+1:]:
                    if old_i > old_j:
                        self.sign *= (-1)**(degree(factors[old_i])*degree(factors[old_j]))
            factors = [factors[k] for k in order]
            merged = False
            for i, left in enumerate(factors):
                for j in range(i+1, len(factors)):
                    right = factors[j]
                    left_tuple = isinstance(left, tuple)
                    right_tuple = isinstance(right, tuple)
                    if not left_tuple and not right_tuple:
                        if left == right: return None
                        continue
                    common = (set(left) if left_tuple else {left}) & set(right)
                    if len(common) >= 2: return None
                    if not common: continue
                    pivot = next(iter(common))
                    # Переносим правый множитель непосредственно за левым.
                    self.sign *= (-1)**(degree(right)*sum(degree(x) for x in factors[i+1:j]))
                    self.sign *= (-1)**right.index(pivot)
                    right_rest = [x for x in right if x != pivot]
                    if left_tuple:
                        self.sign *= (-1)**left.index(pivot)
                        replacement = [tuple([pivot] + [x for x in left if x != pivot] + right_rest)]
                    else:
                        replacement = [pivot] + right_rest
                    factors = factors[:i] + replacement + factors[i+1:j] + factors[j+1:]
                    merged = True
                    break
                if merged: break
            if not merged: return factors

    def canonical(self):
        if self.is_zero() or self.is_one(): return self
        return Simplex(self.components, sign=self.sign, multiplicity=self.multiplicity)

    @staticmethod
    def _permutation_sign(perm):
        sign = 1
        n = len(perm)
        for i in range(n):
            for j in range(i+1, n):
                if perm[i] > perm[j]: sign = -sign
        return sign

    def __eq__(self, other):
        if not isinstance(other, Simplex):
            return False
        if self.is_zero() or other.is_zero():
            return self.is_zero() and other.is_zero()
        c1 = self.canonical()
        c2 = other.canonical()
        return c1.components == c2.components and c1.sign == c2.sign

    def __hash__(self):
        if self.is_zero(): return hash(("Simplex.zero",))
        c = self.canonical()
        return hash((c.components, c.sign, c.is_zero()))

    def __add__(self, other):
        from ..combinations.polysimplex import Polysimplex
        if isinstance(other, Simplex):
            return self.as_polysimplex() + other.as_polysimplex()
        if isinstance(other, Polysimplex):
            return self.as_polysimplex() + other
        return NotImplemented

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return self.as_polysimplex().__neg__()

    def __mul__(self, other):
        from ..combinations.polysimplex import Polysimplex
        if isinstance(other, (int, float)):
            return self.as_polysimplex() * other
        if isinstance(other, Simplex):
            if self.is_zero() or other.is_zero(): return Simplex.zero()
            if self.is_one():
                return Simplex(other.components, sign=self.sign * other.sign,
                               multiplicity=other.multiplicity)
            if other.is_one():
                return Simplex(self.components, sign=self.sign * other.sign,
                               multiplicity=self.multiplicity)
            new_components = list(self.components) + list(other.components)
            new_sign = self.sign * other.sign
            return Simplex(new_components, new_sign).canonical()
        if isinstance(other, Polysimplex):
            return self.as_polysimplex() * other
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __matmul__(self, other):
        from ..combinations.polyform import Polyform
        return Polyform.scalar_product(self.as_polysimplex(), other)

    def __str__(self):
        if self.is_zero(): return "[]"
        if self.is_one(): return "[1]"
        sign = "-" if self.sign == -1 else ""
        inner = ", ".join(str(c) for c in self.components)
        return f"{sign}[{inner}]"

    def __repr__(self):
        return self.__str__()

    def grade(self):
        """Грейд симплекса: сумма грейдов компонент (атом → 1, кортеж → len-1)."""
        if self.is_zero(): return -1   # или None, но оставим -1
        if self.is_one(): return 0
        total = 0
        for comp in self.components:
            if isinstance(comp, tuple): total += len(comp) - 1
            else: total += 1
        return total

    def size(self):
        """Размер симплекса: общее количество атомарных элементов (включая элементы внутри кортежей)."""
        if self.is_zero(): return 0
        total = 0
        for comp in self.components:
            if isinstance(comp, tuple): total += len(comp)
            else: total += 1
        return total

    def boundary(self):
        """Apply graded Leibniz contraction with the affine-weight functional."""
        from ..combinations.polysimplex import Polysimplex
        if self.is_zero() or self.is_one():
            return Polysimplex.zero()
        result = Polysimplex.zero()
        prefix_grade = 0
        for i, comp in enumerate(self.components):
            if isinstance(comp, tuple):
                prefix_grade += len(comp) - 1
                continue
            weight = getattr(comp, 'affine_weight', 1)
            if weight:
                rest = self.components[:i] + self.components[i+1:]
                term = Simplex(rest) if rest else Simplex.one()
                result += Polysimplex([(term, self.sign * (-1)**prefix_grade * weight)])
            prefix_grade += 1
        return result

    def as_polysimplex(self):
        from ..combinations.polysimplex import Polysimplex
        return Polysimplex({self: 1})

    def to_polysimplex(self):
        """Раскрывает кортежи-границы, сохраняя ориентацию и атомарные векторы."""
        from ..combinations.polysimplex import Polysimplex
        if self.is_zero():
            return Polysimplex.zero()
        result = Polysimplex.one()
        for comp in self.components:
            if isinstance(comp, tuple):
                poly = Simplex(list(comp)).boundary()
            else:
                poly = Polysimplex({Simplex([comp]): 1})
            result = result * poly
        return self.sign * result

    def to_polyform(self, index=0):
        from ..combinations.polyform import Polyform
        if self.kind() != 'atomic':
            raise ValueError("Метод to_polyform применим только к атомным симплексам")
        if self.is_one(): return Polyform.one()
        if self.is_zero(): return Polyform.zero()
        # Начинаем с единичной полиформы
        result = Polyform.one()
        for comp in self.components:
            # Каждый компонент — элемент (Point, Vector)
            # Получаем его представление (фрейм) по индексу
            try:
                frame = comp[index]   # Polysimplex
            except (IndexError, AttributeError):
                raise RuntimeError(f"Элемент {comp} не имеет представления с индексом {index}")
            # Преобразуем полисимплекс в полиформу: frame @ frame
            # Поскольку frame — это Polysimplex (линейная комбинация точек), его квадрат даёт полиформу
            poly = frame @ frame   # используем __matmul__, который реализован в Polysimplex
            # Умножаем на текущий результат
            result = result * poly
        return result
