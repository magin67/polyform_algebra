# @title Симплекс – упорядоченный набор элементов. Каждый элемент списка уникален

from src._core.monoid import Monoid

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

    # --------------------------------------------------------------------------
    # Вспомогательные методы для нормализации
    # --------------------------------------------------------------------------
    @staticmethod
    def _normalize_base(base):
        seen = set()
        for e in base:
            if e in seen: return None
            seen.add(e)
        return base

    @staticmethod
    def _normalize_boundaries(boundaries):
        """
        Нормализует список кортежей (границ) по правилам Face:
        - Если две границы имеют >=2 общих элементов -> None (ноль).
        - Если ровно один общий элемент, они сливаются в один кортеж.
        Повторяем, пока возможно.
        Возвращает новый список кортежей или None.
        """
        # Работаем с изменяемыми списками для удобства
        comps = [list(b) for b in boundaries]
        changed = True
        while changed:
            changed = False
            for i in range(len(comps)):
                for j in range(i+1, len(comps)):
                    set_i = set(comps[i])
                    set_j = set(comps[j])
                    common = set_i & set_j
                    if len(common) >= 2: return None
                    if len(common) == 1:
                        elem = next(iter(common))
                        # Циклический сдвиг i: общий элемент в конец
                        idx_i = comps[i].index(elem)
                        rotated_i = comps[i][idx_i+1:] + comps[i][:idx_i] + [elem]
                        # Циклический сдвиг j: общий элемент в начало
                        idx_j = comps[j].index(elem)
                        rotated_j = [elem] + comps[j][idx_j+1:] + comps[j][:idx_j]
                        merged = rotated_i + rotated_j[1:]
                        comps[i] = merged
                        comps.pop(j)
                        changed = True
                        break
                if changed:
                    break
        return [tuple(c) for c in comps]

    @staticmethod
    def _merge_base_and_boundaries(base, boundaries):
        base_list = list(base)
        boundaries_list = list(boundaries)
        changed = True
        while changed:
            changed = False
            for i, b in enumerate(boundaries_list):
                common = [e for e in base_list if e in b]
                if len(common) >= 2: return None
                if len(common) == 1:
                    elem = common[0]
                    # Не удаляем elem из base_list
                    boundaries_list.pop(i)
                    # Добавляем все элементы кортежа, кроме elem (чтобы не дублировать)
                    new_elems = [e for e in b if e != elem]
                    base_list.extend(new_elems)
                    changed = True
                    break
        return base_list, boundaries_list

    def _full_normalize(self, components):
        """
        Полная нормализация списка компонент (смесь атомарных элементов и кортежей).
        Возвращает новый список компонент или None (если симплекс нулевой).
        """
        base = [c for c in components if not isinstance(c, tuple)]
        boundaries = [c for c in components if isinstance(c, tuple)]

        # 1. Нормализация базы (проверка дубликатов)
        base = self._normalize_base(base)
        if base is None: return None

        # 2. Нормализация границ (слияние кортежей)
        boundaries = self._normalize_boundaries(boundaries)
        if boundaries is None: return None

        # 3. Свёртка базы и границ
        merged = self._merge_base_and_boundaries(base, boundaries)
        if merged is None: return None
        base, boundaries = merged

        # 4. После свёртки может оказаться, что в базе снова появились дубликаты? Проверим.
        base = self._normalize_base(base)
        if base is None: return None

        # 5. Также после свёртки границы могли приобрести общие элементы? В теории нет, но перестрахуемся
        boundaries = self._normalize_boundaries(boundaries)
        if boundaries is None: return None

        # Возвращаем объединённый список (база + границы)
        return base + boundaries

    def canonical(self):
        if self.is_zero(): return self
        sign = getattr(self, 'sign', 1)   # запасной вариант

        normalized = self._full_normalize(self.components)
        if normalized is None or normalized == []: return self.zero()

        temp = Simplex(normalized, sign=sign)
        # Сортировка компонент по строковому представлению
        indexed = list(enumerate(temp.components))
        indexed.sort(key=lambda x: str(x[1]))
        perm = [i for i, _ in indexed]
        sign_perm = self._permutation_sign(perm)
        new_sign = temp.sign * sign_perm
        sorted_components = tuple(elem for _, elem in indexed)
        return Simplex(sorted_components, sign=new_sign)

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
            new_components = list(self.components) + list(other.components)
            new_sign = self.sign * other.sign
            return Simplex(new_components, new_sign).canonical()
        return NotImplemented

    def __repr__(self):
        sign = "-" if self.sign == -1 else ""
        return f"{sign}{list(self.components)}"

    def __str__(self):
        if self.is_zero(): return "[]"
        if self.is_one(): return "[1]"
        sign = "-" if self.sign == -1 else ""
        inner = ", ".join(str(c) for c in self.components)
        return f"{sign}[{inner}]"

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
        from src.combinations.polysimplex import Polysimplex
        if self.is_zero():
            return Polysimplex()
        # Если есть компонента-кортеж (граница), возвращаем ноль
        for comp in self.components:
            if isinstance(comp, tuple):
                return Polysimplex()
        n = len(self.components)
        if n == 0:
            return Polysimplex()
        result = Polysimplex()
        for i in range(n):
            comp = self.components[i]
            # Если компонент — вектор (multiplicity=0), пропускаем
            if hasattr(comp, 'multiplicity') and comp.multiplicity == 0:
                continue
            rest = [self.components[j] for j in range(n) if j != i]
            rest_simplex = Simplex(rest, sign=1)
            coeff = (-1) ** i
            result += Polysimplex({rest_simplex: coeff})
        return result

    def to_polysimplex(self):
        from src.combinations.polysimplex import Polysimplex
        if self.is_zero():
            return Polysimplex()
        result = Polysimplex.one()
        for comp in self.components:
            if isinstance(comp, tuple):
                s = Simplex(list(comp), sign=1)
                poly = s.boundary()
            else:
                if hasattr(comp, 'multiplicity') and comp.multiplicity == 0:
                    return Polysimplex()
                poly = Polysimplex({Simplex([comp], sign=1): 1})
            result = result * poly
        # Коррекция знака: если число границ (кортежей) нечётное, умножаем на -1
        num_boundaries = sum(1 for comp in self.components if isinstance(comp, tuple))
        if num_boundaries % 2 == 1: result = -result
        return result

    def to_polyform(self, index=0):
        from src.combinations.polyform import Polyform
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
