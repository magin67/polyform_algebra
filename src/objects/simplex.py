# @title Симплекс – упорядоченный набор элементов. Каждый элемент списка уникален

from _core.monoid import Monoid

from _core.monoid import Monoid

class Simplex(Monoid):
    @classmethod
    def zero(cls):
        obj = super().zero()
        obj.sign = 1
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
                else:
                    # Кортеж длины 1: либо точка, либо вектор
                    elem = comp[0]
                    if hasattr(elem, 'multiplicity') and elem.multiplicity == 0:
                        # Вектор: выносим из кортежа
                        elements = [elem]
                    else:
                        # Точка: кортеж эквивалентен единице, просто игнорируем его
                        elements = []   # пустой симплекс = единица?
                        # Но единица — это особый объект, лучше сделать Simplex.one()
                        # Пока для простоты оставим как есть, но потом нужно будет вернуть one.
                        # Временно: если все элементы удалены, создаём единичный симплекс.
                        # Однако в контексте создания Point из [self] мы сюда не попадём, так как self не кортеж.
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

    # --------------------------------------------------------------------------
    # Вспомогательные методы для нормализации
    # --------------------------------------------------------------------------
    @staticmethod
    def _extract_vectors_from_components(components):
        """
        Возвращает новый список компонент, в котором из кортежей удалены векторы,
        и отдельный список вынесенных векторов (которые добавляются в базу).
        Правило: если в кортеже есть векторы, они выносятся в базу.
        Знак: при выносе вектора знак кортежа меняется, если размер кортежа чётный?
        Пока не реализуем строго, просто удаляем векторы.
        """
        new_comps = []
        vectors = []
        for comp in components:
            if isinstance(comp, tuple):
                # Разделяем кортеж на точки и векторы
                points = [e for e in comp if not (hasattr(e, 'multiplicity') and e.multiplicity == 0)]
                vecs = [e for e in comp if hasattr(e, 'multiplicity') and e.multiplicity == 0]
                vectors.extend(vecs)
                if points:
                    new_comps.append(tuple(points))
            else:
                if hasattr(comp, 'multiplicity') and comp.multiplicity == 0:
                    vectors.append(comp)
                else:
                    new_comps.append(comp)
        # Добавляем вынесенные векторы в конец (порядок пока не важен)
        return new_comps + vectors

    @staticmethod
    def _normalize_base(base):
        seen = set()
        for e in base:
            if e in seen:
                return None
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
                    if len(common) >= 2:
                        return None
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
                if len(common) >= 2:
                    return None
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
        # Шаг 0: предварительно выносим векторы из границ (упрощённо)
        # Это нужно, чтобы в границах остались только точки (или вообще ничего)
        # Пока просто разделим на базу и границы, не вынося векторы из границ.
        # Но для корректности лучше сделать отдельный шаг.
        # Однако для простоты оставим пока так.
        base = [c for c in components if not isinstance(c, tuple)]
        boundaries = [c for c in components if isinstance(c, tuple)]

        # 1. Нормализация базы (проверка дубликатов)
        base = self._normalize_base(base)
        if base is None:
            return None

        # 2. Нормализация границ (слияние кортежей)
        boundaries = self._normalize_boundaries(boundaries)
        if boundaries is None:
            return None

        # 3. Свёртка базы и границ
        merged = self._merge_base_and_boundaries(base, boundaries)
        if merged is None:
            return None
        base, boundaries = merged

        # 4. После свёртки может оказаться, что в базе снова появились дубликаты? Проверим.
        base = self._normalize_base(base)
        if base is None:
            return None

        # 5. Также после свёртки границы могли приобрести общие элементы? В теории нет, но перестрахуемся
        boundaries = self._normalize_boundaries(boundaries)
        if boundaries is None:
            return None

        # Возвращаем объединённый список (база + границы)
        return base + boundaries

    # --------------------------------------------------------------------------
    # Остальные методы класса (canonical, __eq__, __hash__, __mul__, boundary, ...)
    # --------------------------------------------------------------------------
    def canonical(self):
        if self.is_zero(): return self
        sign = getattr(self, 'sign', 1)   # запасной вариант

        normalized = self._full_normalize(self.components)
        if normalized is None:
            return Simplex.zero()
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
                if perm[i] > perm[j]:
                    sign = -sign
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
        if result is not NotImplemented:
            return result
        if isinstance(other, Simplex):
            if set(self.components) & set(other.components):
                return self.zero()
            new_components = list(self.components) + list(other.components)
            new_sign = self.sign * other.sign
            return Simplex(new_components, new_sign).canonical()
        return NotImplemented

    def boundary(self):
        from combinations.polysimplex import PolySimplex
        if self.is_zero():
            return PolySimplex()
        # Если есть компонента-кортеж (граница), возвращаем ноль
        for comp in self.components:
            if isinstance(comp, tuple):
                return PolySimplex()
        n = len(self.components)
        if n == 0:
            return PolySimplex()
        result = PolySimplex()
        for i in range(n):
            comp = self.components[i]
            # Если компонент — вектор (multiplicity=0), пропускаем
            if hasattr(comp, 'multiplicity') and comp.multiplicity == 0:
                continue
            rest = [self.components[j] for j in range(n) if j != i]
            rest_simplex = Simplex(rest, sign=1)
            coeff = (-1) ** i
            result += PolySimplex({rest_simplex: coeff})
        return result

    def to_polysimplex(self):
        from combinations.polysimplex import PolySimplex
        if self.is_zero():
            return PolySimplex()
        result = PolySimplex.one()
        for comp in self.components:
            if isinstance(comp, tuple):
                s = Simplex(list(comp), sign=1)
                poly = s.boundary()
            else:
                if hasattr(comp, 'multiplicity') and comp.multiplicity == 0:
                    return PolySimplex()
                poly = PolySimplex({Simplex([comp], sign=1): 1})
            result = result * poly
        # Коррекция знака: если число границ (кортежей) нечётное, умножаем на -1
        num_boundaries = sum(1 for comp in self.components if isinstance(comp, tuple))
        if num_boundaries % 2 == 1:
            result = -result
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


# class Simplex(Monoid):
#     def __init__(self, elements, sign=1, dual=False, multiplicity=None):
#         self.sign = sign
#         # Выполняем нормализацию компонент (свёртка границ с атомарными)
#         normalized = self._normalize_components(elements)
#         if normalized is None:
#             # Нулевой симплекс
#             super().__init__([], multiplicity=-1, dual=dual)
#             self.sign = 1
#             return
#         # Вычисляем multiplicity, если не передан
#         if multiplicity is None:
#             mult = 1
#             for comp in normalized:
#                 # Для кортежа multiplicity = 0
#                 if isinstance(comp, tuple):
#                     mult *= 0
#                 else:
#                     mult *= getattr(comp, 'multiplicity', 1)
#         else:
#             mult = multiplicity
#         super().__init__(tuple(normalized), multiplicity=mult, dual=dual)

#     def _normalize_components(self, components):
#         """Возвращает новый список компонент после свёртки кортежей с атомарными элементами,
#            или None, если симплекс должен быть нулевым."""
#         comps = list(components)
#         changed = True
#         while changed:
#             changed = False
#             # Найдём индексы атомарных и кортежей
#             atomic = [(i, c) for i, c in enumerate(comps) if not isinstance(c, tuple)]
#             tuples = [(i, c) for i, c in enumerate(comps) if isinstance(c, tuple)]
#             for ti, tup in tuples:
#                 for ai, elem in atomic:
#                     if elem in tup:
#                         # Проверяем количество общих элементов между tup и всеми атомарными
#                         common = [e for e in tup if any(e == a for _, a in atomic)]
#                         if len(common) >= 2:
#                             return None
#                         # Ровно один общий
#                         # Удаляем кортеж и атомарный элемент
#                         # Сначала удаляем атомарный (его индекс может быть меньше или больше)
#                         if ti < ai:
#                             comps.pop(ai)
#                             comps.pop(ti)
#                         else:
#                             comps.pop(ti)
#                             comps.pop(ai)
#                         # Добавляем элементы кортежа, кроме общего
#                         new_comps = [e for e in tup if e != elem]
#                         # Вставляем их на место удалённого кортежа (индекс ti стал меньше)
#                         # После удалений, ti указывает на позицию, где был кортеж (если удаляли сначала атомарный, то ti может измениться)
#                         # Упростим: просто добавим в конец списка (порядок не важен для дальнейшей нормализации)
#                         comps.extend(new_comps)
#                         changed = True
#                         break
#                 if changed:
#                     break
#         return comps

#     def canonical(self):
#         if self.is_zero():
#             return self
#         # Сначала нормализуем компоненты (свёртка границ с атомарными)
#         normalized = self._normalize_components(self.components)
#         if normalized is None:
#             return Simplex.zero()
#         # Создаём временный симплекс с нормализованными компонентами
#         temp = Simplex(normalized, sign=self.sign)
#         # Затем выполняем сортировку (как раньше)
#         indexed = list(enumerate(temp.components))
#         indexed.sort(key=lambda x: str(x[1]))
#         perm = [i for i, _ in indexed]
#         sign_perm = self._permutation_sign(perm)
#         new_sign = temp.sign * sign_perm
#         sorted_components = tuple(elem for _, elem in indexed)
#         return Simplex(sorted_components, sign=new_sign)


#     @staticmethod
#     def _component_key(comp):
#         # Для хеширования и сравнения компонент
#         if isinstance(comp, tuple):
#             return ('tuple', comp)
#         else:
#             return ('elem', comp)

#     def merge_components(self):
#         """Выполняет слияние компонент-кортежей по правилам границ (Face). 
#         Возвращает новый симплекс с объединёнными компонентами."""
#         # Копируем компоненты
#         comps = list(self.components)
#         changed = True
#         while changed:
#             changed = False
#             new_comps = []
#             # Проходим по списку, пытаясь найти пару для слияния
#             i = 0
#             while i < len(comps):
#                 current = comps[i]
#                 # Ищем партнёра среди оставшихся
#                 merged = False
#                 for j in range(i+1, len(comps)):
#                     other = comps[j]
#                     # Слияние возможно только если оба — кортежи (границы)
#                     if isinstance(current, tuple) and isinstance(other, tuple):
#                         # Ищем общие элементы
#                         set_cur = set(current)
#                         set_oth = set(other)
#                         common = set_cur & set_oth
#                         if len(common) >= 2:
#                             # более двух общих элементов -> ноль
#                             return Simplex.zero()
#                         if len(common) == 1:
#                             elem = next(iter(common))
#                             # Проверяем, не вектор ли общий элемент
#                             if hasattr(elem, 'multiplicity') and elem.multiplicity == 0:
#                                 # общий элемент - вектор -> ноль
#                                 return Simplex.zero()
#                             # Выполняем слияние по алгоритму
#                             # Циклический сдвиг текущего: общий элемент в конец
#                             idx_cur = current.index(elem)
#                             current_rot = current[idx_cur+1:] + current[:idx_cur] + (elem,)
#                             # Циклический сдвиг другого: общий элемент в начало
#                             idx_oth = other.index(elem)
#                             other_rot = (elem,) + other[idx_oth+1:] + other[:idx_oth]
#                             # Объединяем
#                             merged_comp = current_rot + other_rot[1:]
#                             # Заменяем current на merged_comp, убираем other
#                             current = merged_comp
#                             # Удаляем other
#                             comps.pop(j)
#                             changed = True
#                             merged = True
#                             break
#                 if merged:
#                     # продолжаем с тем же i, но current уже заменён
#                     comps[i] = current
#                 else:
#                     new_comps.append(current)
#                     i += 1
#             comps = new_comps
#         # Удаляем возможные дубликаты (хотя их быть не должно)
#         result_comps = []
#         seen = set()
#         for c in comps:
#             key = self._component_key(c)
#             if key not in seen:
#                 seen.add(key)
#                 result_comps.append(c)
#         return Simplex(result_comps, sign=self.sign)

#     # def canonical(self):
#     #     if self.is_zero():
#     #         return self
#     #     # Сначала выполняем слияние границ
#     #     merged = self.merge_components()
#     #     # Затем сортируем компоненты по строковому представлению, 
#     #     # но с учётом знака перестановки компонент (как в старом canonical)
#     #     indexed = list(enumerate(merged.components))
#     #     # Определяем ключ сортировки: сортируем по строке, но кортежи должны быть сравнимы
#     #     indexed.sort(key=lambda x: str(x[1]))
#     #     perm = [i for i, _ in indexed]
#     #     sign_perm = self._permutation_sign(perm)
#     #     new_sign = merged.sign * sign_perm
#     #     sorted_components = tuple(elem for _, elem in indexed)
#     #     return Simplex(sorted_components, sign=new_sign)

#     @staticmethod
#     def _permutation_sign(perm):
#         sign = 1
#         n = len(perm)
#         for i in range(n):
#             for j in range(i+1, n):
#                 if perm[i] > perm[j]:
#                     sign = -sign
#         return sign

#     # @staticmethod
#     # def _permutation_sign(perm): # +1 для чётной перестановки, -1 для нечётной
#     #     sign = 1
#     #     n = len(perm)
#     #     for i in range(n):
#     #         for j in range(i+1, n):
#     #             if perm[i] > perm[j]: sign = -sign
#     #     return sign

#     # def canonical(self):
#     #     if self.is_zero(): return self
#     #     # Сортируем компоненты по строковому представлению
#     #     indexed = list(enumerate(self.components))
#     #     indexed.sort(key=lambda x: str(x[1]))
#     #     perm = [i for i, _ in indexed]   # перестановка исходных индексов
#     #     sign_perm = self._permutation_sign(perm)
#     #     new_sign = self.sign * sign_perm
#     #     sorted_components = tuple(elem for _, elem in indexed)
#     #     return Simplex(sorted_components, sign=new_sign)

#     def __eq__(self, other):
#         if not isinstance(other, Simplex): return False
#         c1 = self.canonical()
#         c2 = other.canonical()
#         return c1.components == c2.components and c1.sign == c2.sign

#     def __hash__(self):
#         c = self.canonical()
#         return hash((c.components, c.sign))

#     def __neg__(self):
#         return Simplex(self.components, sign=-self.sign)

#     # def __mul__(self, other):
#     #     result = super().__mul__(other)
#     #     if result is not NotImplemented:
#     #         return result
#     #     if isinstance(other, Simplex):
#     #         if set(self.components) & set(other.components):
#     #             return self.zero()
#     #         new_components = list(self.components) + list(other.components)
#     #         new_sign = self.sign * other.sign
#     #         return Simplex(new_components, new_sign).canonical()
#     #     return NotImplemented

#     def __mul__(self, other):
#         result = super().__mul__(other)
#         if result is not NotImplemented:
#             return result
#         if isinstance(other, Simplex):
#             if set(self.components) & set(other.components):
#                 return self.zero()
#             new_components = list(self.components) + list(other.components)
#             new_sign = self.sign * other.sign
#             return Simplex(new_components, new_sign).canonical()
#         return NotImplemented

#     def _poly_class(self):
#         from combinations.polysimplex import PolySimplex
#         return PolySimplex

#     def permute(self, indices):
#         if len(indices) != len(self.components):
#             raise ValueError("Number of indices must match simplex size")
#         if sorted(indices) != list(range(len(self.components))):
#             raise ValueError("Indices must be a permutation of 0..size-1")
#         new_components = [self.components[i] for i in indices]
#         return Simplex(new_components, sign=self.sign)

#     def boundary(self):
#         from combinations.polysimplex import PolySimplex
#         if self.is_zero(): return PolySimplex()
#         n = len(self.components)
#         if n == 0: return PolySimplex()
#         result = PolySimplex()
#         for i in range(n):
#             mult_i = getattr(self.components[i], 'multiplicity', 0)
#             if mult_i == 0: continue
#             rest = [self.components[j] for j in range(n) if j != i]
#             rest_simplex = Simplex(rest, sign=1)
#             coeff = mult_i * ((-1) ** i)
#             term = PolySimplex({rest_simplex: coeff})
#             result += term
#         return result

#     def __repr__(self):
#         sign = "-" if self.sign == -1 else ""
#         return f"{sign}{list(self.components)}"

#     def __str__(self):
#         if self.is_zero():
#             return "[]"
#         if self.is_one():
#             return "[1]"
#         sign = "-" if self.sign == -1 else ""
#         inner = ", ".join(str(c) for c in self.components)
#         return f"{sign}[{inner}]"