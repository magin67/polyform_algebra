# @title Форма – набор явно заданных непересекающихся кортежей. Элементы неупорядочены.

import sympy as sp

from .monoid import Monoid
from .face import Face
# from .polyform import Polyform

class quForm(Monoid):
    def __init__(self, data=None, dual=False):
        if isinstance(data, Face): # Преобразуем границу в форму: берём компоненты как множества
            data = [set(comp) for comp in data.components]
        # data - список множеств (или списков) элементов (символов)
        # Нормализация: объединение пересекающихся множеств, проверка на >=2 общих
        normalized = self._normalize([frozenset(s) for s in data if s])
        # Сортируем компоненты для каноничности
        normalized.sort(key=lambda t: (len(t), str(t[0]) if t else ''))
        # Сохраняем как кортеж кортежей
        super().__init__(normalized)

    @staticmethod
    def _normalize(sets):
        lst = [set(s) for s in sets]   # преобразуем кортежи в множества
        changed = True
        while changed:
            changed = False
            new_lst = []
            while lst:
                current = lst.pop(0)
                merged = False
                for i, other in enumerate(lst):
                    common = current & other
                    if len(common) >= 2:
                        return []
                    if common:
                        current = current | other
                        lst.pop(i)
                        merged = True
                        changed = True
                        break
                if merged:
                    lst.insert(0, current)
                else:
                    new_lst.append(current)
            lst = new_lst
        return [tuple(sorted(c, key=str)) for c in lst]

    def grade(self) -> int: # Возвращает грейд формы = сумма (len(comp) - 1) по компонентам
        if self.is_zero(): return -1 # Ноль не имеет грейда, вернём -1 как индикатор
        if self.is_one(): return 0
        return sum(len(comp) - 1 for comp in self.components)

    def __mul__(self, other):
        result = super().__mul__(other)
        if result is not NotImplemented: return result
        if not isinstance(other, quForm): return NotImplemented
        all_sets = list(self.components) + list(other.components)
        normalized = self._normalize([set(c) for c in all_sets])
        return quForm([set(c) for c in normalized])

    def _poly_class(self):
        from .polyform import Polyform
        return Polyform

    def __str__(self):
        if self.is_zero(): return "[]"
        if self.is_one(): return "[1]"
        comps = ', '.join('{' + ', '.join(str(e) for e in c) + '}' for c in self.components)
        return f"[{comps}]"

    def __repr__(self): return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, quForm): return False
        return self.components == other.components

    def __hash__(self): return hash(self.components)

    def __mul__(self, other: 'quForm') -> 'quForm':
        if self.is_zero() or other.is_zero(): return quForm.zero()
        if self.is_one(): return other
        if other.is_one(): return self

        all_sets = list(self.components) + list(other.components)
        normalized = self._normalize(all_sets)
        if not normalized: return quForm.zero()
        return quForm([set(c) for c in normalized])

    def __rmul__(self, other):
        from .polyform import Polyform
        if isinstance(other, (int, float, sp.Expr)): return Polyform({self: other})
        raise TypeError("Умножение формы возможно только на другую форму или число")

    def get_elements(self): # Множество всех элементов, входящих в форму
        elems = set()
        for comp in self.components:
            elems.update(comp)
        return elems

    def get_vertices(self): # аналог get_elements
        return self.get_elements()
