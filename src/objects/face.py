# @title Граница – набор непересекающихся компонент

from _core.monoid import Monoid
from .simplex import Simplex

class Face(Monoid):
    def __init__(self, components, dual=False):
        norm = self._normalize(components)
        if norm is None:
            super().__init__([])
        else:
            super().__init__(norm)

    def grade(self):
        return sum(len(c) - 1 for c in self.components)

    def size(self):
        return sum(len(c) for c in self.components)

    def cardinality(self):
        return len(self.components)

    @staticmethod
    def _normalize(components):
        lst = [list(c) for c in components]
        changed = True
        while changed:
            changed = False
            new_lst = []
            while lst:
                current = lst.pop(0)
                merged = False
                for i, other in enumerate(lst):
                    common = set(current) & set(other)
                    if len(common) >= 2:
                        return None
                    if len(common) == 1:
                        common_elem = common.pop()
                        # Циклический сдвиг current: общий элемент в конец
                        idx_cur = current.index(common_elem)
                        current_rot = current[idx_cur+1:] + current[:idx_cur] + [common_elem]
                        # Циклический сдвиг other: общий элемент в начало
                        idx_oth = other.index(common_elem)
                        other_rot = [common_elem] + other[idx_oth+1:] + other[:idx_oth]
                        merged_comp = current_rot + other_rot[1:]
                        current = merged_comp
                        lst.pop(i)
                        merged = True
                        changed = True
                        break
                if merged:
                    lst.insert(0, current)
                else:
                    new_lst.append(current)
            lst = new_lst
        # Удаляем возможные дубликаты (хотя их быть не должно)
        result = []
        for comp in lst:
            seen = set()
            uniq = []
            for e in comp:
                if e not in seen:
                    seen.add(e)
                    uniq.append(e)
            result.append(tuple(uniq))
        return result

    def get_elements(self):
        elems = set()
        for comp in self.components:
            elems.update(comp)
        return elems

    def __mul__(self, other):
        result = super().__mul__(other)
        if result is not NotImplemented: return result
        if not isinstance(other, Face): return NotImplemented
        all_components = list(self.components) + list(other.components)
        normalized = self._normalize(all_components)
        if normalized is None: return self.zero()
        return Face(normalized)

    def _poly_class(self):
        from ..combinations.polyface import PolyFace
        return PolyFace

    def __repr__(self):
        inner = ", ".join(str(tuple(c)) for c in self.components)
        return f"Face([{inner}])"

    def __str__(self):
        inner = ", ".join(str(tuple(c)) for c in self.components)
        return f"[{inner}]"

    def to_quform(self): # Преобразует границу в квадратичную форму (quForm)
        from .quform import quForm
        return quForm(self)

    def to_polysimplex(self): # Преобразует границу в PolySimplex (линейную комбинацию симплексов)
        if self.is_zero(): return self.zero()
        # Берём первую компоненту
        first = Simplex(list(self.components[0]))
        result = first.boundary()
        for comp in self.components[1:]:
            s = Simplex(list(comp))
            result = result * s.boundary()
        return result