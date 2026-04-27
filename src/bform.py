# @title Форма границ - косвенное задание

# from src import Monoid, quForm

from .monoid import Monoid
from .quform import quForm
# from .polyform import Polyform


class bForm(Monoid):
    def __init__(self, vectors, dual=False):
        # Удаляем дубликаты векторов (по id)
        seen = set()
        unique = []
        for v in vectors:
            if id(v) not in seen:
                seen.add(id(v))
                unique.append(v)
        super().__init__(tuple(unique), dual)
        self.multiplicity = 0   # важная установка

    def canonical(self): # Сортируем векторы по их строковому представлению (или хешу)
        sorted_vectors = sorted(self.components, key=lambda v: str(v))
        return bForm(sorted_vectors)

    def __hash__(self):
        return hash(self.canonical().components)

    def __eq__(self, other):
        if not isinstance(other, bForm): return False
        return self.canonical().components == other.canonical().components

    def __mul__(self, other):
        result = super().__mul__(other)
        if result is not NotImplemented: return result
        # Проверка на общие векторы
        if set(self.components) & set(other.components): return self.zero()
        # Объединение
        new_vectors = list(self.components) + list(other.components)
        return bForm(new_vectors, dual=self.dual ^ other.dual)

    def _poly_class(self):
        from .poly_bf import PolyBForm
        return PolyBForm

    def expand(self): # Преобразует bForm в Polyform, перемножая to_polyform() каждого вектора.
        from .polyform import Polyform
        result = Polyform(quForm.one())
        for v in self.components:
            result = result * v.to_polyform()
        return result

    def __repr__(self):
        if self.is_zero(): return "bForm()"
        inner = ", ".join(v.name if hasattr(v, 'name') and v.name else repr(v) for v in self.components)
        return f"bForm([{inner}])"

    def __str__(self):
        if self.is_zero(): return "[0]"
        if self.is_one(): return "[1]"
        inner = ", ".join(v.name if hasattr(v, 'name') and v.name else str(v) for v in self.components)
        return f"[{inner}]"