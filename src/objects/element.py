
from objects.simplex import Simplex
from _core.lincomb import LinearCombination

class Element(Simplex):
    @classmethod
    def _validate_linear_combination(cls, comb, required_multiplicity):
        # if not isinstance(comb, LinearCombination):
        #     raise TypeError("Data must be a LinearCombination")
        if not comb.is_homogeneous():
            raise ValueError("Linear combination is not homogeneous")
        if comb.grade() != 1:
            raise ValueError("Grade must be 1")
        if abs(comb.multiplicity() - required_multiplicity) > 1e-10:
            raise ValueError(f"Multiplicity must be {required_multiplicity}, got {comb.multiplicity()}")
        return comb # Возвращаем комбинацию для сохранения

    def __init__(self, name='', data=None, multiplicity = 1):
        self.data = data
        if name is not None: self.name = str(name)
        elif data is not None: self.name = str(data)
        else: self.name = None
        super().__init__([self], multiplicity = multiplicity)   # компонента – ссылка на себя

    def __hash__(self):
        return hash((self.name, id(self.data)))

    def __eq__(self, other):
        return isinstance(other, Element) and self.name == other.name and self.data is other.data

    def __add__(self, other):
        if not isinstance(other, Element):
            return NotImplemented
        return LinearCombination({self: 1, other: 1})

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return LinearCombination({self: -1})

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return LinearCombination({self: other})
        return super().__mul__(other)

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return self * other
        return NotImplemented

    def __truediv__(self, scalar):
        if isinstance(scalar, (int, float)):
            return self * (1.0 / scalar)
        return NotImplemented

    def __repr__(self):
        if self.multiplicity == 0: return f"Vector({self.name})"
        elif self.multiplicity == 1: return f"Point({self.name})"
        else: return f"Element({self.name})"

    def __str__(self):
        return self.name if self.name is not None else f"p{id(self)}"


class Point(Element):
    def __init__(self, name=None, data=None):
        if isinstance(data, LinearCombination):
            self._validate_linear_combination(data, 1)
        super().__init__(name=name, data=data, multiplicity=1)

    @staticmethod
    def create_list(names_or_objs):
        if isinstance(names_or_objs, (list, tuple)):
            return [Point(name) for name in names_or_objs]
        else:
            return [Point(name) for name in names_or_objs]  # на самом деле нужно обработать случай одного элемента


class Vector(Element):
    def __init__(self, name=None, data=None):
        if isinstance(data, LinearCombination):
            self._validate_linear_combination(data, 0)
        super().__init__(name=name, data=data, multiplicity=0)

    @staticmethod
    def create_list(names_or_objs):
        if isinstance(names_or_objs, (list, tuple)):
            return [Vector(name) for name in names_or_objs]
        else:
            return [Vector(names_or_objs)]