
from .simplex import Simplex
from ..combinations.polysimplex import Polysimplex

'''Элементы - точки, векторы и другие объекты 1-го грейда'''
class Element(Simplex):
    def __init__(self, name=None, frame=None, multiplicity=1, data=None):
        super().__init__([self], multiplicity=multiplicity)
        self.name = name
        self.user_data = data

        # Every representation has a corresponding (initially unknown) basis.
        if frame is None:
            self.frames = [Polysimplex({Simplex([self]): 1})]
        else:
            frames = frame if isinstance(frame, list) else [frame]
            for representation in frames:
                self._validate_frame(representation, multiplicity)
            self.frames = list(frames)
        self.basis_for_frame = [None] * len(self.frames)

    @staticmethod
    def _validate_frame(poly, expected_multiplicity):
        if not isinstance(poly, Polysimplex):
            raise TypeError("Frame must be a PolySimplex")
        if not poly.is_homogeneous():
            raise ValueError("Frame is not homogeneous")
        if poly.grade() != 1:
            raise ValueError("Frame grade must be 1")
        if abs(poly.multiplicity() - expected_multiplicity) > 1e-10:
            raise ValueError(f"Frame multiplicity must be {expected_multiplicity}, got {poly.multiplicity()}")

    def add_frame(self, frame, basis=None):
        self._validate_frame(frame, self.multiplicity)
        self.frames.append(frame)
        self.basis_for_frame.append(basis)

    def boundary(self): # Возвращает полисимплекс, представляющий кратность элемента
        return Polysimplex({Simplex.one(): self.multiplicity})

    def __getitem__(self, idx):
        return self.frames[idx]

    def __hash__(self):
        """Match equality across Element subclasses without recursing into Simplex."""
        if self.name:
            return hash((self.name, self.multiplicity))
        return hash(id(self))

    def __eq__(self, other):
        if not isinstance(other, Element): return False
        if self.name and other.name:
            return self.name == other.name and self.multiplicity == other.multiplicity
        return self is other

    def __str__(self):
        if self.name: return self.name
        return f"{self.__class__.__name__}()"

    def __repr__(self):
        return self.__str__()

    # --- Арифметика ---
    def __add__(self, other):
        if not isinstance(other, Element):
            return NotImplemented
        # Список сохраняет оба слагаемых до сложения их коэффициентов.
        return Polysimplex([(self, 1), (other, 1)])

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return Polysimplex({self: -1})

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return Polysimplex({self: other})
        return super().__mul__(other)

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return self * other
        return NotImplemented

    def __truediv__(self, scalar):
        if isinstance(scalar, (int, float)):
            return self * (1.0 / scalar)
        return NotImplemented

    def in_basis(self, src_index, src_basis, target_basis):
        """Transform a stored representation between explicitly linked bases.

        The source basis must be supplied; it cannot be inferred from a frame.
        """
        if src_index >= len(self.frames):
            raise IndexError("Нет представления с таким индексом")
        frame = self.frames[src_index]
        # A frame alone does not identify its source coordinate system.
        if src_basis is None:
            raise ValueError("Неизвестен базис исходного представления. Укажите src_basis.")
        # link_to() must establish the transition before it can be used.
        if target_basis not in src_basis.transitions:
            raise ValueError(f"Нет перехода из {src_basis} в {target_basis}. Сначала вызовите link_to.")
        # Delegate coordinate conversion to the common implementation.
        return frame.transform(src_basis, target_basis)



'''Точки - это элементы кратности 1'''
class Point(Element):
    def __init__(self, name=None, frame=None, data=None):
        super().__init__(name=name, frame=frame, multiplicity=1, data=data)

    @classmethod
    def create_list(cls, names):
        return [cls(name) for name in names]

'''Векторы - это элементы кратности 0'''
class Vector(Element):
    def __init__(self, name=None, frame=None, data=None):
        super().__init__(name=name, frame=frame, multiplicity=0, data=data)

    @classmethod
    def create_list(cls, names):
        return [cls(name) for name in names]
