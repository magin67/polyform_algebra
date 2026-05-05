
from objects.simplex import Simplex
from combinations.polysimplex import PolySimplex

'''Элементы - точки, векторы и другие объекты 1-го грейда'''
class Element(Simplex):
    def __init__(self, name=None, frame=None, multiplicity=1, data=None):
        self.name = name
        self.user_data = data
        self.frames = []
        if frame is None:
            self._auto_frame = True
        else:
            self._auto_frame = False
            if isinstance(frame, PolySimplex):
                self._validate_frame(frame, multiplicity)
                self.frames.append(frame)
            elif isinstance(frame, list):
                for f in frame:
                    self._validate_frame(f, multiplicity)
                self.frames.extend(frame)
            else:
                raise TypeError("frame must be PolySimplex or list thereof")
        super().__init__([self], multiplicity=multiplicity)
        if self._auto_frame:
            trivial = PolySimplex({self: 1})
            self.frames.append(trivial)

    @staticmethod
    def _validate_frame(poly, expected_multiplicity):
        if not isinstance(poly, PolySimplex):
            raise TypeError("Frame must be a PolySimplex")
        if not poly.is_homogeneous():
            raise ValueError("Frame is not homogeneous")
        if poly.grade() != 1:
            raise ValueError("Frame grade must be 1")
        if abs(poly.multiplicity() - expected_multiplicity) > 1e-10:
            raise ValueError(f"Frame multiplicity must be {expected_multiplicity}, got {poly.multiplicity()}")

    def add_frame(self, poly):
        self._validate_frame(poly, self.multiplicity)
        self.frames.append(poly)

    def __getitem__(self, idx):
        return self.frames[idx]

    def boundary(self): # Возвращает полисимплекс, представляющий кратность элемента
        return PolySimplex({Simplex.one(): self.multiplicity})

    # Безопасные хеш и равенство, чтобы избежать рекурсии через Simplex
    def __hash__(self):
        if self.name:
            return hash((self.__class__.__name__, self.name, self.multiplicity))
        return hash((self.__class__.__name__, id(self), self.multiplicity))

    def __eq__(self, other):
        if not isinstance(other, Element):
            return False
        if self.name and other.name:
            return self.name == other.name and self.multiplicity == other.multiplicity
        return self is other

    def __str__(self):
        if self.name:
            return self.name
        return f"{self.__class__.__name__}()"

    def __repr__(self):
        return self.__str__()

    # --- Арифметика ---
    def __add__(self, other):
        if not isinstance(other, Element):
            return NotImplemented
        return PolySimplex({self: 1, other: 1})

    def __sub__(self, other):
        return self + (-other)

    def __neg__(self):
        return PolySimplex({self: -1})

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return PolySimplex({self: other})
        return super().__mul__(other)

    def __rmul__(self, other):
        if isinstance(other, (int, float)):
            return self * other
        return NotImplemented

    def __truediv__(self, scalar):
        if isinstance(scalar, (int, float)):
            return self * (1.0 / scalar)
        return NotImplemented


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
