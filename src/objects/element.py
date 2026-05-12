
from objects.simplex import Simplex
from combinations.polysimplex import Polysimplex

'''Элементы - точки, векторы и другие объекты 1-го грейда'''
class Element(Simplex):
    def __init__(self, name=None, frame=None, multiplicity=1, data=None):
        self.name = name
        self.user_data = data

        self.frames = []
        self.basis_for_frame = []   # список базисов (или None)
        self.basis_for_frame.append(None) # базис пока неизвестен
        if frame is None:
            # создаём тривиальное представление
            trivial = Polysimplex({Simplex([self]): 1})
            self.frames.append(trivial)
        else:
            if isinstance(frame, Polysimplex):
                self._validate_frame(frame, multiplicity)
                self.frames.append(frame)
            elif isinstance(frame, list):
                for f in frame:
                    self._validate_frame(f, multiplicity)
                self.frames.extend(frame)
            else:
                raise TypeError("frame must be PolySimplex or list thereof")
        super().__init__([self], multiplicity=multiplicity)

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

    # Безопасные хеш и равенство, чтобы избежать рекурсии через Simplex
    def __hash__(self):
        if self.name:
            return hash((self.__class__.__name__, self.name)) #, self.multiplicity
        return hash((self.__class__.__name__, id(self))) #, self.multiplicity

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
        return Polysimplex({self: 1, other: 1})

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
        """
        Возвращает новое представление (Polysimplex) элемента в базисе target_basis,
        используя представление с индексом src_index, которое должно быть выражено в базисе src_basis.
        Если src_basis не указан (None), то считается, что представление уже выражено в target_basis? Нет, так нельзя.
        Требуется, чтобы src_basis и target_basis были связаны через link_to.
        """
        if src_index >= len(self.frames):
            raise IndexError("Нет представления с таким индексом")
        frame = self.frames[src_index]
        # Если src_basis не задан, мы не можем преобразовать
        if src_basis is None:
            raise ValueError("Неизвестен базис исходного представления. Укажите src_basis.")
        # Проверяем, что src_basis и target_basis связаны
        if target_basis not in src_basis.transitions:
            raise ValueError(f"Нет перехода из {src_basis} в {target_basis}. Сначала вызовите link_to.")
        M, _ = src_basis.transitions[target_basis]   # из src в target
        # Преобразуем frame (Polysimplex) из src_basis в target_basis
        new_frame = frame.transform(src_basis, target_basis)
        return new_frame



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
