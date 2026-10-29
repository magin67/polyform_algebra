
from .simplex import Simplex
from ..combinations.polysimplex import Polysimplex

"""Elements, points and affine vectors of exterior grade one."""
class Element(Simplex):
    def __init__(self, name=None, frame=None, multiplicity=None, data=None, *,
                 affine_weight=None):
        """Construct an element with an affine weight, defaulting to one.

        ``multiplicity`` is the legacy spelling, retained for compatibility.
        Supplying both names is allowed only when their values agree.
        """
        if affine_weight is None:
            affine_weight = 1 if multiplicity is None else multiplicity
        elif multiplicity is not None and multiplicity != affine_weight:
            raise ValueError("affine_weight and legacy multiplicity must agree")
        super().__init__([self], multiplicity=affine_weight)
        self.name = name
        self.user_data = data

        # Every representation has a corresponding (initially unknown) basis.
        if frame is None:
            self.frames = [Polysimplex({Simplex([self]): 1})]
        else:
            frames = frame if isinstance(frame, list) else [frame]
            for representation in frames:
                self._validate_frame(representation, affine_weight)
            self.frames = list(frames)
        self.basis_for_frame = [None] * len(self.frames)

    @staticmethod
    def _validate_frame(poly, expected_weight):
        if not isinstance(poly, Polysimplex):
            raise TypeError("Frame must be a PolySimplex")
        if not poly.is_homogeneous():
            raise ValueError("Frame is not homogeneous")
        if poly.grade() != 1:
            raise ValueError("Frame grade must be 1")
        if abs(poly.affine_weight() - expected_weight) > 1e-10:
            raise ValueError(f"Frame affine weight must be {expected_weight}, got {poly.affine_weight()}")

    @property
    def affine_weight(self):
        """Value of the affine-weight functional epsilon on this element.

        Legacy storage is shared with ``multiplicity``; there is no second
        independent weight that could diverge from the old API.
        """
        return self.multiplicity

    @affine_weight.setter
    def affine_weight(self, value):
        self.multiplicity = value

    def add_frame(self, frame, basis=None):
        self._validate_frame(frame, self.affine_weight)
        self.frames.append(frame)
        self.basis_for_frame.append(basis)

    def boundary(self):
        """Return the affine weight as a scalar polysimplex: partial x=epsilon(x)."""
        return Polysimplex({Simplex.one(): self.affine_weight})

    def __getitem__(self, idx):
        return self.frames[idx]

    def __hash__(self):
        """Match equality across Element subclasses without recursing into Simplex."""
        if self.name:
            return hash((self.name, self.affine_weight))
        return hash(id(self))

    def __eq__(self, other):
        if not isinstance(other, Element): return False
        if self.name and other.name:
            return self.name == other.name and self.affine_weight == other.affine_weight
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
        # Preserve both summands until their coefficients are accumulated.
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



class Point(Element):
    """An element of affine weight one; no metric or locality is implied."""
    def __init__(self, name=None, frame=None, data=None):
        super().__init__(name=name, frame=frame, affine_weight=1, data=data)

    @classmethod
    def create_list(cls, names):
        return [cls(name) for name in names]

class Vector(Element):
    """An affine vector, an element of affine weight zero."""
    def __init__(self, name=None, frame=None, data=None):
        super().__init__(name=name, frame=frame, affine_weight=0, data=data)

    @classmethod
    def create_list(cls, names):
        return [cls(name) for name in names]
