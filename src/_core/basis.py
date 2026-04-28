# @title Базис

# from src import Monoid, AffineVector

from .monoid import Monoid

class Basis:
    _default = None # базис по умолчанию

    @classmethod
    def set_default(cls, basis):
        if not isinstance(basis, cls): raise TypeError("Argument must be a Basis instance")
        cls._default = basis

    @classmethod
    def get_default(cls):
        if cls._default is None: raise RuntimeError("Default basis not set")
        return cls._default

    @classmethod
    def has_default(cls):
        return cls._default is not None

    def __init__(self, data, default_multiplicity=1, is_orthogonal=False):
        """
        data: словарь {элемент: кратность} или список элементов (кратность default_multiplicity)
        default_multiplicity: кратность для элементов, заданных списком
        """
        self.elements = []
        self.index = {}
        if isinstance(data, dict): items = data.items()
        else: items = [(elem, default_multiplicity) for elem in data]

        for elem, mult in items:
            monoid = self._to_monoid(elem, mult)
            if monoid in self.index: raise ValueError(f"Duplicate element {monoid} in basis")
            self.elements.append(monoid)
            self.index[monoid] = len(self.elements) - 1
        self.is_orthogonal = is_orthogonal

    def _to_monoid(self, elem, mult):
        # print(f"_to_monoid: {elem}, type={type(elem)}, mult={mult}")
        if isinstance(elem, Monoid):
            if elem.multiplicity != mult: raise ValueError(f"Multiplicity mismatch for {elem}: expected {mult}, got {elem.multiplicity}")
            return elem
        return Monoid([elem], multiplicity=mult) # Создаём моноид с заданной кратностью

    def get_elements(self):
        return self.elements.copy()

    def multiplicity(self, elem):
        return elem.multiplicity   # теперь кратность хранится в самом моноиде

    def size(self):
        return len(self.elements)

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, idx):
        return self.elements[idx]

    def __repr__(self):
        return f"Basis({ {elem: elem.multiplicity for elem in self.elements} })"

    def __str__(self):
        return repr(self)

    def coordinates_of(self, vector): # Возвращает координаты `vector` в базисе. Предполагается, что базис ортогональный.
        if not self.is_orthogonal: raise NotImplementedError("Only orthogonal bases are supported")
        coeffs = {}
        for b in self.get_elements():
            norm2 = b.dot(b)
            if norm2 != 0:
                coeff = b.dot(vector) / norm2
                if coeff != 0:
                    coeffs[b] = coeff
        return coeffs
