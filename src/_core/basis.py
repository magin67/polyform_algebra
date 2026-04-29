# @title Базис

class Basis:
    def __init__(self, elements, default_multiplicity=None):
        """
        elements: список или словарь {элемент: кратность}
        multiplicities: игнорируется, оставлен для совместимости
        """
        if isinstance(elements, dict):
            self.elements = list(elements.keys())
        else:
            self.elements = list(elements)
        # Проверяем хешируемость
        for elem in self.elements:
            if not hasattr(elem, '__hash__'):
                raise TypeError(f"Element {elem} is not hashable")
        self._index = {elem: i for i, elem in enumerate(self.elements)}
        self.is_orthogonal = False   # по умолчанию

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, idx):
        return self.elements[idx]

    def get_elements(self):
        return self.elements.copy()

    def coordinates_of(self, vector):
        if not self.is_orthogonal:
            raise NotImplementedError("Only orthogonal bases are supported")
        coeffs = {}
        for b in self.elements:
            if not hasattr(b, 'dot'):
                raise TypeError(f"Basis element {b} has no 'dot' method")
            norm2 = b.dot(b)
            if norm2 != 0:
                coeff = b.dot(vector) / norm2
                if coeff != 0:
                    coeffs[b] = coeff
        return coeffs