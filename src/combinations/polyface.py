# @title Полиграница – линейная комбинация границ

from .._core.lincomb import LinearCombination
from ..objects.face import Face
from ..combinations.polysimplex import PolySimplex

class PolyFace(LinearCombination):
    def __init__(self, data=None):
        if data is None:
            super().__init__(term_type=Face)
            return
        prepared = {}
        if isinstance(data, dict):
            new_data = {}
            for key, coeff in data.items():
                if coeff == 0:
                    continue
                if not isinstance(key, Face):
                    # Преобразуем key в Face:
                    # если key — кортеж элементов (например, (a,b,c)) -> Face из одной компоненты
                    # если key — кортеж кортежей (например, ((a,b),(c,d))) -> Face из нескольких компонент
                    if isinstance(key, tuple) and all(isinstance(x, tuple) for x in key):
                        face = Face(list(key))
                    elif isinstance(key, tuple):
                        face = Face([key])
                    else:
                        raise TypeError(f"Cannot convert {key} to Face")
                    new_data[face] = new_data.get(face, 0) + coeff
                else:
                    new_data[key] = new_data.get(key, 0) + coeff
            prepared = new_data
        elif isinstance(data, (list, tuple)):
            new_data = {}
            for item in data:
                if len(item) != 2:
                    raise ValueError("Each item must be (face_data, coeff)")
                key, coeff = item
                if coeff == 0:
                    continue
                if not isinstance(key, Face):
                    if isinstance(key, tuple) and all(isinstance(x, tuple) for x in key):
                        face = Face(list(key))
                    elif isinstance(key, tuple):
                        face = Face([key])
                    else:
                        raise TypeError(f"Cannot convert {key} to Face")
                    new_data[face] = new_data.get(face, 0) + coeff
                else:
                    new_data[key] = new_data.get(key, 0) + coeff
            prepared = new_data
        elif isinstance(data, Face):
            prepared = {data: 1}
        elif isinstance(data, PolyFace):
            prepared = data.terms
        else:
            raise TypeError(f"Unsupported type: {type(data)}")
        super().__init__(prepared, term_type=Face)

    def to_polysimplex(self): # Преобразует полигрань в полисимплекс, применяя to_polysimplex к каждой грани.
        result = PolySimplex()
        for face, coeff in self.terms.items():
            ps = face.to_polysimplex()
            if not ps.is_zero():
                result = result + ps * coeff
        return result

    # Методы __str__, __repr__ можно унаследовать от LinearCombination, но для удобства переопределим
    def __str__(self):
        if self.is_zero():
            return "0"
        parts = []
        for face, coeff in self.terms.items():
            face_str = str(face)
            if coeff == 1:
                parts.append(face_str)
            elif coeff == -1:
                parts.append(f"-{face_str}")
            else:
                parts.append(f"{coeff}*{face_str}")
        result = parts[0]
        for part in parts[1:]:
            if part[0] == '-':
                result += f" - {part[1:]}"
            else:
                result += f" + {part}"
        return result

    def __repr__(self):
        if self.is_zero():
            return "PolyFace()"
        return f"PolyFace({self.__str__()})"