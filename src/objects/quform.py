# @title Форма – это квадратичный симплекс

from .simplex import Simplex

class quForm(Simplex):
    """Квадратичная форма от симплекса (знак всегда положительный)."""
    
    def __init__(self, data, sign=1, dual=False, multiplicity=None):
        """
        data: может быть
          - списком компонент (атомов и/или кортежей)
          - кортежем (границей)
          - объектом Simplex
        """
        if isinstance(data, Simplex):
            # Берём компоненты из симплекса, игнорируем его знак
            super().__init__(data.components, sign=1, dual=data.dual, multiplicity=data.multiplicity)
        elif isinstance(data, tuple):
            # Один кортеж как граница
            super().__init__([data], sign=1, dual=dual, multiplicity=multiplicity)
        elif isinstance(data, list):
            super().__init__(data, sign=1, dual=dual, multiplicity=multiplicity)
        else:
            raise TypeError("quForm data must be list, tuple or Simplex")
    
    @classmethod
    def one(cls):
        obj = Simplex.one()
        obj.__class__ = cls
        return obj
    
    @classmethod
    def zero(cls): # создаём нулевой симплекс с multiplicity=-1
        return cls([], multiplicity=-1)

    def is_zero(self):
        return len(self.components) == 0 and self.sign == -1

    def is_one(self):
        return len(self.components) == 0 and self.sign == 1

    def __mul__(self, other):
        if not isinstance(other, quForm):
            return NotImplemented
        result = super().__mul__(other)
        if result.is_zero():
            return quForm.zero()
        return quForm(result)   # передаём симплекс, а не его components

    def canonical(self):
        # Сортируем элементы внутри каждого кортежа
        sorted_components = []
        for comp in self.components:
            if isinstance(comp, tuple):
                # Сортируем элементы кортежа по строковому представлению (или по id, но строки надёжнее)
                sorted_comp = tuple(sorted(comp, key=str))
                sorted_components.append(sorted_comp)
            else:
                sorted_components.append(comp)
        # Сортируем компоненты между собой
        sorted_components.sort(key=lambda x: str(x))
        # Создаём новый объект с sign=1 (игнорируем старый знак)
        return quForm(sorted_components, sign=1)
    
    def __str__(self):
        if self.is_zero(): return "<0>"
        if self.is_one(): return "<1>"
        inner = ", ".join(str(c) for c in self.components)
        return f"f[{inner}]"    

    def __repr__(self):
        return f"quForm({self.components})"
