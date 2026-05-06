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
    def zero(cls):
        return cls([])
    
    def __mul__(self, other):
        if not isinstance(other, quForm):
            return NotImplemented
        result = super().__mul__(other)
        if result.is_zero():
            return quForm.zero()
        # Результат — симплекс. Приводим к quForm, устанавливая знак = 1
        return quForm(result.components, sign=1)
    
    def canonical(self):
        # Получаем каноническую форму от Simplex (с учётом знака)
        can = super().canonical()
        # Принудительно устанавливаем знак = 1
        can.sign = 1
        can.__class__ = quForm
        return can
    
    def __str__(self):
        if self.is_zero():
            return "0"
        if self.is_one():
            return "1"
        inner = ", ".join(str(c) for c in self.components)
        return f"f[{inner}]"    

    def __repr__(self):
        return f"quForm({self.components})"
