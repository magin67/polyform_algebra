
from .simplex import Simplex

class Point(Simplex):
    def __init__(self, name=None, obj=None):
        self.obj = obj
        if name is not None: self.name = str(name)
        elif obj is not None: self.name = str(obj)
        else: self.name = None
        super().__init__([self])   # компонента – ссылка на себя
        self.multiplicity = 1

    def __hash__(self):
        return hash((self.name, id(self.obj)))

    def __eq__(self, other):
        return isinstance(other, Point) and self.name == other.name and self.obj is other.obj

    def __repr__(self):
        return f"Point({self.name})"

    def __str__(self):
        return self.name if self.name is not None else f"p{id(self)}"
    
    @staticmethod
    def create_list(names_or_objs):
        if isinstance(names_or_objs, (list, tuple)):
            return [Point(name) for name in names_or_objs]
        else:
            return [Point(name) for name in names_or_objs]  # на самом деле нужно обработать случай одного элемента
