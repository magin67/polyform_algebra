# @title Полиформа границ - косвенное задание

from .._core.basis import Basis
from .affine_vector import AffineVector
from ..objects.quform import quForm
from ..objects.bform import bForm

from .polyform import Polyform

class PolyBForm(Polyform):
    @classmethod
    def from_eigenvectors(cls, poly, normalized=False, tolerance=1e-10):
        eigen_pairs = poly.eigenvectors(normalized=normalized, tolerance=tolerance)
        terms = {}
        basis_vectors = []
        for ev, vec in eigen_pairs:
            if abs(ev) < tolerance:
                continue
            basis_vectors.append(vec)        # vec — это Vector (или AffineVector)
            bf = bForm([vec])
            if normalized:
                terms[bf] = terms.get(bf, 0) + ev
            else:
                terms[bf] = terms.get(bf, 0) + 1
        obj = cls(terms)
        obj._basis = Basis(basis_vectors)
        obj._basis.is_orthogonal = True
        return obj

    @classmethod
    def from_vectors(cls, vectors):
        return cls(bForm(vectors))

    def __init__(self, data=None):
        super().__init__(data, term_type=bForm)

    def _term_str(self, term, use_repr=False):
        if term.is_one(): return "[1]"
        vectors_str = ", ".join(str(v) for v in term.components)
        return f"[{vectors_str}]"

    def _create_empty(self):
        return PolyBForm()

    def to_laplacian(self):
        return self.expand().to_laplacian()

    def eigenvectors(self, tolerance=1e-10):
        return self.expand().eigenvectors(tolerance)

    def expand(self): # Преобразует PolyBForm в Polyform, раскрывая каждый bForm
        result = Polyform()  # пустая полиформа (0)
        for bf, coeff in self.terms.items():
            # Произведение полиформ всех векторов в bf
            prod = Polyform(quForm.one())  # единица
            for v in bf.components:
                prod = prod * v.to_polyform()
            # Умножаем на коэффициент
            if coeff != 1:
                prod = prod * coeff
            result = result + prod
        return result

    # Метрика
    def exp(self): # Копипаст от полиформы, но пока так
        if not self.is_homogeneous() or self.grade() != 1: return super().exp() # Если не диагональная, используем стандартный (медленный) метод из родителя
        self._skip_basis_update = True
        result = PolyBForm({bForm.one(): 1})   # единица
        for frm, coeff in self.terms.items():
            # Создаём терм (1 + coeff * frm) как полиформу
            term_terms = {bForm.one(): 1}
            if coeff != 0:
                term_terms[frm] = term_terms.get(frm, 0) + coeff
            term = PolyBForm(term_terms)
            result = result * term
        self._skip_basis_update = False
        return result
