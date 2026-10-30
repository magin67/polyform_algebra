"""Ordered bases and links between coordinate representations."""

import numpy as np

from ..objects.simplex import Simplex
from .space import Space

class Basis:
    def __init__(self, elements, is_orthonormal=False, *, space=None):
        self.elements = list(elements)
        if len(set(self.elements)) != len(self.elements):
            raise ValueError("Basis elements must be distinct")
        self._index = {elem: i for i, elem in enumerate(self.elements)}
        self.is_orthonormal = is_orthonormal  # Legacy declaration, not a metric check.
        self.transitions = {}
        self._space = None
        if space is not None:
            self.bind_space(space)

    @property
    def space(self):
        """The owning space, or None for a legacy basis with unspecified context."""
        return self._space

    def bind_space(self, space):
        """Attach an unspecified basis once; never silently identify two spaces."""
        if not isinstance(space, Space):
            raise TypeError("space must be a Space")
        if len(self) != space.dimension:
            raise ValueError("basis size must equal space dimension")
        if self._space is not None and self._space is not space:
            raise ValueError("basis already belongs to a different space")
        self._space = space

    def reordered(self, elements):
        """Return a basis in the same space with an explicit full permutation."""
        elements = list(elements)
        if len(elements) != len(self) or set(elements) != set(self.elements):
            raise ValueError("reordering requires a full permutation of basis elements")
        return Basis(elements, self.is_orthonormal, space=self.space)

    def centroid(self, name=None):
        """Return (mean element, frame) for a nonempty grade-one Element basis.

        The frame has coefficients 1/n and records this Basis. By linearity,
        its affine weight is sum(epsilon(h_i))/n, not necessarily one in a
        mixed or vector basis. This is an arithmetic mean, with no metric or
        isotropic-locality assertion. The mean is not inserted into the basis.
        See PMG terminology, 'Centroid', and Spaces/Coordinates,
        'Bi-coordinates and frames'.
        """
        from ..objects.element import Element
        from ..combinations.polysimplex import Polysimplex

        if not self.elements:
            raise ValueError("centroid requires a nonempty basis")
        if not all(isinstance(elem, Element) for elem in self.elements):
            raise TypeError("centroid currently requires grade-one Element objects")
        if len(set(self.elements)) != len(self.elements):
            raise ValueError("centroid requires distinct basis elements")
        coefficient = 1.0 / len(self.elements)
        frame = Polysimplex([(elem, coefficient) for elem in self.elements])
        centroid = Element(name=name, frame=frame, affine_weight=frame.affine_weight())
        centroid.basis_for_frame[0] = self
        return centroid, frame

    def add_element(self, elem):
        """Append a new basis element and invalidate all incident transitions."""
        if elem not in self._index:
            if self.space is not None:
                raise ValueError("cannot resize a basis of a fixed-dimensional space")
            self._index[elem] = len(self.elements)
            self.elements.append(elem)
            self._invalidate_transitions()

    def remove_element(self, elem):
        """Remove an element and invalidate links in both directions; absent is a no-op."""
        if elem in self._index:
            if self.space is not None:
                raise ValueError("cannot resize a basis of a fixed-dimensional space")
            idx = self._index.pop(elem)
            del self.elements[idx]
            # Rebuild indices after removal.
            for i, e in enumerate(self.elements):
                self._index[e] = i
            self._invalidate_transitions()

    def _invalidate_transitions(self):
        """Remove both directions of every link affected by a basis change."""
        for other in list(self.transitions):
            other.transitions.pop(self, None)
        self.transitions.clear()

    def create_basis_from_frames(self, index=0):
        """
        Создаёт новый базис, состоящий из элементов, входящих во фреймы (представления)
        элементов текущего базиса.
        index: индекс фрейма (обычно 0)
        """
        new_elements = set()
        for elem in self.elements:
            # Проверяем, есть ли у элемента фрейм с указанным индексом
            if hasattr(elem, '__getitem__') and len(elem.frames) > index:
                frame = elem[index]   # Polysimplex
                # Извлекаем все элементы (точки, векторы) из фрейма
                new_elements.update(frame.get_elements('all'))
            else:
                # Если фрейма нет, добавляем сам элемент (может быть точкой или вектором)
                new_elements.add(elem)
        return Basis(list(new_elements))

    def __len__(self):
        return len(self.elements)

    def __getitem__(self, idx):
        return self.elements[idx]

    def __contains__(self, elem):
        return elem in self._index

    def __str__(self):
        return f"Basis({[str(elem) for elem in self.elements]})"

    def __repr__(self):
        return f"Basis({[repr(elem) for elem in self.elements]})"

    def index(self, elem):
        return self._index.get(elem, -1)

    def link_to(self, other_basis, frame_ind=0, tolerance=1e-10):
        """Link coordinate bases; declared spaces must be the same object."""
        from ..combinations.polysimplex import Polysimplex
        if self.space is not other_basis.space:
            raise ValueError("basis links require the same space; bind legacy bases first")
        n = len(self)
        m = len(other_basis)
        if n != m:
            raise ValueError("Размеры базисов должны совпадать")

        # Columns express old basis elements in the target basis: x_new=M x_old.
        M = np.zeros((n, n))
        for i, elem in enumerate(self.elements):
            try:
                frame = elem[frame_ind]  # A representation in other_basis must exist.
            except IndexError:
                raise ValueError(f"Элемент {elem} не имеет фрейма с индексом {frame_ind}")
          
            coeffs = np.zeros(n)
            for term, coeff in frame.terms.items():
                if len(term.components) != 1 or isinstance(term.components[0], tuple):
                    raise TypeError("Представление должно быть линейной комбинацией атомарных элементов")
                pt = term.components[0]
                j = other_basis.index(pt)
                if j == -1:
                    raise ValueError(f"Элемент {pt} не принадлежит other_basis")
                coeffs[j] += coeff
            M[:, i] = coeffs

        if np.linalg.matrix_rank(M, tol=tolerance) < n:
            raise ValueError("Матрица перехода вырождена")
        M_inv = np.linalg.inv(M)

        # Store both directions so invalidation remains symmetric.
        self.transitions[other_basis] = (M, M_inv)
        other_basis.transitions[self] = (M_inv, M)

        # Columns of the inverse express target elements in the source basis.
        for j, b_elem in enumerate(other_basis.elements):
            coeffs = M_inv[:, j]
            terms = {}
            for i, c in enumerate(coeffs):
                if abs(c) > tolerance:
                    elem = self[i]
                    terms[Simplex([elem])] = c
            frame = Polysimplex(terms, term_type=Simplex)
            b_elem.add_frame(frame, basis=self)

