# @title Базис

import numpy as np

from ..objects.simplex import Simplex

class Basis:
    def __init__(self, elements, is_orthonormal=False):
        self.elements = list(elements)          # список базисных элементов (Point, Vector)
        self._index = {elem: i for i, elem in enumerate(self.elements)}
        self.is_orthonormal = is_orthonormal      # флаг, можно вычислить позже
        self.transitions = {}  # словарь для хранения матриц перехода

    def add_element(self, elem):
        if elem not in self._index:
            self._index[elem] = len(self.elements)
            self.elements.append(elem)
            # При добавлении элемента все существующие связи становятся недействительными
            self.transitions.clear()

    def remove_element(self, elem):
        if elem in self._index:
            idx = self._index.pop(elem)
            del self.elements[idx]
            # Обновляем индексы
            for i, e in enumerate(self.elements):
                self._index[e] = i
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

    def link_to(self, other_basis, frame_ind=0, tolerance=1e-10): # Связывание базисов
        from ..combinations.polysimplex import Polysimplex
        n = len(self)
        m = len(other_basis)
        if n != m:
            raise ValueError("Размеры базисов должны совпадать")

        # Строим матрицу перехода из self в other, используя существующие представления
        M = np.zeros((n, n))
        for i, elem in enumerate(self.elements):
            try:
                frame = elem[frame_ind]   # представление в other_basis (должно существовать)
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

        # Сохраняем матрицы
        self.transitions[other_basis] = (M, M_inv)
        other_basis.transitions[self] = (M_inv, M)

        # Добавляем представления элементам other_basis в базисе self (обратные)
        for j, b_elem in enumerate(other_basis.elements):
            coeffs = M_inv[:, j]   # столбец j обратной матрицы
            terms = {}
            for i, c in enumerate(coeffs):
                if abs(c) > tolerance:
                    elem = self[i]
                    terms[Simplex([elem])] = c
            frame = Polysimplex(terms, term_type=Simplex)
            b_elem.add_frame(frame, basis=self)

