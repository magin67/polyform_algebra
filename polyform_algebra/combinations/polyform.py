# @title Полиформа - линейная комбинация форм

from collections import Counter
from numbers import Real
import numpy as np

from .polysimplex import Polysimplex
from ..objects.simplex import Simplex
from ..objects.quform import quForm
from ..objects.element import Element, Point   # новый класс вектора

class Polyform(Polysimplex):
    @classmethod
    def zero(cls):
        return cls()
    
    @classmethod
    def one(cls):
        return cls({quForm.one(): 1})

    @classmethod
    def from_networkx(cls, G, return_points=None):
        """Import a loop-free undirected network with finite non-negative weights.

        Missing weights mean one; zero weights are absent connections. Parallel
        conductances are summed. Comparable node labels retain sorted order;
        otherwise insertion order is used. The returned mapping/list identifies
        points even when different labels have the same string representation.
        Signed algebraic expressions should be built directly, not imported here.
        """
        import networkx as nx

        if not isinstance(G, nx.Graph):
            raise TypeError("Expected a NetworkX graph")
        if G.is_directed():
            raise ValueError("Directed graphs are not supported by this network constructor")
        if return_points not in (None, 'dict', 'list'):
            raise ValueError("return_points must be None, 'dict' or 'list'")
        nodes = list(G.nodes())
        try:
            nodes = sorted(nodes)
        except TypeError:
            pass  # Mixed label types need not define a common ordering.

        # Element equality uses names: reserve all original names before adding
        # suffixes so that 1, '1' and '1#1' cannot accidentally merge vertices.
        names = [str(node) for node in nodes]
        reserved = set(names)
        used = set()
        points = {}
        for node, name in zip(nodes, names):
            candidate = name
            suffix = 1
            while not candidate or candidate in used:
                candidate = f"{name}#{suffix}"
                suffix += 1
                if candidate in reserved:
                    candidate = ''
            used.add(candidate)
            points[node] = Point(candidate)

        terms = Counter()
        for u, v, data in G.edges(data=True):
            if u == v:
                raise ValueError("Self-loops are not supported by this network constructor")
            weight = data.get('weight', 1)
            if isinstance(weight, (bool, np.bool_)) or not isinstance(weight, Real):
                raise TypeError("Edge weights must be real numbers")
            try:
                weight = float(weight)
            except OverflowError as exc:
                raise ValueError("Edge weights must be finite floating-point values") from exc
            if not np.isfinite(weight) or weight < 0:
                raise ValueError("Edge weights must be finite and non-negative")
            if weight:
                form = quForm([(points[u], points[v])])
                terms[form] += weight
                if not np.isfinite(terms[form]):
                    raise ValueError("Sum of parallel conductances must be finite")
        poly = cls(terms, vertices=points.values())
        if return_points == 'dict':
            return poly, points
        if return_points == 'list':
            return poly, list(points.values())
        return poly

    @classmethod
    def from_laplacian(cls, laplacian_matrix, vertex_labels, tolerance=1e-10):
        """Import a real undirected network Laplacian in the supplied Point order.

        Require a finite square symmetric matrix, non-positive off-diagonal
        entries and zero row sums. All checks use an absolute tolerance (no
        relative tolerance). Small asymmetry is averaged and positive roundoff
        off-diagonal entries are clipped to zero. Conductances at or below
        tolerance are omitted. Empty and disconnected
        networks are valid. Labels must be distinct Point objects.
        """
        if (isinstance(tolerance, (bool, np.bool_)) or
                not isinstance(tolerance, Real) or
                not np.isfinite(tolerance) or tolerance < 0):
            raise ValueError("Tolerance must be finite and non-negative")
        try:
            L = np.asarray(laplacian_matrix)
        except (ValueError, TypeError) as exc:
            raise ValueError("Matrix must be a rectangular two-dimensional array") from exc
        if L.ndim != 2 or L.shape[0] != L.shape[1]:
            raise ValueError("Matrix must be two-dimensional and square")
        if L.dtype.kind not in 'iuf':
            raise TypeError("Matrix entries must be real numbers")
        L = L.astype(float)
        if not np.all(np.isfinite(L)):
            raise ValueError("Matrix entries must be finite")
        n = L.shape[0]
        labels = list(vertex_labels)
        if len(labels) != n:
            raise ValueError("Number of vertex labels must match matrix size")
        if any(not isinstance(label, Point) for label in labels):
            raise TypeError("Vertex labels must be Point objects")
        if len(set(labels)) != n or len({str(label) for label in labels}) != n:
            raise ValueError("Vertex labels must have distinct identities and names")
        if not np.allclose(L, L.T, atol=tolerance, rtol=0):
            raise ValueError("Matrix must be symmetric")
        off_diagonal = L[~np.eye(n, dtype=bool)]
        if np.any(off_diagonal > tolerance) or np.any(np.diag(L) < -tolerance):
            raise ValueError("Laplacian conductances must be non-negative")
        with np.errstate(over='ignore', invalid='ignore'):
            row_sums = L.sum(axis=1)
        if not np.all(np.isfinite(row_sums)) or not np.all(np.abs(row_sums) <= tolerance):
            raise ValueError("Laplacian row sums must be zero")
        terms = Counter()
        for i in range(n):
            for j in range(i+1, n):
                weight = max(0.0, -(L[i, j] / 2 + L[j, i] / 2))
                if weight > tolerance:
                    terms[quForm([(labels[i], labels[j])])] += weight
        return cls(terms, vertices=labels)

    @classmethod
    def from_pairs(cls, pairs):
        """Build a polyform, adding coefficients of repeated quadratic terms."""
        return cls([(quForm(key), coeff) for key, coeff in pairs])

    def __init__(self, data=None, *, vertices=None):
        if vertices is None and isinstance(data, Polyform):
            vertices = data._vertices
        self._vertices = None if vertices is None else tuple(dict.fromkeys(vertices))
        super().__init__(data, term_type=quForm)
        if self._vertices is not None and not self.get_elements().issubset(self._vertices):
            raise ValueError("Form contains elements outside the supplied vertices")

    @property
    def vertices(self):
        """Фиксированные вершины графа, включая изолированные; иначе носитель записи."""
        if self._vertices is not None:
            return self._vertices
        return tuple(sorted(self.get_elements(), key=str))

    def _from_terms(self, terms):
        """Rebuild selected terms without losing isolated graph vertices."""
        return self.__class__(terms, vertices=self._vertices)

    def _with_context(self, result, other=None):
        """Переносит вершины, которые могут исчезнуть из ненулевых слагаемых.

        Создаём новый объект, чтобы арифметика не меняла контекст операндов.
        """
        if result is NotImplemented: return result
        explicit = self._vertices is not None or getattr(other, '_vertices', None) is not None
        if not explicit: return result
        vertices = list(self.vertices)
        if isinstance(other, Polyform):
            vertices.extend(other.vertices)
        elif isinstance(other, Simplex):
            vertices.extend(sorted(other.get_elements(), key=str))
        return Polyform(result, vertices=vertices)

    def __add__(self, other):
        return self._with_context(super().__add__(other), other)

    def __neg__(self):
        return self._with_context(super().__neg__())

    def __mul__(self, other):
        return self._with_context(super().__mul__(other), other)

    def __rmul__(self, scalar):
        return self._with_context(super().__rmul__(scalar))

    def __pow__(self, exponent):
        return self._with_context(super().__pow__(exponent))

    def extract_grade(self, grade):
        return self._with_context(super().extract_grade(grade))

    def _spanning_form(self):
        """Выбирает T_n по всем вершинам, а не по текущему старшему грейду.

        У несвязного графа старший грейд ниже n−1; его коэффициент не равен τ.
        """
        # Старое базисное представление с атомарными формами имеет другой объём.
        if any(not term.is_one() and term.kind() != 'boundary' for term in self.terms):
            return None
        vertices = self.vertices
        if not vertices:
            return None
        return quForm([tuple(vertices)]) if len(vertices) > 1 else quForm.one()

    def _validate_query(self, obj):
        if not isinstance(obj, (quForm, Polyform)):
            raise TypeError("Expected quForm or Polyform")
        if not obj.get_elements().issubset(self.vertices):
            raise ValueError("Query contains elements outside the metric vertices")

    @staticmethod
    def add_term(prepared, term, coeff):
        cterm = term.canonical()
        key = quForm(cterm.components, sign=1)
        prepared[key] = prepared.get(key, 0) + coeff * cterm.sign

    def __str__(self):
        from .._core.lincomb import LinearCombination
        return LinearCombination.__str__(self)

    def __repr__(self):
        return f"Polyform({self.__str__()})"

    def is_equivalent(self, other, tolerance=1e-10):
        """Compare quadratic forms after expanding boundaries into atomic blades.

        For X = sum(x_I e_I), its square has coefficients x_I*x_J.
        These internal coefficient pairs are only a comparison device: they do
        not add general bilinear objects to the public quadratic algebra.
        Unlike multiplication probes, this also distinguishes top-grade forms
        and scalars from zero. Element frames are not expanded implicitly.
        """
        if tolerance < 0 or not np.isfinite(tolerance):
            raise ValueError("Tolerance must be finite and non-negative")
        if not isinstance(other, Polyform):
            return False
        coefficients = Counter()
        for form, direction in ((self, 1), (other, -1)):
            for term, weight in form.terms.items():
                # Retain orientation until the outer square is formed.
                expansion = Simplex(list(term.components)).to_polysimplex()
                items = list(expansion.terms.items())
                for i, (left, a) in enumerate(items):
                    for j in range(i, len(items)):
                        right, b = items[j]
                        key = frozenset((left, right))
                        factor = 1 if i == j else 2
                        coefficients[key] += direction * weight * factor * a * b
        return all(abs(value) <= tolerance for value in coefficients.values())

    # --- Методы метрики ---
    def potential(self, obj_polyform):
        """Коэффициент произведения при фиксированной остовной форме T_n."""
        self._validate_query(obj_polyform)
        spanning = self._spanning_form()
        prod = self * obj_polyform
        if spanning is not None:
            return prod.terms.get(spanning, 0)
        # Совместимость с существующими атомарными базисами: фиксируем грейд self.
        top = prod.extract_grade(self.grade())
        lead = top.leading_coefficient()
        if isinstance(lead, list):
            raise ValueError("Metric needs a fixed top form; multiple coefficients found")
        return lead

    def norm(self, obj_polyform):
        """Нормированный потенциал; нулевая остовная компонента недопустима."""
        self._validate_query(obj_polyform)
        spanning = self._spanning_form()
        metric_lead = self.terms.get(spanning, 0) if spanning is not None else self.leading_coefficient()
        if isinstance(metric_lead, list):
            raise ValueError("Metric needs a fixed top form; multiple coefficients found")
        if metric_lead == 0:
            raise ValueError("Metric has zero spanning coefficient; normalization is undefined")
        return self.potential(obj_polyform) / metric_lead

    @staticmethod
    def polar_form(p, q, tolerance=1e-10):
        # Это полярная форма элементов (векторов и/или точек) p, q, для которых известны фреймы (Polysimplex)
        # Собираем базисные элементы
        if isinstance(p, Simplex): p = p.as_polysimplex()
        if isinstance(q, Simplex): q = q.as_polysimplex()

        for obj in (p, q):
            if not isinstance(obj, Polysimplex) or isinstance(obj, Polyform):
                raise TypeError("Expected an external object, not a form")
            if not obj.is_zero() and (not obj.is_homogeneous() or obj.grade() != 1):
                raise ValueError("Polarization currently supports only homogeneous grade 1 objects")
        # Раскрываем (ab) в b-a перед извлечением атомарных коэффициентов.
        def expand(obj):
            result = Polysimplex.zero()
            for term, coeff in obj.terms.items():
                result += coeff * term.to_polysimplex()
            return result
        p, q = expand(p), expand(q)

        basis_elements = sorted(p.get_elements('all') | q.get_elements('all'), key=str)

        coeff_p = p.as_coefficient_dict(basis_elements)
        coeff_q = q.as_coefficient_dict(basis_elements)
        
        sum_p = sum(coeff_p.values())
        sum_q = sum(coeff_q.values())
        
        # Диагональная часть
        diag_terms = {}
        for elem in basis_elements:
            c = sum_p * coeff_q[elem] + sum_q * coeff_p[elem]
            if abs(c) > tolerance:
                form = quForm([elem])
                diag_terms[form] = diag_terms.get(form, 0.0) + c
        
        # Векторная часть
        diff_terms = {}
        n = len(basis_elements)
        for i in range(n):
            for j in range(i+1, n):
                b_i = basis_elements[i]
                b_j = basis_elements[j]
                c = coeff_p[b_i] * coeff_q[b_j] + coeff_q[b_i] * coeff_p[b_j]
                if abs(c) > tolerance:
                    form = quForm([(b_i, b_j)])
                    diff_terms[form] = diff_terms.get(form, 0.0) + c
        
        # Полярная форма = диагональная - разностная
        terms = diag_terms.copy()
        for form, c in diff_terms.items():
            terms[form] = terms.get(form, 0.0) - c
        
        # Округление
        cleaned = {form: Polyform._round_value(c, tolerance) for form, c in terms.items() if abs(c) > tolerance}
        return Polyform(cleaned)

    @staticmethod
    def scalar_product(p, q, tolerance=1e-10):
        # p, q - фреймы (Polysimplex)
        return Polyform.polar_form(p, q, tolerance) * 0.5

    def to_matrix(self):
        """
        Преобразует полиформу 1-го порядка в матрицу.
        Для термов quForm([(a,b)]) добавляет внедиагональные элементы,
        для термов quForm([a]) — диагональные.
        Возвращает (matrix, index), где index — словарь {элемент: индекс}.
        """
        if not self.is_zero() and (not self.is_homogeneous() or self.grade() != 1):
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")

        # Собираем все элементы (точки или векторы) из полиформы
        elements = self.vertices
        n = len(elements)
        idx = {elem: i for i, elem in enumerate(elements)}
        M = np.zeros((n, n))

        for term, coeff in self.terms.items():
            if len(term.components) != 1:
                raise ValueError("Форма должна состоять из одной компоненты")
            comp = term.components[0]
            if isinstance(comp, tuple) and len(comp) == 2: # Граничная форма (b - a)
                a, b = comp
                if a not in idx or b not in idx: continue
                i, j = idx[a], idx[b]
                M[i, j] -= coeff
                M[j, i] -= coeff
                M[i, i] += coeff
                M[j, j] += coeff
            else: # Атомарная форма (точка или вектор)
                if comp not in idx: continue
                i = idx[comp]
                M[i, i] += coeff

        return M, idx

    # --- Диагонализация - приведение к базису из собственных векторов ---

    def to_eigenbasis(self, tolerance=1e-10):
        """
        Преобразует однородную полиформу 1-го порядка в диагональную форму в базисе собственных векторов.
        Возвращает Polyform: Σ λ_i * quForm([v_i]), где v_i — собственный вектор (как объект Vector),
        λ_i — собственное число.
        """
        # import numpy as np
        # from objects.element import Vector
        # from .polysimplex import Polysimplex
        from .._core.basis import Basis

        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Полиформа должна быть однородной и иметь грейд 1")

        # Получаем матрицу (лапласиан) и отображение элементов->индексы
        M, idx = self.to_matrix()   # M — матрица, idx — {element: index}
        n = M.shape[0]
        if n == 0: return Polyform.zero(), Basis([])

        # Вычисляем собственные значения и векторы (отсортированы по возрастанию)
        eig_vals, eig_vecs = np.linalg.eigh(M)

        # Словарь для новой полиформы и список для базиса
        new_terms = {}
        basis_elements = []

        for k in range(n):
            ev = eig_vals[k]
            vec = eig_vecs[:, k]

            # Нормируем вектор (если уже нормирован, то норма = 1, но для нулевого собственного числа он может быть не нормирован)
            norm = np.linalg.norm(vec)
            if norm > tolerance: vec = vec / norm
            else: continue # Если вектор нулевой (practically), пропускаем? Не должно быть.

            # Строим фрейм (Polysimplex) для собственного вектора в исходном базисе
            frame_terms = {}
            for elem, i in idx.items():
                coeff = vec[i]
                if abs(coeff) < tolerance: continue
                frame_terms[elem] = coeff
            frame = Polysimplex(frame_terms)

            # Создаём элемент с этим фреймом и добавляем в базис
            mult = frame.multiplicity()
            elem = Element(name='v' + str(k), frame=frame, multiplicity=mult)
            basis_elements.append(elem)

            # Добавляем в полиформу (если собственное число не нулевое), коэффициент = собственное число.
            if abs(ev) > tolerance:
                qform = quForm([elem])
                new_terms[qform] = new_terms.get(qform, 0.0) + ev

        # Создаём ортонормированный базис
        basis = Basis(basis_elements, is_orthonormal = True)
        # Возвращаем полиформу и базис
        return Polyform(new_terms), basis

    def expand(self, index=0):
        if not self.is_homogeneous() or self.grade() != 1:
            raise ValueError("Expand применим только к однородным полиформам 1-го порядка")
        result = Polyform.zero()
        for term, coeff in self.terms.items():
            # term — quForm([el])
            el = term.components[0]   # это Element
            if not isinstance(el, Element):
                raise ValueError("Терм должен содержать элемент")
            frame = el[index]   # Polysimplex
            bilinear = Polyform.scalar_product(frame, frame)  # получаем полиформу нормы элемента
            result += bilinear * coeff
        return result

    def exp(self, tolerance=1e-10):
        """Exponentiate commuting simple squares, including a scalar component.

        Each supported positive-grade term is decomposable, so its square
        vanishes and exp(c*f) = 1+c*f. The scalar unit is not nilpotent and
        contributes the ordinary real exponential instead.
        """
        from math import exp

        if tolerance < 0 or not np.isfinite(tolerance):
            raise ValueError("Tolerance must be finite and non-negative")
        scalar = self.terms.get(quForm.one(), 0)
        result = Polyform(Polyform.one(), vertices=self._vertices)
        for term, coeff in self.terms.items():
            if term.is_one():
                continue
            if not (term * term).is_zero():
                raise ValueError("Factorized exponential requires square-zero terms")
            result = result * (Polyform.one() + Polyform({term: coeff}))
        return exp(scalar) * result
