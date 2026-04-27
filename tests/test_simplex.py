import pytest
import sympy as sp

from src import Simplex, Face, quForm, PolySimplex, AffineVector, PolyFace, Polyform

# ---------- Переменные ----------
a, b, c, d, e, f, g, h, r, x, y, z = sp.symbols('a b c d e f g h r x y z')

# ---------- Simplex ----------

def test_monoid_grade():
    s = Simplex([a,b])
    assert s.grade() == 2   # размер
    f = Face([(a,b)])
    assert f.grade() == 1   # len(comp)-1
    q = quForm([{a,b}])
    assert q.grade() == 1

def test_simplex_creation():
    s = Simplex([a, b])
    assert s.size() == 2
    assert s.grade() == 2
    assert not s.is_zero()
    assert s.components == (a, b)

def test_simplex_mul():
    s1 = Simplex([a, b])
    s2 = Simplex([c, d])
    s3 = s1 * s2
    assert s3.components == (a, b, c, d)
    # общий элемент -> нулевой симплекс
    s4 = Simplex([a, b]) * Simplex([b, c])
    assert s4.is_zero()

def test_simplex_boundary():
    s = Simplex([a, b])
    bnd = s.boundary()
    assert isinstance(bnd, PolySimplex)
    # граница [a,b] = [b] - [a]
    assert bnd.terms == {Simplex([b]): 1, Simplex([a]): -1}
    s2 = Simplex([a, b, c])
    bnd2 = s2.boundary()
    # [b,c] - [a,c] + [a,b]
    expected = {Simplex([b,c]): 1, Simplex([a,c]): -1, Simplex([a,b]): 1}
    assert bnd2.terms == expected

def test_simplex_permute():
    s = Simplex([a, b, c])
    s_perm = s.permute([1, 0, 2])  # (b, a, c) - нечётная перестановка
    # Ожидаем, что s_perm = -s (после приведения к каноническому виду)
    assert s_perm == -s

# ---------- PolySimplex ----------

def test_polysimplex_eq():
    s1 = PolySimplex({(a, b, c): 1, (d, e): 1})
    s2 = PolySimplex({(b, a, c): -1, (d, e): 1})
    assert s1 == s2

def test_polysimplex_add():
    s1 = PolySimplex({(a, b, c): 1, (d, e): 1})
    s2 = PolySimplex({(b, a, c): 1, (e, d): 1})
    assert (s1 + s2).is_zero()

def test_polysimplex_sub():
    s1 = PolySimplex({(a, b, c): 1, (d, e): 1})
    s2 = PolySimplex({(b, a, c): -1, (e, d): -1})
    assert (s1 - s2).is_zero()

def test_polysimplex_mul():
    s1 = PolySimplex({(a, b): 1})
    s2 = PolySimplex({(c, d): 1})
    prod = s1 * s2
    assert len(prod.terms) == 1
    assert list(prod.terms.keys())[0].components == (a, b, c, d)

def test_polysimplex_boundary():
    s = PolySimplex({Simplex([a, b ,c]): 2})
    bnd = s.boundary()
    assert bnd.terms[Simplex([b, c])] == 2
    assert bnd.terms[Simplex([a, c])] == -2
    assert bnd.terms[Simplex([a, b])] == 2

# ---------- AffineVector ----------

def test_affine_vector_sum_zero():
    with pytest.raises(ValueError):
        AffineVector({a: 1, b: 2})   # сумма не ноль

def test_affine_vector_from_edges():
    v = AffineVector({(a, b): 1, (c, d): 2})
    # b - a + 2(d - c) = -a + b -2c + 2d
    expected_terms = {Simplex([a]): -1, Simplex([b]): 1, Simplex([c]): -2, Simplex([d]): 2}
    assert v.terms == expected_terms

def test_affine_vector_from_vertices():
    v = AffineVector({a: 1, b: -1})
    assert v.terms == {Simplex([a]): 1, Simplex([b]): -1}

def test_affine_vector_dot():
    v = AffineVector({(a, b): 1})   # b - a
    w = AffineVector({(b, c): 1})   # c - b
    assert v.dot(w) == -1   # (-1)*0 + 1*1 + 0*(-1) = 1? Проверим: коэфф при a: -1, b:1, c:0; у w: a:0, b:1, c:-1 => -1*0 + 1*1 + 0*(-1)=1. На самом деле -1? Пересчитаем: v = -a + b, w = -b + c. Скалярное произведение = (-1)*0 + 1*(-1) + 0*1 = -1. Да, -1.
    assert v.dot(v) == 2   # 1^2 + (-1)^2 = 2

def test_affine_vector_bilinear_form():
    v = AffineVector({(a,b): 1})
    w = AffineVector({(b,c): 1})
    qf = v @ w
    # Ожидаем три формы: {a,b}, {a,c}, {b,c} с коэффициентами -0.5, 0.5, -0.5
    assert len(qf.terms) == 3
    forms = {frozenset(k.components[0]) for k in qf.terms}
    expected = {frozenset({a,b}), frozenset({a,c}), frozenset({b,c})}
    assert forms == expected
    # Сумма коэффициентов должна быть -0.5 (проверка необязательна)

def test_affine_vector_to_polyform():
    v = AffineVector({(a, b): 1})
    qf = v.to_polyform()
    assert len(qf.terms) == 1
    form = list(qf.terms.keys())[0]
    # Сравниваем компоненты как множества, так как порядок не важен
    assert set(form.components[0]) == {a, b}
    assert qf.terms[form] == 1

# ---------- Face ----------

def test_face_creation():
    f = Face([(a,b), (c,d)])
    assert f.cardinality() == 2
    assert f.size() == 4
    assert f.grade() == 2
    assert not f.is_zero()

def test_face_mul_single_overlap():
    f1 = Face([(a,b)])
    f2 = Face([(b,c)])
    prod = f1 * f2
    # (a,b)*(b,c) -> (a,b,c)
    assert len(prod.components) == 1
    assert tuple(prod.components[0]) == (a,b,c)

def test_face_mul_no_overlap():
    f1 = Face([(a,b)])
    f2 = Face([(c,d)])
    prod = f1 * f2
    assert prod.cardinality() == 2
    assert prod.components[0] == (a,b)
    assert prod.components[1] == (c,d)

def test_face_mul_two_common():
    f1 = Face([(a,b,c)])
    f2 = Face([(a,b)])
    prod = f1 * f2
    # два общих элемента -> нулевая грань
    assert prod.is_zero()

# ---------- PolyFace ----------

def test_polyface_from_dict():
    pf = PolyFace({(a,b,c): 1, (b,c,d): 2})
    assert len(pf.terms) == 2
    for term in pf.terms:
        assert isinstance(term, Face)

def test_polyface_from_list():
    pf = PolyFace([((a,b,c), 1), ((b,c,d), 2)])
    assert len(pf.terms) == 2

def test_polyface_to_polysimplex():
    face = Face([(a,b,c), (d,e)])
    pf = PolyFace(face)
    ps = pf.to_polysimplex()
    assert len(ps.terms) == 6   # (bnd треугольника) * (bnd ребра) даёт 3*2=6 слагаемых

# ---------- Polyform ----------

def test_polyform_from_components():
    v = AffineVector({(a,b): 1})          # b - a
    form = quForm([{b,c}])                # форма из ребра b-c
    qf = Polyform.from_components([v, form])
    assert not qf.is_zero()
    # Проверим, что результат содержит форму {a,b,c} с ненулевым коэффициентом
    expected_form = quForm([{a,b,c}])
    assert qf.terms.get(expected_form, 0) != 0