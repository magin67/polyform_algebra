import pytest
import sys
import os

# Добавляем корень проекта в путь
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.objects.element import Point, Vector
from src.objects.simplex import Simplex
# from combinations.polysimplex import Polysimplex
# from objects.quform import quForm
# from combinations.polyform import Polyform

# ---------- Точки ----------
a, b, c, d, e, f, g, h, r, x, y, z = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'r', 'x', 'y', 'z'])

# ---------- Point ----------

def test_point_mul():
    smp = a * b
    assert smp == Simplex([a, b])

# ---------- Vector ----------

def test_vector_comps():
    v = Vector('v', frame=b - a - 2*c + 2*d)
    # b - a + 2(d - c) = -a + b -2c + 2d
    expected_terms = {Simplex([a]): -1, Simplex([b]): 1, Simplex([c]): -2, Simplex([d]): 2}
    assert v[0].terms == expected_terms

def test_vector_from_vertices():
    v = Vector('v', frame=(a - b))
    assert v[0].terms == {Simplex([a]): 1, Simplex([b]): -1}

def test_vector_lists():
    v = Vector('v', frame=(a + b - 2*c))
    els, coeffs = v[0].terms_lists()
    assert set(els) == set([a, b, c])
    assert set(coeffs) == set([1, 1, -2])

def test_vector_bilinear_form():
    v = Vector('v', frame=(b - a))
    w = Vector('w', frame=(c - b))
    qf = v[0] @ w[0]
    # Ожидаем три формы: {a,b}, {a,c}, {b,c} с коэффициентами -0.5, 0.5, -0.5
    assert len(qf.terms) == 3
    forms = {frozenset(k.components[0]) for k in qf.terms}
    expected = {frozenset({a, b}), frozenset({a, c}), frozenset({b, c})}
    assert forms == expected
    # Сумма коэффициентов должна быть -0.5 (проверка необязательна)

def test_vector_to_polyform():
    v = Vector('v', frame=(b - a))
    qf = v[0] @ v[0]
    assert len(qf.terms) == 1
    form = list(qf.terms.keys())[0]
    # Сравниваем компоненты как множества, так как порядок не важен
    assert set(form.components[0]) == {a, b}
    assert qf.terms[form] == 1


'''

# ---------- Polyform ----------

def test_quForm_grade():
    q = quForm([(a, b)])
    assert q.grade() == 1


def test_polyform_from_components():
    v = Vector({(a,b): 1})          # b - a
    form = quForm([{b,c}])                # форма из ребра b-c
    qf = Polyform.from_components([v, form])
    assert not qf.is_zero()
    # Проверим, что результат содержит форму {a,b,c} с ненулевым коэффициентом
    expected_form = quForm([{a,b,c}])
    assert qf.terms.get(expected_form, 0) != 0

'''