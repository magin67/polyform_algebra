import pytest
import sys
import os

from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.simplex import Simplex
from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.objects.quform import quForm
from polyform_algebra.combinations.polyform import Polyform

# ---------- Точки ----------
a, b, c, d, e, f, g, h, r, x, y, z = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'r', 'x', 'y', 'z'])

# ---------- И формы ----------
qf_ab = quForm([(a, b)])
qf_bc = quForm([(b, c)])
qf_cd = quForm([(c, d)])
qf_abc = quForm([(a, b, c)])
qf_bcd = quForm([(b, c, d)])
qf_ab_cd = quForm([(a, b), (c, d)])
qf_abcd = quForm([(a, b, c, d)])

# ---------- quForm ----------

def test_quForm_grade():
    assert qf_ab.grade() == 1

def test_quForm_mul_zero():
    assert qf_abc*qf_ab == quForm.zero()  # так как (a,b) входит в (a,b,c), результат должен быть нулевым

def test_quForm_mul_one():
    assert qf_ab*qf_bc == qf_abc

def test_quForm_mul_one2():
    assert qf_bc*qf_ab_cd == qf_abcd


# ---------- Polyform ----------

def test_Polyform_equal():
    L1 = Polyform.from_pairs([
        ([(a, b, c)], 1),
        ([(a, b, d)], 1),
        ([(a, c, d)], 1),
        ([(b, c, d)], 1)
    ])
    L2 = Polyform.from_pairs([
        ([(a, b), (c, d)], 1),
        ([(a, c), (b, d)], 1),
        ([(a, d), (b, c)], 1)
    ])
    assert L1.is_equivalent(L2) == True



'''
def test_polyform_from_components():
    v = Vector({(a,b): 1})          # b - a
    form = quForm([{b,c}])                # форма из ребра b-c
    qf = Polyform.from_components([v, form])
    assert not qf.is_zero()
    # Проверим, что результат содержит форму {a,b,c} с ненулевым коэффициентом
    expected_form = quForm([{a,b,c}])
    assert qf.terms.get(expected_form, 0) != 0

'''