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

# ---------- Polyform ----------

def test_quForm_grade():
    qf = quForm([(a, b)])
    assert qf.grade() == 1




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