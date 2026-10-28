import pytest
import sys
import os

from polyform_algebra.objects.element import Point, Vector
from polyform_algebra.objects.simplex import Simplex
from polyform_algebra.combinations.polysimplex import Polysimplex
from polyform_algebra.objects.quform import quForm
from polyform_algebra.combinations.polyform import Polyform

# ---------- Points ----------
a, b, c, d, e, f, g, h, r, x, y, z = Point.create_list(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'r', 'x', 'y', 'z'])

# ---------- Quadratic forms ----------
qf_ab = quForm([(a, b)])
qf_ac = quForm([(a, c)])
qf_ad = quForm([(a, d)])
qf_bc = quForm([(b, c)])
qf_bd = quForm([(b, d)])
qf_cd = quForm([(c, d)])

qf_abc = quForm([(a, b, c)])
qf_bcd = quForm([(b, c, d)])

qf_ab_cd = quForm([(a, b), (c, d)])
qf_abcd = quForm([(a, b, c, d)])

# ---------- quForm ----------

def test_quForm_grade():
    """The square of an edge boundary has grade one, matching its argument."""
    assert qf_ab.grade() == 1

def test_quForm_mul_zero():
    """An edge contained in a triangle boundary makes their exterior product zero."""
    assert qf_abc*qf_ab == quForm.zero()  # The shared edge forces the product to vanish.

def test_quForm_mul_one():
    """Edges sharing one vertex merge: [(ab)]^2*[(bc)]^2=[(abc)]^2."""
    assert qf_ab*qf_bc == qf_abc

def test_quForm_mul_one2():
    """Joining two disjoint edges by bc produces the four-vertex boundary square."""
    assert qf_bc*qf_ab_cd == qf_abcd


# ---------- Polyform ----------

def test_Polyform_equal():
    """The sum of four triangle squares equals the sum of three paired-edge squares."""
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

def test_Polyform_exp():
    """For a two-edge path, exp(L)=1+L+L^2/2, with one spanning-tree term."""
    pfPath3 = qf_ab + qf_bc
    pfPathExp = quForm.one() + qf_ab + qf_bc + qf_abc
    assert pfPath3.exp() == pfPathExp

# ---------- Metric ----------

def test_Path_Potential():
    """Unit path edges give potentials 1, 2 and 3 for pairs at those distances."""
    pfPath4 = qf_ab + qf_bc + qf_cd
    pfPathExp = pfPath4.exp()
    assert pfPathExp.potential(qf_ab) == 1
    assert pfPathExp.potential(qf_ac) == 2
    assert pfPathExp.potential(qf_ad) == 3
