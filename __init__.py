# polyform_algebra/__init__.py
# Re-export main classes for convenient access

from objects.element import Point, Vector, Element
from objects.simplex import Simplex
from objects.quform import quForm
from combinations.polysimplex import Polysimplex
from combinations.polyform import Polyform
from _core.basis import Basis
from _core.monoid import Monoid
from _core.lincomb import LinearCombination

__all__ = [
    'Point', 'Vector', 'Element',
    'Simplex',
    'quForm',
    'Polysimplex',
    'Polyform',
    'Basis',
    'Monoid',
    'LinearCombination',
]