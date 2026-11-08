# polyform_algebra/__init__.py
# Re-export main classes for convenient access

from .objects.element import Point, Vector, Element
from .objects.simplex import Simplex
from .objects.quform import quForm
from .combinations.polysimplex import Polysimplex
from .combinations.polyform import Polyform
from ._core.basis import Basis
from ._core.space import Space
from ._core.monoid import Monoid
from ._core.lincomb import LinearCombination


def get_version():
    """Return installed distribution metadata; raise PackageNotFoundError if absent."""
    from importlib.metadata import version
    return version("polyform_algebra")


def __getattr__(name):
    """Expose version lazily so an uninstalled checkout can still import classes."""
    if name == "__version__":
        return get_version()
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    'get_version', '__version__',
    'Point', 'Vector', 'Element',
    'Simplex',
    'quForm',
    'Polysimplex',
    'Polyform',
    'Basis',
    'Space',
    'Monoid',
    'LinearCombination',
]
