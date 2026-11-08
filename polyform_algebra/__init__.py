# polyform_algebra/__init__.py
# Re-export main classes for convenient access

from importlib import import_module

# Keep historic imports while allowing the reference core to load on its own.
_LEGACY = {
    'Point': '.objects.element', 'Vector': '.objects.element', 'Element': '.objects.element',
    'Simplex': '.objects.simplex', 'quForm': '.objects.quform',
    'Polysimplex': '.combinations.polysimplex', 'Polyform': '.combinations.polyform',
    'Basis': '._core.basis', 'Space': '._core.space', 'Monoid': '._core.monoid',
    'LinearCombination': '._core.lincomb',
}


def get_version():
    """Return installed distribution metadata; raise PackageNotFoundError if absent."""
    from importlib.metadata import version
    return version("polyform_algebra")


def __getattr__(name):
    """Resolve legacy classes and installed version without eager thematic imports."""
    if name == "__version__":
        return get_version()
    if name in _LEGACY:
        value = getattr(import_module(_LEGACY[name], __name__), name)
        globals()[name] = value
        return value
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


def __dir__():
    return sorted(set(globals()) | set(__all__))
