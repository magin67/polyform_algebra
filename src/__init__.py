# polyform_algebra/__init__.py
# Lazy imports to avoid circular dependencies

def __getattr__(name):
    """Lazy loading of modules to avoid circular imports"""
    _lazy_imports = {
        'Point': ('objects.element', 'Point'),
        'Vector': ('objects.element', 'Vector'),
        'Element': ('objects.element', 'Element'),
        'Simplex': ('objects.simplex', 'Simplex'),
        'quForm': ('objects.quform', 'quForm'),
        'Polysimplex': ('combinations.polysimplex', 'Polysimplex'),
        'Polyform': ('combinations.polyform', 'Polyform'),
        'Basis': ('_core.basis', 'Basis'),
        'Monoid': ('_core.monoid', 'Monoid'),
        'LinearCombination': ('_core.lincomb', 'LinearCombination'),
    }
    
    if name in _lazy_imports:
        module_name, attr_name = _lazy_imports[name]
        from importlib import import_module
        module = import_module(f'.{module_name}', __name__)
        return getattr(module, attr_name)
    
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

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
