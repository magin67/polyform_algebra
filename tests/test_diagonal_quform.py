"""Diagonal coordinate squares versus full formal polarization.

PMG core specification, Formal forms: [X+Y]^2=[X]^2+{X,Y}+[Y]^2.
https://github.com/magin67/pmg-notes/blob/v5/pmg-library/ru/Основания/Спецификация%20ядра%20PMG.md
The diagonal map is a restricted computational helper, not a formal identity.
"""
import pytest
from polyform_algebra import Basis, Element, Point, Vector, Simplex, Polysimplex, Polyform, quForm


@pytest.fixture(params=["to_diagonal_quform", "to_quadratic_form"])
def convert(request):
    """Exercise both the descriptive API and its backward-compatible entry point."""
    return lambda obj, *args, **kwargs: getattr(obj, request.param)(*args, **kwargs)


def test_diagonal_square_retains_squared_signed_coordinates(convert):
    """Coordinates 2 and -3 produce coefficients 4 and 9, without a mixed term."""
    u,v=Vector.create_list("uv")
    x=2*u-3*v
    original=dict(x.terms)
    basis=Basis([v,u],is_orthonormal=True)
    assert convert(x,basis)==4*quForm([u])+9*quForm([v])
    assert x.terms==original
    assert basis.elements==[v,u]


def test_diagonal_part_differs_from_polarization(convert):
    """The full square includes {a,b}; merely declaring a basis does not remove it."""
    a,b=Point.create_list("ab")
    basis=Basis([a,b],is_orthonormal=True)
    x=a+b
    diagonal=convert(x,basis)
    full=x@x
    assert not diagonal.is_equivalent(full)
    mixed=Polyform.polar_form(a.as_polysimplex(),b.as_polysimplex())
    assert (full-diagonal).is_equivalent(mixed)
    assert mixed == 2*(a@b)


def test_zero_and_single_coordinate(convert):
    """Zero has zero diagonal square; a single atomic coordinate needs no mixed terms."""
    a=Point("a")
    basis=Basis([a],is_orthonormal=True)
    assert convert(Polysimplex.zero(),basis).is_zero()
    assert convert(Polysimplex.zero(),Basis([],is_orthonormal=True)).is_zero()
    assert convert(2*a,basis)==4*quForm([a])


@pytest.mark.parametrize("kind", ["scalar","higher","boundary","mixed","form","zero_form"])
def test_non_coordinate_inputs_are_rejected(convert,kind):
    """The helper only accepts external atomic grade-one coordinates (or external zero)."""
    a,b=Point.create_list("ab")
    inputs={
        "scalar":Polysimplex.one(),
        "higher":(a*b).as_polysimplex(),
        "boundary":Simplex([(a,b)]).as_polysimplex(),
        "mixed":Polysimplex.one()+a,
        "form":quForm([a]).as_polyform(),
        "zero_form":Polyform.zero(),
    }
    with pytest.raises(TypeError):
        convert(inputs[kind],Basis([a,b],is_orthonormal=True))


def test_quadratic_term_type_is_not_external_coordinates(convert):
    """A generic Polysimplex container with quadratic terms must also be rejected."""
    a=Point("a")
    x=Polysimplex({quForm([a]):2},term_type=quForm)
    with pytest.raises(TypeError,match="external coordinates"):
        convert(x,Basis([a],is_orthonormal=True))


def test_basis_membership_is_checked_without_expanding_frames(convert):
    """An element with a valid frame is still not a coordinate in a different basis."""
    a,b=Point.create_list("ab")
    v=Vector("v",frame=b-a)
    basis=Basis([a,b],is_orthonormal=True)
    with pytest.raises(ValueError,match="does not belong"):
        convert(v.as_polysimplex(),basis)
    assert convert(v[0],basis)==quForm([a])+quForm([b])


def test_equal_element_instance_is_a_valid_coordinate(convert):
    """Membership follows library equality by name and affine weight, not object identity."""
    a=Point("a")
    equal=Element("a",multiplicity=1)
    assert convert(equal.as_polysimplex(),Basis([a],is_orthonormal=True))==quForm([a]).as_polyform()


@pytest.mark.parametrize("kind,error", [("missing",TypeError),("undeclared",ValueError),
    ("duplicate",ValueError),("non_atomic",TypeError)])
def test_invalid_bases_are_rejected_even_for_zero(convert,kind,error):
    """Zero input does not bypass validation of the requested coordinate system."""
    a=Point("a")
    bases={"missing":None,"undeclared":Basis([a]),
           "duplicate":Basis([a,Point("a")],is_orthonormal=True),
           "non_atomic":Basis([Simplex.one()],is_orthonormal=True)}
    with pytest.raises(error):
        convert(Polysimplex.zero(),bases[kind])


@pytest.mark.parametrize("tolerance", [-1,float("inf"),float("nan"),True,"small"])
def test_invalid_tolerances_are_rejected(convert,tolerance):
    """Rounding requires a finite non-negative real tolerance for both API names."""
    with pytest.raises(ValueError,match="Tolerance"):
        convert(Polysimplex.zero(),Basis([],is_orthonormal=True),tolerance=tolerance)


def test_tolerance_applies_to_squared_coefficients(convert):
    """A coordinate 0.01 has squared coefficient 0.0001, below a tolerance of 0.001."""
    a=Point("a")
    basis=Basis([a],is_orthonormal=True)
    assert convert(0.01*a,basis,0.001).is_zero()
    assert convert(0.01*a,basis,1e-6)==0.0001*quForm([a])


def test_square_overflow_is_reported_explicitly(convert):
    """A finite coordinate may overflow on squaring; do not return an infinite form."""
    a=Point("a")
    with pytest.raises(ValueError,match="finite"):
        convert(1e200*a,Basis([a],is_orthonormal=True))


def test_diagonal_eigenbasis_evaluation_preserves_path_norm(convert):
    """A unit path of three edges has endpoint resistance three in either representation."""
    a,b,c,d=Point.create_list("abcd")
    lap=Polyform.from_pairs([([(a,b)],1), ([(b,c)],1), ([(c,d)],1)])
    eig,eigenbasis=lap.to_eigenbasis()
    points=eigenbasis.create_basis_from_frames()
    eigenbasis.link_to(points)
    vector=d-a
    coordinates=vector.transform(points,eigenbasis)
    diagonal=convert(coordinates,eigenbasis)
    assert lap.exp().norm(vector@vector)==pytest.approx(3)
    assert eig.exp().norm(diagonal)==pytest.approx(3)
