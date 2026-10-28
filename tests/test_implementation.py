"""Implementation invariants independent of a particular PMG metric."""
import pytest
import networkx as nx
from polyform_algebra import Point, Vector, Element, Simplex, Polysimplex, Polyform, quForm


def test_repeated_pairs_accumulate_coefficients():
    """Repeated quadratic terms add their weights; opposite weights cancel."""
    a,b=Point.create_list("ab")
    assert Polyform.from_pairs([([(a,b)],2), ([(a,b)],3)]) == 5*quForm([(a,b)])
    assert Polyform.from_pairs([([(a,b)],2), ([(b,a)],-2)]).is_zero()


@pytest.mark.parametrize("method", ["extract_max_coeff", "extract_max_coeff_all_grades", "extract_diagonal_components"])
def test_component_selection_preserves_polyform_context(method):
    """Inherited selectors must accept Polyform and retain isolated vertices."""
    graph=nx.empty_graph(3)
    graph.add_edge(0,1,weight=2)
    form=Polyform.from_networkx(graph)
    result=getattr(form,method)()
    assert isinstance(result,Polyform)
    assert result.vertices == form.vertices
    if method == "extract_diagonal_components":
        # A lone edge is not a diagonal component under the existing size criterion.
        assert result.is_zero()
    else:
        assert result == form
    empty=getattr(0*form,method)()
    assert empty.is_zero()
    assert empty.vertices == form.vertices


def test_maximum_selection_keeps_ties_and_filters_grade():
    """Absolute maxima retain ties; signed maxima select the largest coefficient."""
    a,b,c=Point.create_list("abc")
    form=Polyform.from_pairs([([(a,b)],-3), ([(b,c)],3), ([(a,b,c)],2)])
    assert form.extract_max_coeff() == -3*quForm([(a,b)])+3*quForm([(b,c)])
    assert form.extract_max_coeff(by_abs=False) == 3*quForm([(b,c)])
    assert form.extract_max_coeff(grade=2) == 2*quForm([(a,b,c)])


@pytest.mark.parametrize("cls,weight",[(Point,1),(Vector,0)])
def test_equal_elements_have_equal_hashes(cls,weight):
    """Equality is based on name and affine weight, not the concrete subclass."""
    left=cls("x")
    right=Element("x",multiplicity=weight)
    assert left == right
    assert hash(left) == hash(right)
    assert len({left,right}) == 1
    assert left+right == 2*left


def test_frame_basis_lists_stay_aligned():
    """Each coordinate representation has a matching (possibly unknown) basis."""
    a,b=Point.create_list("ab")
    frames=[a.as_polysimplex(),b.as_polysimplex()]
    point=Point("p",frame=frames)
    assert point.frames == frames
    assert point.basis_for_frame == [None,None]
    point.add_frame(frames[0])
    assert len(point.frames) == len(point.basis_for_frame) == 3
    assert len(frames) == 2


def test_polysimplex_copy_normalizes_orientation():
    """Copying and pair construction follow the same normalization path."""
    a,b=Point.create_list("ab")
    original=Polysimplex([(Simplex([b,a]),2)])
    copied=Polysimplex(original)
    assert copied == -2*(a*b)
    assert copied is not original
    assert copied.terms is not original.terms


def test_selection_preserves_custom_term_type():
    """Generic Polysimplex selectors also support explicitly chosen quadratic terms."""
    a=Point("a")
    source=Polysimplex({quForm([a]):2},term_type=quForm)
    selected=source.extract_diagonal_components()
    assert selected.term_type is quForm
    assert selected.terms == source.terms
