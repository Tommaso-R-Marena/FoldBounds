"""Tight Cartesian-box geometry, strict cutoff and exact full-world GNM tests."""
from fractions import Fraction
from dataclasses import replace

import pytest
pytest.importorskip('sympy')

from foldbounds.gnm_joint_rank import classify_contact_boxes
from foldbounds.gnm_exact_coordinate_boxes import (
    exact_coordinate_box_graphs,audit_exact_box_gnm_ranks,
)
from experiments.real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca


def test_exact_cartesian_box_minmax_strict_cutoff():
    # Residues indexed 0 and 3: one nonlocal pair. Each coordinate can move
    # independently by ±50 milliangstrom. Their relative x-coordinate
    # therefore ranges over d±100, an EXACT attainable interval.
    def edge_status(d):
        coords=((0,0,0),(50000,0,0),(50000,0,0),(d,0,0))
        graph=exact_coordinate_box_graphs(coords,Fraction(8),Fraction(1,20))
        return ((0,3) in graph.certainly_on,(0,3) in graph.possibly_on)
    assert edge_status(7800)==(True,True)
    assert edge_status(7900)==(False,True)  # maximum exactly 8000, no guarantee
    assert edge_status(8100)==(False,False) # minimum exactly 8000, impossible contact
    assert edge_status(8150)==(False,False)
    assert edge_status(8000)==(False,True)
    with pytest.raises(ValueError):
        exact_coordinate_box_graphs(((0,0,0),(1,0,0)),Fraction(8),0.05)
    with pytest.raises(ValueError):
        exact_coordinate_box_graphs(((0,0,0),(1,0,0)),Fraction(0),Fraction(1,20))
    with pytest.raises(ValueError):
        exact_coordinate_box_graphs(((0,0,0),(1,0,0)),Fraction(8),Fraction(-1,20))
    with pytest.raises(ValueError):
        exact_coordinate_box_graphs(((0,0,False),(1,0,0)),Fraction(8),Fraction(1,20))


def test_exact_1d3z_box_is_nested_inside_coarse_conservative_graph():
    coords=load_ca(DEFAULT_FIXTURE)
    coarse_lo,coarse_hi=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    new=exact_coordinate_box_graphs(coords,Fraction(8),Fraction(1,20))
    assert len(coords)==69
    assert set(coarse_lo)<=set(new.certainly_on)
    assert set(new.possibly_on)<=set(coarse_hi)
    assert len(coarse_lo)==len(new.certainly_on)==68+152
    assert len(coarse_hi)==68+167
    assert len(new.possibly_on)==68+164
    assert len(new.per_pair_optional)==12
    removed=set(coarse_hi)-set(new.possibly_on)
    assert removed=={(23,50),(41,49),(48,57)}
    # Exact pairwise box min² is already above 8Å squared for each.
    expected={
        (23,50):64514030,
        (41,49):64585089,
        (48,57):64223700,
    }
    for (i,j),min_sq in expected.items():
        r=[abs(coords[i][k]-coords[j][k]) for k in range(3)]
        assert sum(max(0,x-100)**2 for x in r)==min_sq
        assert min_sq>8000**2


def test_real_1d3z_complete_exact_pairwise_box_rank_atlas():
    coords=load_ca(DEFAULT_FIXTURE)
    a=exact_coordinate_box_graphs(coords,Fraction(8),Fraction(1,20))
    result=audit_exact_box_gnm_ranks(a)
    assert (result.n,result.certain_nonlocal,result.possible_nonlocal)==(69,152,164)
    assert (result.optional_contacts,result.exact_graph_worlds)==(12,4096)
    assert (result.globally_invariant,result.opposite_sign_witnessed)==(2096,250)
    assert (result.noninvariant_with_ties_only,result.always_tied)==(0,0)
    assert result.reduced_inverses_checked==4096
    assert len(result.flipping_witnesses)==250
    assert all(0<=p<4096 and 0<=q<4096 and p!=q
               for i,j,p,q in result.flipping_witnesses)


@pytest.mark.parametrize('epsilon,expected_optional,expected_stable,expected_reverse',[
    (Fraction(0),0,2346,0),
    (Fraction(1,100),4,2268,78),
    (Fraction(1,40),7,2220,126),
])
def test_real_1d3z_tighter_box_across_other_uncertainty_widths(
        epsilon,expected_optional,expected_stable,expected_reverse):
    coords=load_ca(DEFAULT_FIXTURE)
    a=exact_coordinate_box_graphs(coords,Fraction(8),epsilon)
    r=audit_exact_box_gnm_ranks(a)
    assert r.optional_contacts==expected_optional
    assert r.exact_graph_worlds==(1<<expected_optional)
    assert (r.globally_invariant,r.opposite_sign_witnessed)==(
        expected_stable,expected_reverse)
    assert r.noninvariant_with_ties_only==r.always_tied==0


def test_resource_and_forged_graph_fail_closed():
    coords=load_ca(DEFAULT_FIXTURE)
    a=exact_coordinate_box_graphs(coords,Fraction(8),Fraction(1,20))
    with pytest.raises(ValueError,match='budget'):
        audit_exact_box_gnm_ranks(a,max_worlds=2048)
    with pytest.raises(ValueError,match='optional contacts'):
        audit_exact_box_gnm_ranks(replace(a,per_pair_optional=()))
    with pytest.raises(ValueError,match='mandatory'):
        audit_exact_box_gnm_ranks(replace(a,possibly_on=()))
