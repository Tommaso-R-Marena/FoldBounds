"""Exact independent Boolean-world oracle and adversarial tree attack tests."""
from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import json

import pytest
pytest.importorskip('sympy')

from foldbounds.gnm_joint_rank import (
    IntegerInverse, scaled_grounded_laplacian, untrusted_propose_inverse,
)
from foldbounds.gnm_reduced_woodbury import (
    ReducedWorldCache, _verify_small_inverse, exact_joint_sign,
)
from foldbounds.gnm_complete_worlds import verify_complete_contact_worlds


def independent_fraction_inverse(matrix):
    """No SymPy, no Woodbury: independent exact Fraction Gaussian elimination."""
    n=len(matrix)
    a=[[Fraction(matrix[i][j]) for j in range(n)] +
       [Fraction(int(i==j)) for j in range(n)] for i in range(n)]
    for k in range(n):
        pivot=next((i for i in range(k,n) if a[i][k]),None)
        if pivot is None:raise ArithmeticError('singular')
        a[pivot],a[k]=a[k],a[pivot]
        d=a[k][k]
        a[k]=[value/d for value in a[k]]
        for i in range(n):
            if i!=k:
                factor=a[i][k]
                if factor:
                    a[i]=[x-factor*y for x,y in zip(a[i],a[k])]
    return tuple(tuple(row[n:]) for row in a)


SMALL_N=7
SMALL_CERTAIN=((0,1),(0,2),(1,2),(1,3),(2,3),(2,5),(3,4),(4,5),(5,6))
SMALL_POSSIBLE=tuple(sorted(set(SMALL_CERTAIN)|{(0,6),(1,4),(1,6)}))


def build_small_certificate():
    n=SMALL_N
    cache=ReducedWorldCache(n,SMALL_CERTAIN,SMALL_POSSIBLE)
    opt=cache.optional
    matrix_signs={pair:set() for pair in combinations(range(n),2)}
    for mask in range(1<<len(opt)):
        edges=tuple(sorted(set(SMALL_CERTAIN)|{
            edge for j,edge in enumerate(opt) if mask&(1<<j)}))
        inverse=independent_fraction_inverse(scaled_grounded_laplacian(n,edges))
        covariance=cache.get(mask)
        for i in range(n):
            for j in range(n):
                exact=covariance.entry_numerator(i,j)
                assert Fraction(exact,cache.base.denominator*covariance.delta)==inverse[i][j]
        for i,j in combinations(range(n),2):
            difference=inverse[i][i]-inverse[j][j]
            matrix_signs[(i,j)].add((difference>0)-(difference<0))
    stable=[];variable=[]
    for (i,j),signs in matrix_signs.items():
        if signs=={1} or signs=={-1}:
            stable.append([i,j,next(iter(signs))])
        elif {-1,1}<=signs:
            variable.append([i,j])
        else:
            raise AssertionError('synthetic test has unclassifiable equality-only world')
    # The verifier independently re-proves no-split root bounds, so the
    # proposal must contain ONLY stable signs absent from that baseline.
    root_low=cache.get(0)
    root_high=cache.get((1<<len(opt))-1)
    extra_stable=[(i,j,sgn) for i,j,sgn in stable
                  if exact_joint_sign(root_low,root_high,i,j)==0]
    assert stable and extra_stable and variable
    # Full tree explicitly enumerates all optional worlds (a conservative
    # untrusted proposal), unlike the much shorter 1D3Z rank-two tree.
    def full(index=0):
        if index==len(opt):return {'leaf':True}
        return {'split':index,'zero':full(index+1),'one':full(index+1)}
    return {
        'schema':'foldbounds-gnm-full-contact-worlds-v1',
        'optional_edges':[list(e) for e in opt],
        'extra_stable_orders':[list(x) for x in extra_stable],
        'changing_pairs':variable,
        'witness_masks':list(range(1<<len(opt))),
        'tree':full(),
    }


def test_small_all_graphs_vs_independent_fraction_oracle():
    proposal=build_small_certificate()
    result=verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,
                                          proposal,selected_prior=())
    assert result.all_worlds==8
    assert result.tree_leaves==8
    assert result.invariant_pairs+result.reversing_pairs==21
    assert result.reversing_pairs==12
    assert result.tree_nodes==15


def test_small_malformed_and_forged_proofs_fail_closed(monkeypatch):
    original=build_small_certificate()
    bad=deepcopy(original)
    bad['tree']['zero']={'leaf':True}
    with pytest.raises(ArithmeticError,match='unjustified cut'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    bad=deepcopy(original)
    bad['tree']['one']=None
    with pytest.raises(ValueError,match='tree node'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    bad=deepcopy(original)
    bad['extra_stable_orders'][0][2]*=-1
    with pytest.raises(ArithmeticError,match='unjustified cut'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    bad=deepcopy(original)
    bad['changing_pairs'][0]=bad['extra_stable_orders'][0][:2]
    with pytest.raises(ValueError,match='cannot be invariant'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    bad=deepcopy(original)
    bad['witness_masks']=[0]
    with pytest.raises(ArithmeticError,match='both signs'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    bad=deepcopy(original)
    bad['optional_edges'].reverse()
    with pytest.raises(ValueError,match='geometry mismatch'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,bad)
    with pytest.raises(ValueError,match='tree budget'):
        verify_complete_contact_worlds(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE,
                                       original,max_tree_nodes=2)
    from foldbounds import gnm_reduced_woodbury as w
    saved=w._small_inverse_proposal
    def corrupt(M):
        N,d=saved(M)
        return (tuple(tuple(v+(i==j==0) for j,v in enumerate(row))
                      for i,row in enumerate(N)),d)
    monkeypatch.setattr(w,'_small_inverse_proposal',corrupt)
    cache=ReducedWorldCache(SMALL_N,SMALL_CERTAIN,SMALL_POSSIBLE)
    with pytest.raises(ArithmeticError,match='failed exact integer check'):
        cache.get(1)


def test_exact_rank_two_predicate_against_all_small_worlds():
    n=SMALL_N
    cache=ReducedWorldCache(n,SMALL_CERTAIN,SMALL_POSSIBLE)
    lo=cache.get(0);hi=cache.get((1<<len(cache.optional))-1)
    for i,j in combinations(range(n),2):
        sign=exact_joint_sign(lo,hi,i,j)
        if sign:
            for mask in range(1<<len(cache.optional)):
                v=cache.get(mask)
                a=v.diagonal_numerators()
                assert (a[i]-a[j])*sign>0


def test_real_1d3z_exact_complete_boolean_atlas():
    from foldbounds.gnm_joint_rank import classify_contact_boxes
    from experiments.real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca
    n=69
    coordinates=load_ca(DEFAULT_FIXTURE)
    assert len(coordinates)==n
    certain,possible=classify_contact_boxes(coordinates,Fraction(8),Fraction(1,20))
    path=Path(__file__).resolve().parents[1]/'research/gnm_1d3z_complete_worlds.json'
    proposal=json.loads(path.read_text(encoding='utf-8'))
    result=verify_complete_contact_worlds(n,certain,possible,proposal,
                  selected_prior=((41,49),(37,40),(5,11)))
    assert result.optional_contacts==15
    assert result.all_worlds==32768
    assert result.original_certified==1947
    assert result.extra_certified==33
    assert result.invariant_pairs==1980
    assert result.reversing_pairs==366
    assert result.tree_nodes==175
    assert result.tree_leaves==88
    assert result.tree_max_depth==10
    assert result.counterexample_networks==19
    assert result.exact_inverse_graphs==204


def test_real_1d3z_reduced_covariances_vs_direct_integer_inverse():
    # Independent full 69x69 inverse for a few fully selected graphs.
    from foldbounds.gnm_joint_rank import classify_contact_boxes
    from experiments.real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca
    coords=load_ca(DEFAULT_FIXTURE)
    certain,possible=classify_contact_boxes(coords,Fraction(8),Fraction(1,20))
    cache=ReducedWorldCache(len(coords),certain,possible)
    for mask in (0,1,30714,32767):
        edges=tuple(sorted(set(certain)|{
            e for b,e in enumerate(cache.optional) if mask&(1<<b)}))
        direct=untrusted_propose_inverse(len(coords),edges)
        world=cache.get(mask)
        for i in range(69):
            actual=world.entry_numerator(i,i)
            assert (actual*direct.denominator ==
                    direct.numerator[i][i]*world.delta*cache.base.denominator)
        for i,j in ((0,1),(3,55),(22,44),(40,60)):
            actual=world.entry_numerator(i,j)
            assert (actual*direct.denominator ==
                    direct.numerator[i][j]*world.delta*cache.base.denominator)

@pytest.mark.parametrize('epsilon,expected_optional,expected_stable,expected_reverse',[
    (Fraction(0),0,2346,0),
    (Fraction(1,100),4,2268,78),
    (Fraction(1,40),7,2220,126),
])
def test_complete_exact_world_enumeration_at_independent_uncertainty_widths(
    epsilon,expected_optional,expected_stable,expected_reverse,
):
    """Direct exact full enumeration of smaller contact uncertainty families.

    These uncertainties were not selected by fitting any NMR-motion labels.
    """
    from foldbounds.gnm_joint_rank import classify_contact_boxes
    from experiments.real_1d3z_joint_rank import DEFAULT_FIXTURE,load_ca
    coords=load_ca(DEFAULT_FIXTURE)
    certain,possible=classify_contact_boxes(coords,Fraction(8),epsilon)
    cache=ReducedWorldCache(len(coords),certain,possible)
    assert len(cache.optional)==expected_optional
    n=len(coords)
    positives=[False]*(n*(n-1)//2)
    negatives=[False]*len(positives)
    tied=[False]*len(positives)
    for mask in range(1<<len(cache.optional)):
        diag=cache.get(mask).diagonal_numerators()
        k=0
        for i in range(n):
            for j in range(i+1,n):
                d=diag[i]-diag[j]
                if d>0:positives[k]=True
                elif d<0:negatives[k]=True
                else:tied[k]=True
                k+=1
    stable=sum((p and not neg and not eq) or (neg and not p and not eq)
               for p,neg,eq in zip(positives,negatives,tied))
    reversing=sum(p and neg for p,neg in zip(positives,negatives))
    assert sum(tied)==0
    assert stable==expected_stable
    assert reversing==expected_reverse
    assert stable+reversing==2346
