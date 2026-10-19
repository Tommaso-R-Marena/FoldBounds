"""Independent SAW enumeration and exact small-tree optimality oracles."""
from dataclasses import replace
from fractions import Fraction
from itertools import combinations
import pytest
from foldbounds.lattice import ROOT, _neighbors, CertificateError
from foldbounds.sparse_monotone import (
    SparseMonotoneCertificate, optimize_sparse_monotone, verify_sparse_monotone)

CONES=(frozenset(((1,0),(0,1))),frozenset(((1,0),(0,-1))),
       frozenset(((-1,0),(0,1))),frozenset(((-1,0),(0,-1))))

def all_walks(n):
    paths=[ROOT]
    for _ in range(2,n):
        paths=[p+((p[-1][0]+dx,p[-1][1]+dy),)
               for p in paths for dx,dy in ((1,0),(-1,0),(0,1),(0,-1))
               if (p[-1][0]+dx,p[-1][1]+dy) not in p]
    return paths

def independent_member(full,prefix):
    if full[:len(prefix)]!=prefix:return False
    moves=frozenset((full[i][0]-full[i-1][0],full[i][1]-full[i-1][1])
                    for i in range(len(prefix),len(full)))
    return any(moves<=cone for cone in CONES)

def independent_mass(cert):
    event=complement=Fraction(0)
    count=0
    for full in all_walks(len(cert.sequence)):
        if not any(independent_member(full,leaf) for leaf in cert.leaves):
            continue
        score=sum(cert.sequence[i]=='H' and cert.sequence[j]=='H'
                  and abs(full[i][0]-full[j][0])+abs(full[i][1]-full[j][1])==1
                  for i in range(len(full)) for j in range(i+2,len(full)))
        w=cert.q**score
        if score>=cert.event_at_least:event+=w
        else:complement+=w
        count+=1
    return event,complement,count

@pytest.mark.parametrize('seq',['HHHHHHHH','HPHHPHHP','PPPPPPPP','HPHHPHHPHH'])
@pytest.mark.parametrize('threshold',[0,1,2])
@pytest.mark.parametrize('budget',[1,3,8])
def test_sparse_exact_against_independent_complete_enumeration(seq,threshold,budget):
    p=optimize_sparse_monotone(seq,ROOT,Fraction(3,2),threshold,
                               leaf_budget=budget,max_core_length=min(7,len(seq)))
    checked=verify_sparse_monotone(p.certificate)
    assert (checked.event,checked.complement,checked.distinct_paths)==independent_mass(p.certificate)
    assert checked.leaf_count<=budget

def test_all_30073_complete_mixed_twelve_mer_walks():
    seq='HPHHPHHPHHHH'
    plan=optimize_sparse_monotone(seq,ROOT,Fraction(3,2),2,
                                  leaf_budget=25,max_core_length=9)
    assert len(all_walks(len(seq)))==30073
    assert plan.mass.event==Fraction(98289,32)
    assert plan.mass.complement==Fraction(1251,2)
    assert plan.mass.distinct_paths==1632
    assert (plan.mass.event,plan.mass.complement,plan.mass.distinct_paths)==independent_mass(plan.certificate)

@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(2)])
def test_reject_shadowing_and_duplicate_leaves(q):
    leaf=ROOT+((2,0),)
    c=SparseMonotoneCertificate('HHHHHH',ROOT,q,1,(ROOT,leaf))
    for selection in ((ROOT,leaf),(leaf,ROOT),(leaf,leaf)):
        with pytest.raises(CertificateError):
            verify_sparse_monotone(replace(c,leaves=selection))

@pytest.mark.parametrize('bad',[0,-1,True,1.5,None])
def test_invalid_budget_rejected(bad):
    with pytest.raises(ValueError):
        optimize_sparse_monotone('HHHH',ROOT,Fraction(2),1,
                                 leaf_budget=bad,max_core_length=4)

@pytest.mark.parametrize('obj',['event','complement'])
def test_optimality_vs_bruteforce_all_antichains(obj):
    from foldbounds.monotone import monotone_family_mass
    seq='HPHHPHHP';q=Fraction(3,2);threshold=2;depth=5;budget=3
    nodes=[ROOT]
    for i in range(2,depth):
        nodes.extend(ch for p in list(nodes) if len(p)==i for ch in _neighbors(p))
    best=Fraction(0)
    for k in range(budget+1):
        for leaves in combinations(nodes,k):
            if any(a[:len(b)]==b or b[:len(a)]==a
                   for j,a in enumerate(leaves) for b in leaves[j+1:]):
                continue
            mass=sum((getattr(monotone_family_mass(seq,p,q,threshold),obj)
                      for p in leaves),Fraction(0))
            best=max(best,mass)
    p=optimize_sparse_monotone(seq,ROOT,q,threshold,
                               leaf_budget=budget,max_core_length=depth,objective=obj)
    assert getattr(p.mass,obj)==best

@pytest.mark.parametrize('threshold',[-1,0,1,2,9])
def test_empty_antichain_is_valid_zero_lower_mass(threshold):
    m=verify_sparse_monotone(SparseMonotoneCertificate('H'*6,ROOT,Fraction(2),threshold,()))
    assert m.event==m.complement==0
    assert m.distinct_paths==0

def test_more_budget_never_weakens_event_lower():
    kw=dict(sequence='HPHHPHHPHHHH',prefix=ROOT,q=Fraction(3,2),
            event_at_least=2,max_core_length=9)
    a=optimize_sparse_monotone(**kw,leaf_budget=10)
    b=optimize_sparse_monotone(**kw,leaf_budget=25)
    assert a.mass.event<=b.mass.event
    assert b.mass.event>Fraction(4455,4)

def independent_turn_union_member(full,leaf):
    if full[:len(leaf)]!=leaf:
        return False
    moves=tuple((full[j][0]-full[j-1][0],full[j][1]-full[j-1][1])
                for j in range(len(leaf),len(full)))
    is_monotone=any(set(moves)<=cone for cone in CONES)
    two_turns=sum(a!=b for a,b in zip(moves,moves[1:]))<=2
    return is_monotone or two_turns

def independent_turn_union_mass(cert):
    event=complement=Fraction(0)
    count=0
    for full in all_walks(len(cert.sequence)):
        if not any(independent_turn_union_member(full,p) for p in cert.leaves):
            continue
        k=sum(cert.sequence[i]=='H' and cert.sequence[j]=='H'
              and abs(full[i][0]-full[j][0])+abs(full[i][1]-full[j][1])==1
              for i in range(len(full)) for j in range(i+2,len(full)))
        weight=cert.q**k
        if k>=cert.event_at_least:event+=weight
        else:complement+=weight
        count+=1
    return event,complement,count

def test_full_mixed_twelve_mer_two_turn_union_25_leaf_oracle():
    seq='HPHHPHHPHHHH';q=Fraction(3,2)
    plan=optimize_sparse_monotone(seq,ROOT,q,2,leaf_budget=25,
                                  max_core_length=8,include_nonmonotone_turns=True)
    assert plan.mass.event==Fraction(71631,16)
    assert plan.mass.complement==Fraction(1746)
    assert plan.mass.distinct_paths==2876
    assert (plan.mass.event,plan.mass.complement,plan.mass.distinct_paths)==independent_turn_union_mass(plan.certificate)

def test_turn_union_optimum_dominates_monotone_only_at_same_budget():
    kw=dict(sequence='HPHHPHHPHHHH',prefix=ROOT,q=Fraction(3,2),
            event_at_least=2,leaf_budget=25,max_core_length=8)
    plain=optimize_sparse_monotone(**kw)
    expanded=optimize_sparse_monotone(**kw,include_nonmonotone_turns=True)
    assert expanded.mass.event>plain.mass.event
    assert plain.mass.event==Fraction(98289,32)
    assert expanded.mass.event==Fraction(71631,16)

@pytest.mark.parametrize('bad',[0,1,None,'true'])
def test_turn_union_mode_requires_boolean(bad):
    kw=dict(sequence='HHHHHH',prefix=ROOT,q=Fraction(2),
            event_at_least=1,leaf_budget=2,max_core_length=4)
    with pytest.raises(ValueError):
        optimize_sparse_monotone(**kw,include_nonmonotone_turns=bad)
    good=SparseMonotoneCertificate('HHHHHH',ROOT,Fraction(2),1,(),bad)
    with pytest.raises(ValueError):
        verify_sparse_monotone(good)

def test_compact_length_64_exact_restricted_family_mass():
    # Regression for huge non-enumerated family; not an independent full 64-mer
    # enumeration oracle. Mathematical soundness is tested independently at
    # smaller lengths, and the formula counts exactly this restricted family.
    p=optimize_sparse_monotone('H'*64,ROOT,Fraction(2),1,
                               leaf_budget=25,max_core_length=6)
    assert p.candidate_prefixes==109
    assert p.mass.leaf_count==25
    assert p.mass.distinct_paths==20032011142543966146
    assert p.mass.event==66437101702969556806
    assert p.mass.complement==2882303761517117437

def test_200_leaf_union_recovers_over_half_exact_mixed_event_mass():
    seq='HPHHPHHPHHHH';q=Fraction(3,2);t=2
    p=optimize_sparse_monotone(seq,ROOT,q,t,leaf_budget=200,
                               max_core_length=8,include_nonmonotone_turns=True)
    assert p.mass.event==Fraction(421407,32)
    assert p.mass.distinct_paths==10824
    assert (p.mass.event,p.mass.complement,p.mass.distinct_paths)==independent_turn_union_mass(p.certificate)
    assert p.mass.event>Fraction(23166,2)
    comp=optimize_sparse_monotone(seq,ROOT,q,t,leaf_budget=200,
                                 max_core_length=8,objective='complement',
                                 include_nonmonotone_turns=True)
    assert comp.mass.complement==19426
    assert (comp.mass.event,comp.mass.complement,comp.mass.distinct_paths)==independent_turn_union_mass(comp.certificate)
