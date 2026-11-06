"""Independent full-SA W oracle and adversarial proof replay."""
from fractions import Fraction
from random import Random
import pytest
from dataclasses import replace
from foldbounds.horizon_certificate import (
    ROOT, HorizonLeaf, HorizonCertificate, build_horizon_certificate,
    verify_horizon_certificate, _component_possible)

STEPS=((1,0),(-1,0),(0,1),(0,-1))
def walks(n):
    result=[ROOT]
    for _ in range(2,n):
        children=[]
        for path in result:
            for dx,dy in STEPS:
                p=(path[-1][0]+dx,path[-1][1]+dy)
                if p not in path:children.append(path+(p,))
        result=children
    return result

def score(seq,path):
    return sum(seq[i]=='H' and seq[j]=='H' and j-i>=2 and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(path)) for j in range(i+1,len(path)))

@pytest.mark.parametrize('n',range(2,10))
def test_exhaustive_spectra_and_all_horizons(n):
    allwalks=walks(n)
    for seq in ('H'*n,('HP'*n)[:n],('PH'*n)[:n],'P'*n):
        for prefix in (ROOT,allwalks[len(allwalks)//2][:min(5,n)]):
            legal=[p for p in allwalks if p[:len(prefix)]==prefix]
            for h in (0,2,n):
                c=build_horizon_certificate(seq,prefix,horizon=h)
                for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                    b=verify_horizon_certificate(c,q,1)
                    z=sum((q**score(seq,w) for w in legal),Fraction())
                    e=sum((q**score(seq,w) for w in legal if score(seq,w)>=1),Fraction())
                    assert b.exact_mass<=z<=b.partition_upper
                    assert b.exact_event<=e<=b.event_upper
                    assert b.exact_complement<=z-e<=b.complement_upper
                    if h>=n-len(prefix):
                        assert b.partition_upper==z and b.exact_mass==z


def test_random_prefixes_and_unfinished_budgets():
    rng=Random(7729)
    for _ in range(100):
        n=rng.randrange(6,11); seq=''.join(rng.choice('HP') for _ in range(n))
        allwalks=walks(n); pre=rng.choice(allwalks)[:rng.randrange(2,n)]
        cert=build_horizon_certificate(seq,pre,horizon=6,max_leaves=2)
        b=verify_horizon_certificate(cert,Fraction(3,2),rng.randrange(0,4))
        full=[w for w in allwalks if w[:len(pre)]==pre]
        total=sum((Fraction(3,2)**score(seq,w) for w in full),Fraction())
        assert b.exact_mass<=total<=b.partition_upper


def test_adversarial_forgery_rejected():
    cert=build_horizon_certificate('HHHHHHH',horizon=3)
    assert len(cert.leaves)>1
    e=cert.leaves[0]
    for malformed in (
        replace(cert,leaves=cert.leaves[1:]),
        replace(cert,leaves=cert.leaves+(e,)),
        replace(cert,leaves=(replace(e,count_upper=0),)+cert.leaves[1:]),
        replace(cert,leaves=(replace(e,score_upper=-1),)+cert.leaves[1:]),
        replace(cert,leaves=(replace(e,kind='dead',count_upper=0),)+cert.leaves[1:]),
    ):
        with pytest.raises((ValueError,AssertionError)):
            verify_horizon_certificate(malformed,Fraction(2))


def test_horizon_monotone_cap_for_trap_rich_prefixes():
    n=12
    for prefix in walks(8)[::47][:30]:
        before=verify_horizon_certificate(build_horizon_certificate('H'*n,prefix,horizon=0),Fraction(2))
        after=verify_horizon_certificate(build_horizon_certificate('H'*n,prefix,horizon=4),Fraction(2))
        # Taking min of two separately certified upper bounds is always sound.
        assert min(before.partition_upper,after.partition_upper)>0 or after.dead_leaves>0


def test_parity_component_obstruction_finds_nontrivial_deadends():
    # Explicitly look for a deep geometric trap; it need not be immediate.
    seen=0
    for w in walks(10):
        if not _component_possible(w,14):
            seen+=1
            cert=build_horizon_certificate('H'*14,w,horizon=3)
            b=verify_horizon_certificate(cert,Fraction(2))
            assert b.partition_upper==0
    assert seen>0



def test_depth_selection_never_weakens_credible_bound():
    from foldbounds.horizon_certificate import select_best_horizon_certificate
    for n in (12,20,30):
        seq='H'*n; q=Fraction(2)
        root=verify_horizon_certificate(build_horizon_certificate(seq,horizon=0),q)
        for max_h in (2,4,6,8):
            cert=select_best_horizon_certificate(seq,q=q,horizon=max_h)
            upper=verify_horizon_certificate(cert,q).partition_upper
            assert upper<=root.partition_upper
            assert cert.horizon<=max_h


def test_exact_nodes_are_not_forgeable():
    c=build_horizon_certificate('HHHHH',horizon=3)
    e=next(x for x in c.leaves if x.kind=='exact')
    for bad in (replace(e,count_upper=2),replace(e,score_lower=-1),replace(e,score_upper=3)):
        altered=replace(c,leaves=tuple(bad if x==e else x for x in c.leaves))
        with pytest.raises(ValueError):verify_horizon_certificate(altered,Fraction(2))
