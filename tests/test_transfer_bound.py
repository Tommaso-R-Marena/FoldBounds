"""Independent exact-oracle checks for prefix-obstacle transfer bounds."""
from collections import defaultdict
from fractions import Fraction
from itertools import product
import pytest

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead
from foldbounds.transfer_bound import nonbacktracking_mass_upper, relaxed_cross_spectrum
from foldbounds.adaptive import adaptive_refine, verify_refinement_run


def independent_walks(n):
    paths=[((0,0),(1,0))]
    for _ in range(2,n):
        result=[]
        for path in paths:
            x,y=path[-1]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                p=(x+dx,y+dy)
                if p not in path:
                    result.append(path+(p,))
        paths=result
    return paths


def independent_contacts(sequence,path):
    return sum(sequence[i]=='H' and sequence[j]=='H' and
               abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
               for i in range(len(sequence))
               for j in range(i+2,len(sequence)))


@pytest.mark.parametrize('n',range(2,8))
def test_exhaustive_prefix_upper_enclosures(n):
    # No production contact scoring used by this complete oracle.
    full=independent_walks(n)
    for seq_tuple in product('HP',repeat=n):
        sequence=''.join(seq_tuple)
        by_prefix=defaultdict(lambda:defaultdict(int))
        for walk in full:
            k=independent_contacts(sequence,walk)
            for m in range(2,n):
                by_prefix[walk[:m]][k]+=1
        for path,histogram in by_prefix.items():
            for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                z=sum((num*q**k for k,num in histogram.items()),Fraction(0))
                for threshold in (-1,0,1,2,4):
                    b=nonbacktracking_mass_upper(sequence,path,q,threshold)
                    e=sum((num*q**k for k,num in histogram.items()
                           if k>=threshold),Fraction(0))
                    assert z<=b.partition
                    assert e<=b.event
                    assert z-e<=b.complement


@pytest.mark.parametrize('sequence',['HHHH','HHHHHH','HHHHHHHH',
                                     'HPHHPHHP','HHPHHHPHHH'])
@pytest.mark.parametrize('q',[Fraction(1),Fraction(3,2),Fraction(2)])
@pytest.mark.parametrize('threshold',[0,1,2,4])
@pytest.mark.parametrize('horizon',[0,2])
def test_integrated_transfer_never_weakens_existing_bounds(sequence,q,threshold,horizon):
    cert=build_frontier(sequence,q,threshold,2)
    old=verify_lookahead(cert,lookahead=horizon,use_matching=True)
    new=verify_lookahead(cert,lookahead=horizon,use_matching=True,use_transfer=True)
    spectrum=exact_spectrum(sequence)
    z=sum((num*q**k for k,num in spectrum.items()),Fraction(0))
    e=sum((num*q**k for k,num in spectrum.items() if k>=threshold),Fraction(0))
    assert old.partition_lower<=new.partition_lower<=z<=new.partition_upper<=old.partition_upper
    assert old.event_lower<=new.event_lower<=new.event_upper<=old.event_upper
    assert old.complement_lower<=new.complement_lower<=new.complement_upper<=old.complement_upper
    assert old.probability_lower<=new.probability_lower<=e/z<=new.probability_upper<=old.probability_upper


def test_four_residue_reference_is_exact():
    b=nonbacktracking_mass_upper('HHHH',((0,0),(1,0)),Fraction(2),1)
    assert b.partition==11 and b.event==4 and b.complement==7
    assert b.future_future_cap==0


def test_no_room_prefix_gives_zero_upper_mass():
    prefix=((0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(0,1),(1,1))
    b=nonbacktracking_mass_upper('H'*10,prefix,Fraction(2),1)
    assert b.partition==b.event==b.complement==0


def test_exact_full_length_and_cache_are_deterministic():
    path=((0,0),(1,0),(1,1),(0,1))
    a=nonbacktracking_mass_upper('HHHH',path,Fraction(3,2),1)
    assert a.partition==a.event==Fraction(3,2)
    assert a.complement==0
    p=((0,0),(1,0))
    assert relaxed_cross_spectrum('HHHHHH',p)==relaxed_cross_spectrum('HHHHHH',p)


@pytest.mark.parametrize('bad',[Fraction(0),Fraction(1,2),2.0])
def test_reject_invalid_q(bad):
    with pytest.raises(ValueError):
        nonbacktracking_mass_upper('HHHH',((0,0),(1,0)),bad,1)


@pytest.mark.parametrize('policy',['breadth','mass','probability'])
def test_adaptive_replay_captures_transfer_selection(policy):
    cert=build_frontier('HHPHHHPH',Fraction(2),1,2)
    run=adaptive_refine(cert,steps=5,use_matching=True,use_transfer=True,
                        lookahead=0,policy=policy)
    assert run.use_transfer is True
    assert verify_refinement_run(run)==run.bounds
    spectrum=exact_spectrum(cert.sequence)
    z=sum((num*cert.q**k for k,num in spectrum.items()),Fraction(0))
    e=sum((num*cert.q**k for k,num in spectrum.items()
           if k>=cert.event_at_least),Fraction(0))
    assert run.bounds.partition_lower<=z<=run.bounds.partition_upper
    assert run.bounds.probability_lower<=e/z<=run.bounds.probability_upper



def test_invalid_prefix_is_rejected():
    with pytest.raises(ValueError):
        nonbacktracking_mass_upper('HHHH',((0,0),(1,0),(0,0)),Fraction(2),1)
