"""Independent 2D HP geometry/thermodynamic checks for mixed-square witnesses."""
from fractions import Fraction
from itertools import product
import random
import pytest
from foldbounds.ground_state_transfer import build_ground_state_certificate
from foldbounds.heterogeneous_ground import (
    compact_hp_score_floor, certify_mixed_hp_population,
)

STEPS=((1,0),(0,1),(-1,0),(0,-1))


def _all_square_paths(side):
    cells={(x,y) for x in range(side) for y in range(side)}
    found=[]
    for start in cells:
        east=(start[0]+1,start[1])
        if east not in cells:continue
        path=[start,east]
        used=set(path)
        def dfs():
            if len(path)==side*side:
                found.append(tuple(path));return
            x,y=path[-1]
            for dx,dy in STEPS:
                p=(x+dx,y+dy)
                if p in cells and p not in used:
                    used.add(p);path.append(p);dfs();path.pop();used.remove(p)
        dfs()
    return found


def _score(sequence,path):
    return sum(
        sequence[i]==sequence[j]=='H'
        and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
        for i in range(len(path)) for j in range(i+2,len(path))
    )


def _all_saw_spectrum(sequence):
    n=len(sequence)
    path=[(0,0),(1,0)]
    used=set(path)
    histogram={}
    def dfs():
        if len(path)==n:
            score=_score(sequence,path)
            histogram[score]=histogram.get(score,0)+1
            return
        x,y=path[-1]
        for dx,dy in STEPS:
            p=(x+dx,y+dy)
            if p not in used:
                used.add(p);path.append(p);dfs();path.pop();used.remove(p)
    dfs()
    return histogram


def test_independent_hamiltonian_contact_loss_all_patterns_small():
    for side in (2,3,4):
        paths=_all_square_paths(side)
        rng=random.Random(side)
        if side==2:
            sequences=[''.join(s) for s in product('HP',repeat=side**2)]
        else:
            sequences=['H'*(side**2),'P'*(side**2)]
            sequences+=[''.join(rng.choice('HP') for _ in range(side**2))
                        for _ in range(30)]
        for sequence in sequences:
            minimum=compact_hp_score_floor(sequence,side)
            assert all(_score(sequence,path)>=minimum for path in paths)


def test_mixed_event_probabilities_contain_independent_exact_oracle():
    cert=build_ground_state_certificate(3)
    for sequence in ('HHHHHHHHH','HPHHHHHHH','HHHHHPHHH',
                     'PPHHHHHHH','HHHHHHHPP'):
        floor=compact_hp_score_floor(sequence,3)
        if floor==0:continue
        spectrum=_all_saw_spectrum(sequence)
        for q in (Fraction(1),Fraction(3,2),Fraction(4),Fraction(10)):
            z=sum((c*q**k for k,c in spectrum.items()),Fraction(0))
            for threshold in range(1,floor+1):
                actual=sum((c*q**k for k,c in spectrum.items()
                            if k>=threshold),Fraction(0))/z
                result=certify_mixed_hp_population(cert,sequence,q,threshold)
                assert result.probability_lower<=actual


def test_64_residue_examples_can_be_checked_without_full_replay():
    # The geometry bound itself does not depend on enumerating 64-site walks.
    sequence=list('H'*64)
    for i in (10,20):sequence[i]='P'
    assert compact_hp_score_floor(''.join(sequence),8)==45
    for i in (30,40):sequence[i]='P'
    assert compact_hp_score_floor(''.join(sequence),8)==41


def test_malformed_sequences_thresholds_and_corrections_fail_closed():
    c=build_ground_state_certificate(3)
    with pytest.raises(ValueError):compact_hp_score_floor('H'*8,3)
    with pytest.raises(ValueError):compact_hp_score_floor('X'*9,3)
    with pytest.raises(ValueError):compact_hp_score_floor('H'*9,True)
    with pytest.raises(ValueError):certify_mixed_hp_population(c,'P'*9,Fraction(10),1)
    with pytest.raises(ValueError):certify_mixed_hp_population(c,'H'*9,Fraction(1,2),1)
    with pytest.raises(ValueError):certify_mixed_hp_population(c,'H'*9,Fraction(10),5)
    with pytest.raises(ValueError):certify_mixed_hp_population(c,'H'*9,Fraction(10),1,
        correction_lower=Fraction(5),correction_upper=Fraction(2))
