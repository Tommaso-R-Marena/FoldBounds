"""Independent exhaustive-oracle checks for disjoint extremal-quadrant shells."""
from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from foldbounds.lattice import ROOT, CertificateError, build_frontier, exact_spectrum, _neighbors
from foldbounds.shell_family import extremal_shell_mass

QUADRANTS = ((1,1),(1,-1),(-1,1),(-1,-1))


def independent_walks(n):
    paths=[ROOT]
    for j in range(2,n):
        paths=[child for path in paths for child in _neighbors(path)]
    return paths


def independent_member(full, prefix, core_length):
    if full[:len(prefix)] != prefix:
        return False
    if core_length == len(full):
        return True
    core=full[:core_length]
    tip=core[-1]
    for sx,sy in QUADRANTS:
        potential=lambda point: sx*point[0]+sy*point[1]
        if potential(tip) != max(map(potential,core)):
            continue
        if all(potential(full[j])-potential(full[j-1])==1
               for j in range(core_length,len(full))):
            return True
    return False


def independent_contacts(seq,path):
    return sum(
        seq[i]=='H' and seq[j]=='H'
        and abs(path[i][0]-path[j][0])+abs(path[i][1]-path[j][1])==1
        for i in range(len(path)) for j in range(i+2,len(path))
    )


@pytest.mark.parametrize('n', range(2,10))
def test_all_small_models_match_independent_enumeration(n):
    walks=independent_walks(n)
    seqs=['H'*n,'P'*n,('HP'*n)[:n],('HPHHPHHPHH'*n)[:n]]
    if n<=4:
        seqs=[''.join(x) for x in product('HP',repeat=n)]
    for core_length in range(2,n+1):
        for prefix in (ROOT, *(p for p in {w[:min(4,n)] for w in walks} if p!=ROOT)):
            if len(prefix)>core_length:
                continue
            selected=[p for p in walks if independent_member(p,prefix,core_length)]
            for sequence in seqs:
                counts={}
                for walk in selected:
                    score=independent_contacts(sequence,walk)
                    counts[score]=counts.get(score,0)+1
                for q in (Fraction(1),Fraction(3,2),Fraction(2)):
                    for threshold in (0,1,2,3):
                        actual=extremal_shell_mass(sequence,prefix,q,threshold,
                                                   core_length=core_length)
                        event=sum((freq*q**score for score,freq in counts.items()
                                   if score>=threshold),Fraction(0))
                        other=sum((freq*q**score for score,freq in counts.items()
                                   if score<threshold),Fraction(0))
                        assert (actual.event,actual.complement,actual.distinct_paths)==(
                            event,other,len(selected))


def test_large_incomplete_family_strictly_improves_event_mass():
    sequence='HPHHPHHPHHHH'
    result=extremal_shell_mass(sequence,ROOT,Fraction(3,2),2,core_length=10)
    assert result.event==6867
    assert result.event>Fraction(459,8) # previous two-turn family
    assert result.distinct_paths==17227
    assert result.seed_count==4067
    assert result.event < sum((freq*Fraction(3,2)**score for score,freq in
                               exact_spectrum(sequence).items() if score>=2),Fraction(0))


def test_all_h_18_verified_subset_mass():
    result=extremal_shell_mass('H'*18,ROOT,Fraction(2),1,core_length=12)
    assert result.event==8264272
    assert result.distinct_paths==2022389


@pytest.mark.parametrize('invalid',[True,2.0,-1,1,'5',None])
def test_core_length_validated(invalid):
    with pytest.raises(ValueError):
        extremal_shell_mass('H'*8,ROOT,Fraction(2),1,core_length=invalid)


@pytest.mark.parametrize('bad_q',[Fraction(0),Fraction(1,2),2])
def test_exact_q_only(bad_q):
    with pytest.raises(ValueError):
        extremal_shell_mass('H'*8,ROOT,bad_q,1,core_length=5)


def test_bad_path_rejected():
    with pytest.raises(CertificateError):
        extremal_shell_mass('H'*8,ROOT+((0,0),),Fraction(2),1,core_length=5)


def test_boolean_threshold_rejected():
    with pytest.raises(ValueError):
        extremal_shell_mass('H'*8,ROOT,Fraction(2),True,core_length=5)


def test_cached_calls_do_not_bypass_type_validation():
    valid=extremal_shell_mass('H'*8,ROOT,Fraction(2),1,core_length=5)
    assert valid.distinct_paths>0
    with pytest.raises(ValueError):
        extremal_shell_mass('H'*8,ROOT,Fraction(2),True,core_length=5)
    with pytest.raises(ValueError):
        extremal_shell_mass('H'*8,ROOT,Fraction(2),1,core_length=True)