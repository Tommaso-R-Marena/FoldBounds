"""Adversarial and independent 3D inverse-design tests.

The full n=16 three-cut replay is an opt-in integration because its 26M-node
threshold-10 checker is expensive. Small-chain tests always re-check all
three cut trees and independently enumerate all cubic SAWs without pruning.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import product
import os
from pathlib import Path
import random
import shutil
import subprocess

import pytest

from foldbounds.cubic_cut_certificate import CheckedNativeGap, contact_score, certify_native_gap
from foldbounds.cubic_contact_radius_sharpness import prove_sharp_contact_error_radius
from foldbounds.cubic_robust_inverse_design import (
    _d2, _hh_mask, _transverse_orbit, cubic_candidate_pairs,
    verify_ground_map_cut, design_robust_unique_orbit,
    validate_complete_robust_design,
)

SEQUENCE = "HHHPPHHHHHHHPPHH"
EXCITED_WITNESS = (
    (0,0,0),(1,0,0),(2,0,0),(3,0,0),
    (3,1,0),(2,1,0),(1,1,0),(0,1,0),
    (0,1,1),(1,1,1),(2,1,1),(3,1,1),
    (3,0,1),(2,0,1),(1,0,1),(0,0,1),
)

@pytest.fixture(scope='module')
def compilers(tmp_path_factory):
    if not shutil.which('g++'):
        pytest.skip('independent C++17 proof checking requires g++')
    source=Path(__file__).resolve().parents[1]/'experiments'
    output=tmp_path_factory.mktemp('cubic_inverse_checks')
    tools={}
    for name in ('cubic_cut_producer','cubic_cut_checker','cubic_ground_map_checker'):
        dst=output/name
        subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-Werror',
                        str(source/(name+'.cpp')),'-o',str(dst)],
                       check=True,timeout=45,capture_output=True,text=True)
        tools[name]=dst
    return tools


def _cut(producer,output,sequence,threshold):
    subprocess.run([str(producer),sequence,str(threshold),str(output),'geometry'],
                   check=True,timeout=60,capture_output=True,text=True)
    return output


@pytest.fixture(scope='module')
def small_proof(compilers,tmp_path_factory):
    output=tmp_path_factory.mktemp('small_cubic_inverse')
    paths=tuple(_cut(compilers['cubic_cut_producer'],output/f'c{t}.cut',"H"*8,t)
                for t in (4,5,6))
    excited=((0,0,0),(1,0,0),(2,0,0),(3,0,0),
             (3,1,0),(2,1,0),(1,1,0),(0,1,0))
    assert contact_score('H'*8,excited)==3
    design=validate_complete_robust_design(
        paths,compilers['cubic_cut_checker'],compilers['cubic_ground_map_checker'],
        excited,max_seconds=60
    )
    return paths,design


def _all_cubic_walks(n):
    """Independent full-walk DFS, without cut trees or geometry symmetry."""
    path=[(0,0,0),(1,0,0)]
    occupied=set(path)
    steps=((0,0,-1),(0,-1,0),(-1,0,0),(0,0,1),(0,1,0),(1,0,0))
    def walk():
        if len(path)==n:
            yield tuple(path)
            return
        x,y,z=path[-1]
        for dx,dy,dz in steps:
            nxt=(x+dx,y+dy,z+dz)
            if nxt not in occupied:
                occupied.add(nxt);path.append(nxt)
                yield from walk()
                path.pop();occupied.remove(nxt)
    yield from walk()


def test_small_3d_certificate_proves_unique_orbit_and_algebra(small_proof):
    _,design=small_proof
    assert design.sequence=='H'*8
    assert (design.ground_contact_count,design.ground_conformations,design.gap)==(5,24,2)
    assert design.target_orbit_size==8
    assert design.design_energy==42
    assert sum(count for _,_,count in design.remaining_high)==16
    assert len(_transverse_orbit(design.target_path))==8
    assert _hh_mask(design.sequence,design.target_path)==design.designated_contact_mask
    assert _d2(design.target_path,*design.target_distance_pair)==design.target_squared_distance
    for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(5),Fraction(10)):
        assert 0<design.lower_orbit_probability(q)<=1
        assert design.lower_orbit_probability(q,radius=1) <= design.lower_orbit_probability(q,radius=0)
        assert design.lower_orbit_probability(q)<=design.lower_orbit_probability(Fraction(10))


def test_small_exact_adversarial_energy_and_probability_oracle(small_proof):
    _,design=small_proof
    s=design.sequence
    walk_data=[]
    for walk in _all_cubic_walks(len(s)):
        mask=_hh_mask(s,walk)
        b=int(_d2(walk,*design.target_distance_pair)==design.target_squared_distance)
        walk_data.append((mask,b))
    assert len(walk_data)==13_565
    max_nom=design.design_energy
    chosen=(design.designated_contact_mask,1)
    assert sum(1 for row in walk_data if row==chosen)==8
    for mask,bonus in walk_data:
        energy=8*(mask&design.designated_contact_mask).bit_count()+4*(mask&~design.designated_contact_mask).bit_count()+2*bonus
        if (mask,bonus)==chosen:
            assert energy==max_nom
        else:
            assert energy<=max_nom-2
    # Per-pair errors are globally shared for each contact label, not
    # chosen independently per conformation. Exhaustive random tests are
    # secondary to the exact analytic safety argument.
    pairs=cubic_candidate_pairs(len(s))
    rng=random.Random(1611)
    for _ in range(25):
        eta={bit:rng.choice((-1,0,1)) for bit,(i,j) in enumerate(pairs)
             if s[i]==s[j]=='H'}
        def energy(mask,bonus):
            base=8*(mask&design.designated_contact_mask).bit_count()+4*(mask&~design.designated_contact_mask).bit_count()+2*bonus
            return base+sum(delta for bit,delta in eta.items() if mask&(1<<bit))
        scores=Counter((energy(mask,b), (mask,b)) for mask,b in walk_data)
        maximum=max(e for e,_ in scores)
        winning=set((mask,b) for e,(mask,b) in scores if e==maximum)
        assert winning=={chosen}
        assert scores[maximum,chosen]==8
        for q in (Fraction(2),Fraction(3)):
            partition=sum((q**energy(mask,b) for mask,b in walk_data),Fraction(0))
            p=8*q**maximum/partition
            assert p>=design.lower_orbit_probability(q,radius=1)


def test_untrusted_cut_tree_and_header_tampering_fail_closed(small_proof,compilers,tmp_path):
    files,design=small_proof
    valid=files[1]
    checked=verify_ground_map_cut(valid,compilers['cubic_ground_map_checker'])
    assert checked.high_count==24
    assert sum(n for _,n in checked.multiplicities)==24
    assert len(checked.multiplicities)==3
    truncated=tmp_path/'truncated.cut'
    truncated.write_bytes(valid.read_bytes()[:-1])
    with pytest.raises(ValueError):
        verify_ground_map_cut(truncated,compilers['cubic_ground_map_checker'])
    malicious=tmp_path/'badcut.cut'
    data=bytearray(valid.read_bytes())
    data[11+8+32] |= 1 # false root cut, cannot justify K<5
    malicious.write_bytes(data)
    with pytest.raises(ValueError,match='rejected'):
        verify_ground_map_cut(malicious,compilers['cubic_ground_map_checker'])
    wrong=CheckedNativeGap('H'*8,5,23,2,5**6)
    with pytest.raises(ValueError):
        design_robust_unique_orbit(wrong,checked)
    forged=type(checked)(checked.sequence,checked.threshold,checked.cutoff_sha256,
                         checked.high_count,checked.upper_bins,
                         checked.entries[:-1]+((checked.entries[0][0],((0,0,0),)*8),))
    with pytest.raises(ValueError):
        design_robust_unique_orbit(CheckedNativeGap('H'*8,5,24,2,5**6),forged)


def test_sharp_uncertainty_radius_two_has_explicit_global_tie(small_proof,compilers):
    """Even the optimal robust interval is derived, with an attacking vector."""
    files,design=small_proof
    excited=((0,0,0),(1,0,0),(2,0,0),(3,0,0),
             (3,1,0),(2,1,0),(1,1,0),(0,1,0))
    gap=certify_native_gap(files,compilers['cubic_cut_checker'],excited)
    maps=verify_ground_map_cut(files[1],compilers['cubic_ground_map_checker'])
    attack=prove_sharp_contact_error_radius(design,maps)
    assert attack.critical_error_radius==Fraction(2)
    assert attack.missing_target_pairs>=1
    assert attack.selected_energy_at_boundary==attack.competitor_energy_at_boundary
    assert attack.high_sector_max_checked==attack.selected_energy_at_boundary
    assert attack.global_non_ground_energy_margin==12
    assert sum(value==-2 for _,value in attack.adversarial_pair_errors)==attack.missing_target_pairs
    assert sum(value==+2 for _,value in attack.adversarial_pair_errors)==attack.missing_target_pairs
    # Independently enumerate ALL cubic n8 geometries and verify the attack
    # creates a real GLOBAL tie but no more energetic competitor.
    perturb=dict(attack.adversarial_pair_errors)
    selected_mask=design.designated_contact_mask
    def energy(mask,bonus):
        return (8*(mask&selected_mask).bit_count()
                +4*(mask&~selected_mask).bit_count()
                +2*bonus+sum(eta for bit,eta in perturb.items() if mask&(1<<bit)))
    values=Counter()
    for walk in _all_cubic_walks(8):
        mask=_hh_mask('H'*8,walk)
        b=int(_d2(walk,*design.target_distance_pair)==design.target_squared_distance)
        values[(energy(mask,b),mask,b)]+=1
    global_energy=max(e for e,_,_ in values)
    assert global_energy==attack.selected_energy_at_boundary
    assert values[(global_energy,selected_mask,1)]==8
    assert values[(global_energy,attack.competitor_contact_mask,1)]>0
    # Radius strictly beyond 2 really invalidates the target uniqueness:
    # increasing each extremal attack from magnitude 2 to 3 makes that
    # competitor strictly better than the chosen orbit.
    stronger={bit:3*(eta//2) for bit,eta in perturb.items()}
    def perturbed(mask,b):
        return (8*(mask&selected_mask).bit_count()
                +4*(mask&~selected_mask).bit_count()
                +2*b+sum(v for bit,v in stronger.items() if mask&(1<<bit)))
    assert perturbed(attack.competitor_contact_mask,1)>perturbed(selected_mask,1)


def test_invalid_errors_and_bad_restraints_fail_closed(small_proof):
    _,design=small_proof
    for q in (2,2.0,Fraction(0),Fraction(1,2),True):
        with pytest.raises(ValueError):
            design.lower_orbit_probability(q)
    for r in (-1,2,0.5,True):
        with pytest.raises(ValueError):
            design.lower_orbit_probability(Fraction(2),radius=r)
    for t in (Fraction(0),Fraction(1),0.9,True):
        with pytest.raises(ValueError):
            design.minimum_integer_q(t)
    assert design.minimum_integer_q(Fraction(9,10))>=1


@pytest.mark.skipif(os.getenv('FOLDBOUNDS_FULL_ROBUST_16')!='1',
                    reason='set FOLDBOUNDS_FULL_ROBUST_16=1 for full three-cut 16-mer independent replay')
def test_16mer_verified_unique_design_and_robust_occupancy(compilers,tmp_path):
    producer=compilers['cubic_cut_producer']
    files=tuple(_cut(producer,tmp_path/f'n16_{t}.cut',SEQUENCE,t)
                for t in (10,11,12))
    selected=validate_complete_robust_design(
        files,compilers['cubic_cut_checker'],compilers['cubic_ground_map_checker'],
        EXCITED_WITNESS,max_seconds=1000,
        selected_map=15_668_042_105_480,distance_pair=(0,12),squared_distance=2,
    )
    assert (selected.ground_conformations,selected.ground_contact_count,selected.gap)==(2032,11,2)
    assert selected.target_orbit_size==8
    assert len(selected.remaining_high)==14
    assert sum(n for _,_,n in selected.remaining_high)==2024
    assert selected.design_energy==90
    assert selected.minimum_integer_q(Fraction(9,10),radius=0)==5
    assert selected.minimum_integer_q(Fraction(9,10),radius=1)==6
    assert selected.minimum_integer_q(Fraction(99,100),radius=1)==11
    assert selected.lower_orbit_probability(Fraction(6),radius=1)>Fraction(96,100)
    checked=verify_ground_map_cut(files[1],compilers['cubic_ground_map_checker'],max_seconds=600)
    sharp=prove_sharp_contact_error_radius(selected,checked)
    assert sharp.critical_error_radius==2
    assert sharp.missing_target_pairs==2
    assert sharp.selected_energy_at_boundary==sharp.competitor_energy_at_boundary==86
    assert sharp.high_sector_max_checked==86
    assert sharp.global_non_ground_energy_margin==12
