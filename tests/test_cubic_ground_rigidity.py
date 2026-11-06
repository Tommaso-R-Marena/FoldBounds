"""Independent small Hamiltonian oracles and adversarial geometric tests."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
from fractions import Fraction
from itertools import product
from pathlib import Path
import shutil
import subprocess

import pytest

from foldbounds.filled_box_contact_rigidity import (
    FilledBoxCertificateError, certify_filled_box_contacts,
    extract_nonbackbone_contacts, verify_filled_box_certificate,
    contact_map_confidence_upper, split_proper_and_mirror_orbits,
    certify_pairwise_contact_energy, verify_pairwise_contact_energy,
)
from foldbounds.cubic_ground_thermo import (
    CUBE_GROUND_COUNTS, CUBIC_FREE_BLOCK_COUNTS,
    optimal_cubic_walk_cap, verify_cubic_walk_cap,
    cubic_ground_population, verify_cubic_ground_population,
)

ROOT = Path(__file__).resolve().parents[1]


def independent_hamiltonian_paths(dim, m):
    """No production helper: brute visit each occupied grid coordinate once."""
    space = set(product(range(m), repeat=dim))
    deltas = []
    for axis in range(dim):
        for sign in (-1, 1):
            vec = [0]*dim
            vec[axis] = sign
            deltas.append(tuple(vec))

    def recurse(path, seen):
        if len(path) == len(space):
            root = path[0]
            yield tuple(tuple(v[k] - root[k] for k in range(dim)) for v in path)
            return
        current = path[-1]
        for dv in deltas:
            next_ = tuple(current[k]+dv[k] for k in range(dim))
            if next_ in space and next_ not in seen:
                yield from recurse(path+[next_], seen|{next_})

    for p in sorted(space):
        next_ = (p[0]+1,) + p[1:]
        if next_ in space:
            yield from recurse([p, next_], {p, next_})


def independent_contacts(path):
    pairs = []
    for j in range(len(path)):
        for i in range(j - 1):
            if sum(abs(a - b) for a, b in zip(path[i], path[j])) == 1:
                pairs.append((i, j))
    return tuple(sorted(pairs))


@pytest.mark.parametrize("dim,side,walks", [(2,2,2), (2,3,10), (3,2,24)])
def test_complete_small_models_reconstruct_all_maps(dim, side, walks):
    histogram = Counter()
    for walk in independent_hamiltonian_paths(dim, side):
        contacts = independent_contacts(walk)
        histogram[contacts] += 1
        assert contacts == extract_nonbackbone_contacts(walk)
        cert = certify_filled_box_contacts(side, dim, contacts)
        assert verify_filled_box_certificate(cert)
        assert walk in cert.orbit
        assert all(independent_contacts(x) == contacts for x in cert.orbit)
        assert len(cert.orbit) == (2 if dim == 2 else 8)
    assert sum(histogram.values()) == walks
    assert set(histogram.values()) == ({2} if dim==2 else {8})
    assert len(histogram) == walks // (2 if dim == 2 else 8)


def test_malformed_contact_graphs_fail_closed():
    original = next(independent_hamiltonian_paths(3,2))
    base = independent_contacts(original)
    assert base
    invalid = (
        base+base[:1], base[1:], base+((1, 3),),
        ((0, 1),)+base, ((0, 2),)+base, ((0, 10),)+base,
        ((True,3),)+base,
    )
    for candidate in invalid:
        with pytest.raises(FilledBoxCertificateError):
            certify_filled_box_contacts(2,3,candidate)
    with pytest.raises(FilledBoxCertificateError):
        certify_filled_box_contacts(2,True,base)
    with pytest.raises(FilledBoxCertificateError):
        certify_filled_box_contacts(3,3,base)
    with pytest.raises(FilledBoxCertificateError):
        certify_filled_box_contacts(2,3,"not a contact list")


def test_tampered_certificate_and_paths_rejected():
    walk=next(independent_hamiltonian_paths(3,2))
    cert=certify_filled_box_contacts(2,3,independent_contacts(walk))
    assert not verify_filled_box_certificate(replace(cert,orbit=cert.orbit[:4]))
    assert not verify_filled_box_certificate(replace(cert,representative=cert.orbit[-1]))
    assert not verify_filled_box_certificate(replace(cert,nonbackbone_contacts=()))
    assert not verify_filled_box_certificate("forged certificate")
    bad=list(walk); bad[2]=bad[1]
    with pytest.raises(FilledBoxCertificateError):
        extract_nonbackbone_contacts(tuple(bad))
    bad=list(walk); bad[0]=tuple(float(a) for a in bad[0])
    with pytest.raises(FilledBoxCertificateError):
        extract_nonbackbone_contacts(tuple(bad))


def test_exact_short_block_optimality_and_tamper():
    def brute(r):
        if r==0:return 1
        return min(CUBIC_FREE_BLOCK_COUNTS[k] * brute(r-k)
                   for k in range(1,min(r,len(CUBIC_FREE_BLOCK_COUNTS)-1)+1))
    for r in range(0,13):
        cert=optimal_cubic_walk_cap(r)
        assert verify_cubic_walk_cap(cert)
        assert cert.bound == brute(r)
    cert=optimal_cubic_walk_cap(25)
    assert cert.bound == 137795893141213125
    assert cert.bound < 5**25
    assert verify_cubic_walk_cap(cert)
    assert not verify_cubic_walk_cap(replace(cert,bound=cert.bound-1))
    assert not verify_cubic_walk_cap(replace(cert,blocks=cert.blocks[:-1]))
    assert not verify_cubic_walk_cap(replace(cert,extra_steps=24))
    assert not verify_cubic_walk_cap(replace(cert,blocks=(True,)+cert.blocks))


def test_exact_population_bounds_64_different_conditions():
    nominal=cubic_ground_population(3,10,16)
    robust=cubic_ground_population(3,10,16,
                                   correction_lower=Fraction(1,2),correction_upper=2)
    assert nominal.ground_paths==826768
    assert nominal.ground_contact_count==28
    assert nominal.ground_contact_map_count==103346
    assert nominal.total_walk_upper.bound==137795893141213125
    assert nominal.event_probability_lower > Fraction(98,100)
    assert robust.event_probability_lower > Fraction(93,100)
    assert nominal.one_ground_contact_map_probability_upper==Fraction(1,103346)
    assert nominal.single_fold_probability_upper==Fraction(1,826768)
    assert robust.one_ground_contact_map_probability_upper < Fraction(4,100000)
    assert verify_cubic_ground_population(nominal)
    assert verify_cubic_ground_population(robust)
    assert not verify_cubic_ground_population(replace(robust,ground_paths=999))
    assert not verify_cubic_ground_population(replace(robust,event_probability_lower=1))


def test_continuum_threshold_and_rational_gates():
    for q in (1,Fraction(11,10),3,10,11,100,Fraction(1001,100)):
        b=cubic_ground_population(3,q,16,correction_lower=Fraction(1,2),
                                  correction_upper=2)
        assert verify_cubic_ground_population(b)
    for q in (Fraction(0),Fraction(1,2),1.1,True):
        with pytest.raises(ValueError):
            cubic_ground_population(3,q,16)
    for t in (True,-1,0,29,35):
        with pytest.raises(ValueError):
            cubic_ground_population(3,10,t)
    for lower,upper in ((0,2),(3,2),(Fraction(1,2),0),(1,1.5)):
        with pytest.raises(ValueError):
            cubic_ground_population(3,10,16,correction_lower=lower,correction_upper=upper)
    a=cubic_ground_population(3,10,16).event_probability_lower
    b=cubic_ground_population(3,11,16).event_probability_lower
    assert b > a


def test_adversarial_map_mass_upper_is_sharp_for_ground_states_only():
    C=24;g=8
    assert contact_map_confidence_upper(g,C)==Fraction(1,3)
    assert contact_map_confidence_upper(g,C,correction_lower=Fraction(1,2),
                                        correction_upper=2)==Fraction(2,3)
    assert contact_map_confidence_upper(1,C)==Fraction(1,24)
    for x in (0,-1,True,25):
        with pytest.raises(ValueError):
            contact_map_confidence_upper(x,C)
    with pytest.raises(ValueError):
        contact_map_confidence_upper(g,C,correction_lower=0)


def test_independent_all_cubic_8site_gibbs_exhaustive():
    """All full-space anchored n=8 SAWs: event probability encloses exact."""
    seq_len=8
    path=((0,0,0),(1,0,0))
    deltas=list(product((-1,0,1),repeat=3))
    ds=[v for v in deltas if sum(abs(c) for c in v)==1]
    spectrum=Counter()
    def visit(path):
        if len(path)==seq_len:
            spectrum[len(independent_contacts(path))]+=1
            return
        for d in ds:
            v=tuple(path[-1][i]+d[i] for i in range(3))
            if v not in path:visit(path+(v,))
    visit(path)
    assert sum(spectrum.values())==13565
    assert spectrum[5]==24
    for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(10)):
        Z=sum((count*q**score for score,count in spectrum.items()),Fraction(0))
        for t in (1,2,3,4,5):
            event=sum((count*q**score for score,count in spectrum.items()
                       if score>=t),Fraction(0))
            result=cubic_ground_population(2,q,t)
            assert result.event_probability_lower <= event/Z


def _compile(tmp_path, name):
    if shutil.which('g++') is None:
        pytest.skip('g++ not available for independent C++ audit')
    source=ROOT/'experiments'/name
    executable=tmp_path/name.replace('.cpp','')
    subprocess.run(['g++','-O3','-std=c++17','-Wall','-Wextra','-Werror',
                    str(source),'-o',str(executable)],check=True,
                   stdout=subprocess.PIPE,stderr=subprocess.PIPE,timeout=60)
    return executable


def test_independent_cpp_cubic_ground_count(tmp_path):
    executable=_compile(tmp_path,'cubic_count_independent.cpp')
    outputs={}
    for side in (2,3):
        out=subprocess.run([str(executable),str(side)],capture_output=True,
                           text=True,check=True,timeout=40)
        outputs[side]=int(out.stdout.strip())
    assert outputs==CUBE_GROUND_COUNTS


def test_independent_cpp_all_contact_masks_and_reconstruction(tmp_path):
    executable=_compile(tmp_path,'cubic_ground_contact_map_audit.cpp')
    for side,expected_map_count in ((2,3),(3,103346)):
        out=subprocess.run([str(executable),str(side)],capture_output=True,
                           text=True,check=True,timeout=45).stdout.splitlines()
        assert f'first_E_paths={CUBE_GROUND_COUNTS[side]}' in out[0]
        assert f'distinct_contact_masks={expected_map_count}' in out[0]
        assert 'multiplicity,count_of_masks'==out[1]
        assert f'8,{expected_map_count}' == out[2]
        samples=[tuple(map(int,x.split()[1:])) for x in out[3:] if x.startswith('SAMPLE')]
        assert len(samples)>=8
        for sample in samples:
            cells=tuple((v%side,(v//side)%side,v//(side*side)) for v in sample)
            root=cells[0]
            anchored=tuple(tuple(c[i]-root[i] for i in range(3)) for c in cells)
            assert anchored[0]==(0,0,0) and anchored[1]==(1,0,0)
            cert=certify_filled_box_contacts(side,3,independent_contacts(anchored))
            assert anchored in cert.orbit
            assert verify_filled_box_certificate(cert)


def test_independent_cpp_short_block_multiplicities(tmp_path):
    executable=_compile(tmp_path,'cubic_short_block_oracle.cpp')
    out=subprocess.run([str(executable),'10'],capture_output=True,
                       text=True,check=True,timeout=40).stdout.splitlines()
    values=tuple(int(row.split(',')[1]) for row in out)
    assert values == CUBIC_FREE_BLOCK_COUNTS


def test_distance_only_models_cannot_select_chirality():
    """Two proper-orientation classes have exactly the same pair energy."""
    for side in (2,):
        walk = next(independent_hamiltonian_paths(3, side))
        cert = certify_filled_box_contacts(side, 3, independent_contacts(walk))
        right, left = split_proper_and_mirror_orbits(cert)
        assert len(right)==len(left)==4
        assert set(right).isdisjoint(left)
        weights = {pair: ((i*31+j*7)%13)-6
                   for i in range(len(walk))
                   for j in range(i+2,len(walk))
                   for pair in [(i,j)]}
        def energy(path):
            contacts=independent_contacts(path)
            return sum(weights[p] for p in contacts)
        assert len({energy(path) for path in right+left})==1
        assert all(independent_contacts(p)==cert.nonbackbone_contacts
                   for p in right+left)
        assert contact_map_confidence_upper(4,8)==Fraction(1,2)
        assert contact_map_confidence_upper(4,8,correction_lower=Fraction(1,2),
                                            correction_upper=2)==Fraction(4,5)
        assert not verify_filled_box_certificate(replace(cert,orbit=right))
        with pytest.raises(FilledBoxCertificateError):
            split_proper_and_mirror_orbits(replace(cert,orbit=right))


def test_large_full_box_inverse_proves_no_exhaustive_search_needed():
    """Full 4³ and 8² boxes reconstructed using only their labeled graphs."""
    for side, dimension in ((4, 3), (8, 2)):
        base=[]
        for y in range(side):
            for x in (range(side) if y % 2 == 0 else range(side - 1,-1,-1)):
                base.append((x,y))
        if dimension == 2:
            path=tuple(base)
        else:
            path=tuple((x,y,z)
                       for z in range(side)
                       for (x,y) in (base if z % 2 == 0 else base[::-1]))
        contacts=independent_contacts(path)
        assert len(contacts)==(
            (side-1)**2 if dimension==2 else 2*side**3-3*side**2+1)
        cert=certify_filled_box_contacts(side,dimension,contacts)
        assert verify_filled_box_certificate(cert)
        assert path in cert.orbit
        assert len(cert.orbit)==(2 if dimension==2 else 8)
        assert all(independent_contacts(other)==contacts for other in cert.orbit)
        if dimension==3:
            assert tuple(map(len,split_proper_and_mirror_orbits(cert)))==(4,4)


def test_arbitrary_rational_pair_potentials_have_exact_orbit_energy():
    for dim,side in ((2,3),(3,2),(3,3)):
        if (dim,side)==(3,3):
            base = [(x,y) for y in range(3)
                    for x in (range(3) if y%2==0 else range(2,-1,-1))]
            walk=tuple((x,y,z) for z in range(3)
                       for (x,y) in (base if z%2==0 else base[::-1]))
        else:
            walk=next(independent_hamiltonian_paths(dim,side))
        cert=certify_filled_box_contacts(side,dim,independent_contacts(walk))
        n=len(walk)
        energies={(i,j):Fraction((i*11+j*7)%13-6, (i+j)%3+1)
                  for i in range(n) for j in range(i+2,n)}
        witness=certify_pairwise_contact_energy(cert,energies)
        assert verify_pairwise_contact_energy(cert,energies,witness)
        assert witness.orbit_size==(2 if dim==2 else 8)
        assert witness.energy == sum((energies[p] for p in cert.nonbackbone_contacts),Fraction(0))
        assert not verify_pairwise_contact_energy(cert,energies,replace(witness,energy=0)) or witness.energy==0
        with pytest.raises(FilledBoxCertificateError):
            certify_pairwise_contact_energy(cert,{(True,3):1})
        with pytest.raises(FilledBoxCertificateError):
            certify_pairwise_contact_energy(cert,{(0,3):1.25})
        with pytest.raises(FilledBoxCertificateError):
            certify_pairwise_contact_energy(cert,{(0,1):1})
