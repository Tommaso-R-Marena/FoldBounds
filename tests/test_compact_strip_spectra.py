"""Independent rooted-walk oracles and adversarial checks for compact tiles."""
from __future__ import annotations

from collections import Counter
from dataclasses import replace
from fractions import Fraction
from itertools import product
from math import isqrt
from random import Random

import pytest

from foldbounds.compact_strip_spectra import (
    CompactStripCertificate,
    best_compact_strip_mass,
    build_compact_strip_certificate,
    compact_strip_spectrum,
    compact_plus_directed_mass,
    tile_paths,
    validate_tile_library,
    verify_compact_strip_certificate,
)


def independent_tile_paths(height):
    """Separate bitmask traversal (different representation/order)."""
    width = 4
    all_used = (1 << (height * width)) - 1
    end = 3
    paths = []

    def visit(idx, mask, sequence):
        if mask == all_used:
            if idx == end:
                paths.append(tuple((v % width, v // width) for v in sequence))
            return
        x, y = idx % width, idx // width
        for nx, ny in ((x - 1, y), (x, y - 1), (x + 1, y), (x, y + 1)):
            if 0 <= nx < width and 0 <= ny < height:
                nxt = ny * width + nx
                if mask & (1 << nxt) or (nxt == end and mask.bit_count() + 1 < height * width):
                    continue
                visit(nxt, mask | 1 << nxt, sequence + [nxt])

    visit(0, 1, [0])
    return tuple(paths)


def construct_path(parts, height, remaining):
    out = []
    for block, path in enumerate(parts):
        out.extend((x + 4 * block, y) for x, y in path)
    if not parts:
        return tuple((i, 0) for i in range(remaining))
    out.extend((4 * len(parts) + i, 0) for i in range(remaining))
    return tuple(out)


def naive_score(sequence, path):
    return sum(sequence[i] == sequence[j] == 'H' and j - i > 1 and
               abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1
               for i in range(len(path)) for j in range(i + 2, len(path)))


def brute_tiled_spectrum(sequence, height, disjoint=False):
    n = len(sequence)
    blocks, tail = divmod(n, 4 * height)
    if not blocks:
        return {} if disjoint else {0: 1}
    ps = independent_tile_paths(height)
    counts = Counter()
    for variants in product(ps, repeat=blocks):
        if variants[0][1] != (1, 0):
            continue
        path = construct_path(variants, height, tail)
        assert path[:2] == ((0, 0), (1, 0))
        assert len(path) == n and len(set(path)) == n
        assert all(abs(a[0]-b[0]) + abs(a[1]-b[1]) == 1 for a,b in zip(path,path[1:]))
        if disjoint:
            steps = {(b[0]-a[0], b[1]-a[1]) for a,b in zip(path,path[1:])}
            if not {(-1,0),(0,1),(0,-1)}.issubset(steps):
                continue
        counts[naive_score(sequence,path)] += 1
    return dict(counts)


@pytest.mark.parametrize('height',range(3,9))
def test_tile_library_independently_exact(height):
    actual = tile_paths(height)
    separate = independent_tile_paths(height)
    assert set(actual) == set(separate)
    assert validate_tile_library(height, actual)
    assert len(actual) >= 4
    assert any(p[1] == (1,0) for p in actual)
    assert len(actual) == (4, 8, 23, 55, 144, 360)[height-3]


@pytest.mark.parametrize('height', (3,4,5))
def test_exact_spectrum_independent_enumeration(height):
    n = 8 * height
    rng = Random(height + 20261105)
    sequences = ['H'*n, 'P'*n, ('HPPH'*n)[:n], ('HP'*n)[:n],
                 ''.join(rng.choice('HP') for _ in range(n))]
    for sequence in sequences:
        for tail in (0, 3):
            seq = sequence + 'HPH'[:tail]
            for disjoint in (False, True):
                expected = brute_tiled_spectrum(seq, height, disjoint)
                result = dict(compact_strip_spectrum(
                    seq,height=height,disjoint_from_partially_directed=disjoint))
                assert result == expected
                cert = build_compact_strip_certificate(
                    seq,height=height,disjoint_from_partially_directed=disjoint)
                for q in (Fraction(1),Fraction(3,2),Fraction(5)):
                    t = verify_compact_strip_certificate(cert,q,1)
                    assert t.partition_lower == sum(
                        (c*q**k for k,c in expected.items()), Fraction(0))
                    assert t.event_lower == sum(
                        (c*q**k for k,c in expected.items() if k >= 1), Fraction(0))
                    assert t.complement_lower == t.partition_lower - t.event_lower
                    assert t.witnessed_paths == sum(expected.values())


def all_saws(n):
    paths = [((0,0),(1,0))]
    for _ in range(2,n):
        nxt=[]
        for path in paths:
            x,y=path[-1]
            for nx,ny in ((x+1,y),(x,y+1),(x-1,y),(x,y-1)):
                if (nx,ny) not in path:
                    nxt.append(path+((nx,ny),))
        paths=nxt
    return paths


def test_3x4_tiling_is_true_subset_of_full_saw_ensemble():
    n=12
    paths = all_saws(n)
    for seq in ('H'*n, ('HP'*n)[:n], 'P'*n):
        full = Counter(naive_score(seq,p) for p in paths)
        tile = dict(compact_strip_spectrum(seq,height=3))
        for k,c in tile.items():
            assert c <= full[k]
        for q in (Fraction(1),Fraction(2),Fraction(7,2)):
            z=sum((c*q**k for k,c in full.items()),Fraction(0))
            b=verify_compact_strip_certificate(build_compact_strip_certificate(seq,height=3),q)
            assert 0 < b.partition_lower <= z


def test_all_H_compact_maximal_contact_64():
    n=64
    cert = build_compact_strip_certificate('H'*n,height=8)
    assert cert.spectrum == ((49,52920),)
    assert build_compact_strip_certificate(
        'H'*n,height=8,disjoint_from_partially_directed=True).spectrum == cert.spectrum
    max_contacts = n + 1 - isqrt(4*n)
    assert 49 == max_contacts
    b=verify_compact_strip_certificate(cert,Fraction(10))
    assert b.partition_lower == 52920 * 10**49
    reflected=verify_compact_strip_certificate(build_compact_strip_certificate(
        'H'*n,height=8,include_reflection=True),Fraction(10))
    assert reflected.witnessed_paths == 105840
    assert reflected.partition_lower == 2*b.partition_lower
    assert b.event_lower == b.partition_lower


def test_many_tiles_non_enumerative_and_allH_closed_form():
    h=8
    for blocks in (1,2,3,4,12):
        n=4*h*blocks
        c=build_compact_strip_certificate('H'*n,height=h)
        assert c.spectrum == (((h-1)*(4*blocks-1),147*360**(blocks-1)),)
        d=build_compact_strip_certificate('H'*n,height=h,
                                          disjoint_from_partially_directed=True)
        assert sum(x[1] for x in d.spectrum) <= c.spectrum[0][1]
        assert d.spectrum[0][0] == c.spectrum[0][0]


def test_not_double_counted_by_partial_directed_families():
    from foldbounds.partial_directed import is_partially_directed_path
    for h in (3,4,5):
        ps = independent_tile_paths(h)
        for variants in product(ps,repeat=2):
            if variants[0][1] != (1,0):
                continue
            p=construct_path(variants,h,0)
            used={ (b[0]-a[0],b[1]-a[1]) for a,b in zip(p,p[1:]) }
            if {(-1,0),(0,1),(0,-1)}.issubset(used):
                assert not is_partially_directed_path(p)


def test_hybrid_improves_existing_model_relative_lower_bound():
    from foldbounds.partial_directed import partially_directed_spectrum
    n=64;q=10
    pdsaw=sum(c*q**k for k,c in partially_directed_spectrum('H'*n))
    combined=compact_plus_directed_mass('H'*n,Fraction(q),height=8)
    assert combined.partition_lower > 40*pdsaw
    assert combined.partition_lower == pdsaw + 105840*q**49


def test_invalid_inputs_and_tampering_rejected():
    cert=build_compact_strip_certificate('H'*32,height=8)
    with pytest.raises(ValueError):
        verify_compact_strip_certificate(replace(cert,spectrum=((49,12345),)),Fraction(2))
    with pytest.raises(ValueError):
        verify_compact_strip_certificate(replace(cert,height=3),Fraction(2))
    with pytest.raises(ValueError):
        verify_compact_strip_certificate(replace(cert,disjoint_from_partially_directed=0),Fraction(2))
    for bad in (Fraction(1,2),2.0,True):
        with pytest.raises(ValueError):verify_compact_strip_certificate(cert,bad)
    for h in (False,2,9,5.0):
        with pytest.raises(ValueError):tile_paths(h)
    with pytest.raises(ValueError):compact_strip_spectrum('HHPX'*8)
    with pytest.raises(ValueError):compact_strip_spectrum('HH'*16,disjoint_from_partially_directed=1)
    with pytest.raises(ValueError):validate_tile_library(4,(tile_paths(4)[0],tile_paths(4)[0]))


def test_best_choice_safe_without_addition():
    n=64; q=Fraction(5)
    height, mass=best_compact_strip_mass('H'*n,q)
    assert height in (4,5,6,7,8)
    assert mass.partition_lower>=verify_compact_strip_certificate(
        build_compact_strip_certificate('H'*n,height=8),q).partition_lower


def test_exact_continuum_strong_coupling_certificate():
    from foldbounds.compact_strip_spectra import (
        build_tail_dominance_certificate, verify_tail_dominance_certificate,
    )
    from foldbounds.partial_directed import partially_directed_spectrum
    cert=build_tail_dominance_certificate('H'*64,height=8,q_min=Fraction(10))
    assert verify_tail_dominance_certificate(cert)
    assert cert.ratio_lower > 80
    for q in (Fraction(10),Fraction(13,2)*2,Fraction(20),Fraction(100),Fraction(1000)):
        pdsaw=sum((c*q**k for k,c in partially_directed_spectrum('H'*64)),Fraction(0))
        hybrid=compact_plus_directed_mass('H'*64,q,height=8)
        assert hybrid.partition_lower / pdsaw >= cert.ratio_lower
    with pytest.raises(ValueError):
        verify_tail_dominance_certificate(replace(cert,ratio_lower=cert.ratio_lower+1))
    with pytest.raises(ValueError):
        build_tail_dominance_certificate('H'*16,height=8,q_min=Fraction(10))


def test_mirror_produces_exact_disjoint_heterogeneous_subset():
    # An independent, explicit coordinate reflection with root-edge and
    # contact checks (not just algebraically multiplying a count).
    h=3
    ps=independent_tile_paths(h)
    n=2*h*4+2
    seq=('HPHHPHHP'*4)[:n]
    for disjoint in (False,True):
        original=brute_tiled_spectrum(seq,h,disjoint)
        mirrored=Counter()
        for variants in product(ps,repeat=2):
            if variants[0][1] != (1,0):continue
            path=construct_path(variants,h,2)
            moves={(b[0]-a[0],b[1]-a[1]) for a,b in zip(path,path[1:])}
            if disjoint and not {(-1,0),(0,1),(0,-1)}.issubset(moves):continue
            reflected=tuple((x,-y) for x,y in path)
            assert reflected[:2]==((0,0),(1,0))
            assert len(set(reflected))==n
            assert reflected!=path
            assert naive_score(seq,reflected)==naive_score(seq,path)
            mirrored[naive_score(seq,reflected)]+=1
        assert dict(mirrored)==original
        cert=build_compact_strip_certificate(seq,height=h,
             disjoint_from_partially_directed=disjoint, include_reflection=True)
        assert dict(cert.spectrum)=={k:2*v for k,v in original.items()}
        verify_compact_strip_certificate(cert,Fraction(7,3))
