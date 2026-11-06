"""Independent Python bridge for the newly Lean-defined exact 3D geometry.

Tests deliberately recompute explicit six-direction nearest-neighbor
geometry and index parity without importing the Lean specification
or trusting a model-supplied contact table. The Lean proof is checked
separately by the pinned CI kernel workflow.
"""
from itertools import product

from foldbounds.cubic_tiling import OFFSETS, ROOT
from foldbounds.local_rewire_confidence import (
    StructuralModel, find_rewire_patches, hp_pair_weights, make_macro_target
)
from fractions import Fraction


def integer_neighbor(p, q):
    return sum(abs(a-b) for a,b in zip(p,q)) == 1


def cubic_color(point):
    return sum(point) % 2


def enumerate_rooted_3d(n):
    out = [ROOT]
    for _ in range(n-2):
        updated=[]
        for path in out:
            for dx,dy,dz in OFFSETS:
                nxt=(path[-1][0]+dx,path[-1][1]+dy,path[-1][2]+dz)
                if nxt not in path:
                    updated.append(path+(nxt,))
        out=updated
    return out


def test_six_lean_unit_directions_match_independent_l1_geometry():
    lean_directions=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
    assert set(lean_directions)==set(OFFSETS)
    for p in product(range(-2,3),repeat=3):
        expected={
            (p[0]+d[0],p[1]+d[1],p[2]+d[2])
            for d in lean_directions
        }
        observed={
            q for q in product(range(-3,4),repeat=3)
            if integer_neighbor(p,q)
        }
        assert observed==expected
        assert all(cubic_color(p)!=cubic_color(q) for q in expected)
        assert all(cubic_color(q)==(cubic_color(p)+1)%2 for q in expected)


def test_all_short_canonical_walks_have_exact_index_spatial_parity():
    for n in range(2,9):
        for path in enumerate_rooted_3d(n):
            assert path[:2]==ROOT
            assert len(path)==len(set(path))
            assert all(cubic_color(p)==i%2 for i,p in enumerate(path))
            assert all(integer_neighbor(a,b) for a,b in zip(path,path[1:]))


def test_all_short_saw_contact_edges_cross_index_parity():
    for n in range(2,9):
        for path in enumerate_rooted_3d(n):
            contact_edges=[
                (i,j) for i in range(n) for j in range(i+2,n)
                if integer_neighbor(path[i],path[j])
            ]
            assert all(i%2!=j%2 for i,j in contact_edges)


def test_full_short_hp_parity_obstruction_for_every_sequence():
    for n in range(2,8):
        paths=enumerate_rooted_3d(n)
        for residues in product("HP",repeat=n):
            if len({i%2 for i,ch in enumerate(residues) if ch=="H"})>1:
                continue
            for path in paths:
                assert all(
                    not integer_neighbor(path[i],path[j])
                    for i in range(n) for j in range(i+2,n)
                    if residues[i]==residues[j]=="H"
                )


def test_each_positive_word_is_injective_and_a_valid_sa_witness():
    positive=((1,0,0),(0,1,0),(0,0,1))
    for k in range(0,8):
        all_paths=set()
        for word in product(range(3),repeat=k):
            coords=list(ROOT)
            for digit in word:
                d=positive[digit]
                coords.append(tuple(a+b for a,b in zip(coords[-1],d)))
            walk=tuple(coords)
            assert len(walk)==k+2
            assert len(set(walk))==k+2
            assert all(integer_neighbor(a,b) for a,b in zip(walk,walk[1:]))
            assert all(sum(a)<sum(b) for a,b in zip(walk[1:],walk[2:]))
            all_paths.add(walk)
        assert len(all_paths)==3**k


def test_eight_site_local_rewire_splice_retains_global_occupancy_and_uniq():
    """Separate Python witness of the premises used in CubicSplice.lean."""
    path=make_macro_target(3,1,1)
    n=len(path)
    assert n==24
    model=StructuralModel("H"*n,path,hp_pair_weights(Fraction(3)))
    patches=find_rewire_patches(model)
    assert patches
    for patch in patches:
        assert patch.choices>=2
        original=path[patch.start:patch.start+8]
        assert len(original)==8
        assert len(set(original))==8
        complete=set()
        for variant in patch.variants:
            assert len(variant)==8
            assert sorted(variant)==sorted(original)
            whole=path[:patch.start]+variant+path[patch.start+8:]
            assert len(whole)==n
            assert set(whole)==set(path)
            assert len(set(whole))==n
            assert all(integer_neighbor(a,b) for a,b in zip(whole,whole[1:]))
            complete.add(whole)
        assert len(complete)==patch.choices
