"""Independent exact-integer geometry regression on four public PDB coordinate models."""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import sys
import random
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"data"))
from experimental_ca_packed import STRUCTURES, decode
from foldbounds.sparse_distance_certificate import (
    make_certificate,reconstruct_exact,verify_certificate,robust_contacts,
    make_signed_certificate,reconstruct_signed_exact,verify_signed_certificate,
    point_digest,select_anchors)

EXPECTED={
    "1A8O":(66,(53,22,6,59),1953,88,78,1858,17,411),
    "1LCD":(51,(26,0,32,15),1128,94,87,1023,18,117),
    "2BEG":(26,(25,11,0,14),253,1,0,251,2,13),
    "7DDO":(597,(117,0,408,590),176121,1235,1123,174769,229,5693),
}
EXPECTED_SOURCE_SHA={
    "1A8O":"f768d71103ff07b07c052871229613a695bcc70d",
    "1LCD":"4c0d82a7b0ae16a7051e296e6b9c7dc003f95255",
    "2BEG":"ac199080295896f566ef6dcdc192a24bc2a6bf70",
    "7DDO":"11279b6a1552b3e3cbccb1ef4469e7207d539d98",
}
EXPECTED_COORD_DIGEST={
    "1A8O":"a6bdc57cb2e3716e5170e847783df5bf7fb7a7f615a4df9903b169996d308c64",
    "1LCD":"f4049fc02bfded77c785910c1443e9036dfe5dc88befb14608c52519a4c84e47",
    "2BEG":"6fd18f622ca805c549bdce574de81a47ccd9598bd769d755266e26813d4c9479",
}

@pytest.mark.parametrize("pdb",tuple(EXPECTED))
def test_public_pdb_models_exact(pdb):
    n,ids,pairs,contacts,yes,no,uncertain,naive=EXPECTED[pdb]
    meta=STRUCTURES[pdb];p=decode(meta["packed"])
    assert meta["sha"]==EXPECTED_SOURCE_SHA[pdb]
    assert len(p)==meta["n"]==n
    if pdb in EXPECTED_COORD_DIGEST:
        assert point_digest(p)==EXPECTED_COORD_DIGEST[pdb]
    cert=make_certificate(p)
    assert cert.anchor_indices==ids
    assert reconstruct_exact(cert)==p
    assert verify_certificate(cert)
    classified=robust_contacts(cert,cutoff_mA=8000,epsilon_mA=10)
    assert (classified.all_pairs,classified.nominal_contacts,
            classified.guaranteed_contacts,classified.guaranteed_noncontacts,
            classified.undecided_pairs)==(pairs,contacts,yes,no,uncertain)
    if pdb!="7DDO":
        base=make_certificate(p,(0,1,2,3))
        assert robust_contacts(base,epsilon_mA=10).undecided_pairs==naive
    else:
        base=make_certificate(p,(0,1,2,3))
        assert robust_contacts(base,epsilon_mA=10).undecided_pairs==naive
    assert classified.guaranteed_contacts<=classified.nominal_contacts
    assert classified.guaranteed_contacts+classified.guaranteed_noncontacts+classified.undecided_pairs==pairs
    sign=make_signed_certificate(p)
    assert reconstruct_signed_exact(sign)==p
    assert verify_signed_certificate(sign)
    assert all(t in (-1,1) for i,t in enumerate(sign.orientation_signs) if i not in sign.anchor_indices)
    assert sum(v is not None for v in sign.squared_to_anchors)==n-3
    assert 3*n-6 < 4*n-10 < n*(n-1)//2

def test_sparse_four_anchor_tampering_fail_closed():
    p=((0,0,0),(2000,0,0),(0,3000,0),(0,0,4000),
       (1000,1000,1000),(8000,8000,8000),(11000,-3000,5000))
    c=make_certificate(p,(0,1,2,3))
    assert reconstruct_exact(c)==p
    row=list(c.squared_to_anchors)
    row[4]=tuple(x+1 for x in row[4])
    assert not verify_certificate(replace(c,squared_to_anchors=tuple(row)))
    assert not verify_certificate(replace(c,recovered_sha256="0"*64))
    assert not verify_certificate(replace(c,anchor_indices=(0,0,1,2)))
    assert not verify_certificate(replace(c,anchors=(p[0],)*4))
    zero=robust_contacts(c,epsilon_mA=0,min_sequence_separation=1)
    assert zero.undecided_pairs==0
    prev=zero
    for eps in (1,10,50,100):
        now=robust_contacts(c,epsilon_mA=eps,min_sequence_separation=1)
        assert now.guaranteed_contacts<=prev.guaranteed_contacts
        assert now.guaranteed_noncontacts<=prev.guaranteed_noncontacts
        prev=now

def test_signed_three_anchor_tampering_and_mirror():
    p=((0,0,0),(2000,0,0),(0,3000,0),(0,0,4000),(1000,1000,1000))
    c=make_signed_certificate(p,(0,1,2))
    assert reconstruct_signed_exact(c)==p
    signs=list(c.orientation_signs)
    signs[3]*=-1
    assert not verify_signed_certificate(replace(c,orientation_signs=tuple(signs)))
    row=list(c.squared_to_anchors)
    row[4]=(row[4][0]+1,row[4][1],row[4][2])
    assert not verify_signed_certificate(replace(c,squared_to_anchors=tuple(row)))
    assert not verify_signed_certificate(replace(c,anchor_indices=(0,0,1)))

def test_random_integer_geometry_oracle():
    rng=random.Random(2701)
    for _ in range(64):
        p=tuple(tuple(rng.randrange(-30000,30000) for _ in range(3)) for _ in range(10))
        four=make_certificate(p)
        three=make_signed_certificate(p)
        assert reconstruct_exact(four)==p
        assert reconstruct_signed_exact(three)==p
        c0=robust_contacts(four,epsilon_mA=0)
        assert c0.undecided_pairs==0

def test_bad_coplanar_input_fails_closed():
    p=((0,0,0),(1,0,0),(2,0,0),(3,0,0))
    with pytest.raises(ValueError):select_anchors(p)


@pytest.mark.parametrize("pdb",tuple(EXPECTED))
def test_exact_jacobian_full_rank_and_local_constraint_minimality(pdb):
    """Each 3-distance site map has exact Jacobian determinant 8*V != 0."""
    p=decode(STRUCTURES[pdb]["packed"])
    cert=make_signed_certificate(p)
    a0,a1,a2=cert.anchors
    def sub(x,y):return tuple(v-w for v,w in zip(x,y))
    def dot(x,y):return sum(v*w for v,w in zip(x,y))
    def cross(x,y):
        return (x[1]*y[2]-x[2]*y[1],
                x[2]*y[0]-x[0]*y[2],
                x[0]*y[1]-x[1]*y[0])
    u,v=sub(a1,a0),sub(a2,a0)
    normal=cross(u,v)
    assert dot(normal,normal)>0
    for i,x in enumerate(p):
        if i in cert.anchor_indices:continue
        jacobian_rows=[tuple(2*z for z in sub(x,a)) for a in (a0,a1,a2)]
        detj=dot(jacobian_rows[0],cross(jacobian_rows[1],jacobian_rows[2]))
        signed_volume=dot(normal,sub(x,a0))
        assert detj==8*signed_volume
        assert detj!=0
