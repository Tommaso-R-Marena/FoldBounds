"""Independent checker, malicious-cut rejection and nonenumerative 3D tail proof."""
from __future__ import annotations

from fractions import Fraction
import os
from pathlib import Path
import shutil
import subprocess

import pytest

from foldbounds.cubic_contact_atlas import build_cubic_atlas
from foldbounds.cubic_cut_certificate import (
    certify_native_gap, contact_score, inspect_cut_header,
    verify_cut_certificate, CheckedNativeGap,
)

SOURCE = Path(__file__).resolve().parents[1]/"experiments"
SEQUENCE = "HHHPPHHHHHHHPPHH"
GROUND_WITNESS = (
    (0,0,0),(1,0,0),(2,0,0),(3,0,0),
    (3,1,0),(2,1,0),(1,1,0),(0,1,0),
    (0,1,1),(0,0,1),(1,0,1),(2,0,1),
    (3,0,1),(3,1,1),(2,1,1),(1,1,1),
)
EXCITED_WITNESS = (
    (0,0,0),(1,0,0),(2,0,0),(3,0,0),
    (3,1,0),(2,1,0),(1,1,0),(0,1,0),
    (0,1,1),(1,1,1),(2,1,1),(3,1,1),
    (3,0,1),(2,0,1),(1,0,1),(0,0,1),
)


@pytest.fixture(scope="module")
def executables(tmp_path_factory):
    compiler = shutil.which("g++")
    if not compiler:
        pytest.skip("a C++17 compiler is required for independent cut-tree checking")
    directory = tmp_path_factory.mktemp("cubic_bnb_exec")
    files = []
    for name in ("cubic_cut_producer", "cubic_cut_checker"):
        binary = directory/name
        subprocess.run(
            [compiler,"-std=c++17","-O3","-Wall","-Wextra","-Werror",
             str(SOURCE/f"{name}.cpp"),"-o",str(binary)],
            check=True,timeout=60,capture_output=True,text=True,
        )
        files.append(binary)
    return tuple(files)


def generate(producer, path, sequence, threshold, mode="geometry"):
    subprocess.run(
        [str(producer),sequence,str(threshold),str(path),mode],
        check=True,timeout=90,capture_output=True,text=True,
    )
    return path


def test_independent_exact_3d_histogram_small(executables,tmp_path):
    producer,checker = executables
    for n in (4,6,8):
        atlas = build_cubic_atlas(n)
        for sequence in ("H"*n,"P"*n,"HP"*(n//2),"H"*(n//2)+"P"*(n//2)):
            spectrum=atlas.spectrum(sequence)
            for threshold in range(1,min(5,max(spectrum)+2)):
                file=tmp_path/f"small_{n}_{sequence}_{threshold}.cut"
                generate(producer,file,sequence,threshold)
                checked=verify_cut_certificate(file,checker)
                assert checked.qualifying_conformations == sum(
                    count for k,count in spectrum.items() if k>=threshold
                )


def test_12mer_unique_orbit_tail_matches_independently_audited_atlas(executables,tmp_path):
    producer,checker=executables
    atlas=build_cubic_atlas(12)
    for sequence in ("HHHPHPHHPHPH","HHPHPHHPHPHH","HPHPHHPHPHHH"):
        for threshold in (7,8):
            file=tmp_path/f"n12_{sequence}_{threshold}.cut"
            generate(producer,file,sequence,threshold)
            result=verify_cut_certificate(file,checker)
            assert result.qualifying_conformations==sum(
                count for k,count in atlas.spectrum(sequence).items()
                if k>=threshold
            )
            assert result.qualifying_conformations == (8 if threshold==7 else 0)
    geometry=tmp_path/"geo.cut"
    capacity=tmp_path/"cap.cut"
    generate(producer,geometry,"HHHPHPHHPHPH",7)
    generate(producer,capacity,"HHHPHPHHPHPH",7,"capacity")
    geo=verify_cut_certificate(geometry,checker)
    cap=verify_cut_certificate(capacity,checker)
    assert geo.qualifying_conformations == cap.qualifying_conformations == 8
    assert geo.total_nodes < cap.total_nodes


def test_n16_ground_count_replayed_in_hosted_ci(executables,tmp_path):
    producer,checker=executables
    proof_path=tmp_path/"n16_ground.cut"
    generate(producer,proof_path,SEQUENCE,11)
    certified=verify_cut_certificate(proof_path,checker,max_seconds=150)
    assert certified.qualifying_conformations==2032
    assert certified.total_nodes==10_851_705
    assert certified.cuts==8_444_857
    assert certified.terminal_leaves==2032
    assert len(certified.upper_bins)>=1
    assert proof_path.stat().st_size < 1_500_000


@pytest.fixture(scope="module")
def certificates_16(executables,tmp_path_factory):
    if os.getenv("FOLDBOUNDS_FULL_CUBIC_CUT_AUDIT") != "1":
        pytest.skip("run full 3D gap test via FOLDBOUNDS_FULL_CUBIC_CUT_AUDIT=1")
    producer,checker=executables
    directory=tmp_path_factory.mktemp("n16_nonexhaustive_proof")
    files=[]
    for t in (10,11,12):
        path=directory/f"threshold_{t}.cut"
        generate(producer,path,SEQUENCE,t)
        files.append(path)
    return tuple(files)


def test_verified_native_gap_16_and_all_q_occupancy(certificates_16,executables):
    _,checker=executables
    assert contact_score(SEQUENCE,GROUND_WITNESS)==11
    assert contact_score(SEQUENCE,EXCITED_WITNESS)==9
    proof=certify_native_gap(certificates_16,checker,EXCITED_WITNESS,max_seconds=180)
    assert proof.sequence==SEQUENCE
    assert proof.ground_contacts==11
    assert proof.ground_conformations==2032
    assert proof.gap==2
    assert proof.max_total_conformations==5**14
    summaries=tuple(inspect_cut_header(p) for p in certificates_16)
    assert tuple(s.qualifying_conformations for s in summaries)==(2032,2032,0)
    assert summaries[1].total_nodes==10_851_705
    assert summaries[1].cuts==8_444_857
    assert summaries[1].total_nodes < 5**14 // 100
    assert certificates_16[1].stat().st_size < 1_500_000
    assert proof.minimum_integer_weight(Fraction(1,2))==939
    assert proof.minimum_integer_weight(Fraction(9,10))==2816
    assert proof.minimum_integer_weight(Fraction(99,100))==9336
    for q in (2,100,1000,2000,5000,5200,10000):
        coarse=Fraction(2032*q*q,2032*q*q+5**14-2032)
        weighted=2032*q**11
        upper=min([Fraction((5**14-2032)*q**9)]+[
            sum((v*q**k for k,v in envelope),Fraction(0))
            for envelope in proof.upper_non_ground_envelopes
        ])
        assert proof.lower_native_probability(Fraction(q))==weighted/(weighted+upper)
        assert proof.lower_native_probability(Fraction(q)) >= coarse
    assert proof.lower_native_probability(Fraction(2816)) >= Fraction(9,10)
    assert proof.lower_native_probability(Fraction(2815)) < Fraction(9,10)
    assert len(proof.upper_non_ground_envelopes)==2


def test_cut_tampering_and_truncation_are_rejected(executables,tmp_path):
    producer,checker=executables
    valid=generate(producer,tmp_path/"valid.cut","HHHPHPHHPHPH",7)
    assert verify_cut_certificate(valid,checker).qualifying_conformations==8
    original=valid.read_bytes()
    payload_offset=11+len("HHHPHPHHPHPH")+32
    mutated=bytearray(original)
    mutated[payload_offset] |= 1  # claim at root that 7 contacts is impossible
    malicious=tmp_path/"malicious_root.cut"
    malicious.write_bytes(mutated)
    with pytest.raises(ValueError,match="REJECTED"):
        verify_cut_certificate(malicious,checker)
    mutated=bytearray(original)
    high_header_offset=11+len("HHHPHPHHPHPH")+16
    mutated[high_header_offset] ^= 1  # falsify qualifying tail count
    malicious=tmp_path/"malicious_count.cut"
    malicious.write_bytes(mutated)
    with pytest.raises(ValueError,match="REJECTED|invalid count header"):
        verify_cut_certificate(malicious,checker)
    malicious=tmp_path/"truncated.cut"
    malicious.write_bytes(original[:-1])
    with pytest.raises(ValueError,match="length"):
        verify_cut_certificate(malicious,checker)
    malicious=tmp_path/"extra_bytes.cut"
    malicious.write_bytes(original+b"\0")
    with pytest.raises(ValueError,match="length"):
        verify_cut_certificate(malicious,checker)
    with pytest.raises(FileNotFoundError):
        verify_cut_certificate(valid,tmp_path/"nonexistent_checker")


def test_forged_native_witness_and_invalid_q(executables,tmp_path):
    producer,checker=executables
    path=generate(producer,tmp_path/"valid.cut","HHHPHPHHPHPH",7)
    with pytest.raises(ValueError):
        certify_native_gap((path,path,path),checker,EXCITED_WITNESS)
    for bad in (0,1.0,True,Fraction(0)):
        with pytest.raises(ValueError):
            CheckedNativeGap(SEQUENCE,11,2032,2,5**14).lower_native_probability(bad)
    with pytest.raises(ValueError):
        contact_score(SEQUENCE,tuple([(0,0,0)]*16))


def test_small_positive_native_gap_certificate_and_all_q_bound(executables,tmp_path):
    """Fast end-to-end positive gap2 proof, always exercised in hosted CI."""
    producer, checker = executables
    s = "H"*8
    atlas = build_cubic_atlas(8)
    assert atlas.spectrum(s)[5] == 24
    assert 4 not in atlas.spectrum(s)
    excited_path = (
        (0,0,0),(1,0,0),(2,0,0),(3,0,0),
        (3,1,0),(2,1,0),(1,1,0),(0,1,0),
    )
    assert contact_score(s, excited_path) == 3
    files = tuple(generate(producer,tmp_path/f"small_t{t}.cut",s,t)
                  for t in (4,5,6))
    proof = certify_native_gap(files,checker,excited_path)
    assert (proof.ground_contacts,proof.ground_conformations,proof.gap) == (5,24,2)
    assert proof.max_total_conformations == 5**6
    assert proof.upper_non_ground_envelopes[0]
    for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(100)):
        z=atlas.partition(s,q)
        actual = 24*q**5/z
        assert proof.lower_native_probability(q) <= actual
    assert proof.minimum_integer_weight(Fraction(9,10)) >= 1
