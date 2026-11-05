"""3D cubic HP geometry, independent full-walk audits and exact all-q designs."""
from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
import shutil
import subprocess

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.cubic_contact_atlas import (
    CubicContactAtlas, audit_cubic_atlas, audit_cubic_landscape,
    build_cubic_atlas, build_cubic_landscape,
)
from foldbounds.exact_phase import certify_sequence_comparison
from foldbounds.full_skyline_phase import certify_complete_skyline
from foldbounds.universal_landscape import build_landscape


def _seq(mask: int, n: int) -> str:
    return "".join("H" if mask & (1 << i) else "P" for i in range(n))


def _class_winners(landscape, h_count: int, q: Fraction) -> tuple[str, ...]:
    masks = tuple(sum(1 << i for i in c)
                  for c in combinations(range(landscape.length), h_count))
    values = landscape.all_partitions(q)
    maximum = max(values[m] for m in masks)
    return tuple(_seq(m, landscape.length) for m in masks
                 if values[m] == maximum)


def _certified_cell(cert, q: Fraction) -> tuple[str, ...]:
    t = q - 1
    if any(root.lo <= t <= root.hi for root in cert.events):
        raise ValueError("witness inside critical root enclosure")
    return cert.cell_winners[sum(root.hi < t for root in cert.events)]


@pytest.mark.parametrize("n,total", [
    (2, 1), (3, 5), (4, 25), (5, 121), (6, 589),
    (7, 2821), (8, 13565), (9, 64661), (10, 308981),
])
def test_cubic_counts_and_full_independent_geometry(n, total):
    atlas = build_cubic_atlas(n)
    assert atlas.conformations == total
    assert audit_cubic_atlas(atlas, max_length=n)
    assert atlas.spectrum("P"*n) == {0: total}
    assert atlas.multiplicities[0][1] % 4 == 1
    assert all(count % 4 == 0 for mask, count in atlas.multiplicities if mask)


@pytest.fixture(scope="module")
def audited_11():
    atlas = build_cubic_atlas(11)
    assert atlas.conformations == 1_468_313
    assert len(atlas.multiplicities) == 2_180
    landscape = build_cubic_landscape(atlas, require_audit=False)
    # Independently re-enumerates ALL 1,468,313 walks and checks all
    # 2048 sequence partition-polynomial coefficient vectors against
    # a completely separate direct-score/binomial transform.
    assert audit_cubic_landscape(
        landscape, atlas, max_length=11, max_comparisons=8_000_000
    )
    return atlas, landscape


def test_all_2048_sequence_polynomials_independently_audited(audited_11):
    atlas, landscape = audited_11
    for mask in range(1 << 11):
        seq = _seq(mask, 11)
        assert landscape.contact_spectrum(seq) == atlas.spectrum(seq)
    for q in (Fraction(1), Fraction(3,2), Fraction(2), Fraction(7)):
        for seq in ("H"*11, "P"*11, "HPPHPPHPPHP", "HPPHPHPPHPP"):
            assert landscape.partition(seq, q) == atlas.partition(seq, q)


def test_all_3d_compositions_have_exact_continuum_global_design(audited_11):
    _, landscape = audited_11
    envelopes = tuple(certify_complete_skyline(landscape, k)
                      for k in range(12))
    assert sum(c.population for c in envelopes) == 2_048
    assert sum(len(c.switches) for c in envelopes) == 11
    assert len(envelopes[4].switches) == 2
    assert len(envelopes[9].switches) == 3
    for q in (Fraction(2), Fraction(10), Fraction(30)):
        for cert in envelopes:
            assert _certified_cell(cert, q) == _class_winners(
                landscape, cert.hydrophobic_count, q
            )


def test_3d_reentrant_sequence_ranking_but_2d_fixed_order(audited_11):
    cubic, l3d = audited_11
    first, second = "HPPHPPHPPHP", "HPPHPHPPHPP"
    assert first.count("H") == second.count("H") == 4
    c3 = certify_sequence_comparison(l3d, first, second,
                                     width=Fraction(1, 1_000_000))
    assert c3.positive_root_count == 2
    assert c3.cell_signs == (1, -1, 1)
    assert 3 < c3.roots[0].lo + 1 < 4
    assert 11 < c3.roots[1].lo + 1 < 12
    assert c3.roots[0].hi - c3.roots[0].lo <= Fraction(1, 1_000_000)
    assert c3.roots[1].hi - c3.roots[1].lo <= Fraction(1, 1_000_000)
    assert l3d.polynomial(first) == (1468313, 537196, 63748, 3488, 204)
    assert l3d.polynomial(second) == (1468313, 498176, 76496, 4744)
    for q, sign in ((Fraction(2), 1), (Fraction(5), -1),
                    (Fraction(10), -1), (Fraction(20), 1)):
        difference = cubic.partition(first, q) - cubic.partition(second, q)
        assert (difference > 0) - (difference < 0) == sign

    two_d = build_atlas(11)
    assert audit_atlas(two_d)
    l2d = build_landscape(two_d)
    c2 = certify_sequence_comparison(l2d, first, second)
    assert c2.positive_root_count == 0
    assert c2.cell_signs == (1,)
    assert l2d.polynomial(first) == (11025, 4510, 550, 24, 2)
    assert l2d.polynomial(second) == (11025, 3972, 592, 28)
    q = Fraction(10)
    winners_2d = _class_winners(l2d, 4, q)
    winners_3d = _class_winners(l3d, 4, q)
    assert winners_2d == ("HPPHPPHPPHP", "PHPPHPPHPPH")
    assert winners_3d == ("HPPHPHPPHPP", "PPHPPHPHPPH")
    assert set(winners_2d).isdisjoint(winners_3d)


@pytest.fixture(scope="module")
def cubic_12():
    # The n=12 geometry is independently audited *again* by separate C++
    # full enumeration in test_cubic12_independent_cpp_histogram below.
    atlas = build_cubic_atlas(12, max_nodes=3_000_000)
    assert atlas.conformations == 6_989_025
    assert len(atlas.multiplicities) == 8_738
    landscape = build_cubic_landscape(
        atlas, require_audit=False, max_length=12,
        max_edge_subsets=3_000_000
    )
    return atlas, landscape


def test_cubic12_independent_cpp_histogram(cubic_12, tmp_path):
    compiler = shutil.which("g++")
    if not compiler:
        pytest.skip("independent C++ audit requires g++")
    path = Path(__file__).resolve().parents[1] / "experiments" / "independent_cubic_audit.cpp"
    executable = tmp_path / "cubic_audit"
    subprocess.run(
        [compiler, "-O3", "-std=c++17", str(path), "-o", str(executable)],
        check=True, capture_output=True, text=True, timeout=45
    )
    run = subprocess.run(
        [str(executable), "12"], check=True,
        capture_output=True, text=True, timeout=45
    )
    independent = tuple(tuple(map(int, line.split()))
                        for line in run.stdout.splitlines())
    atlas, _ = cubic_12
    assert independent == atlas.multiplicities
    assert sum(c for _, c in independent) == 6_989_025


def test_cubic12_complete_all_4096_sequence_phase_classes(cubic_12):
    atlas, landscape = cubic_12
    result = tuple(certify_complete_skyline(landscape, k)
                   for k in range(13))
    assert sum(c.population for c in result) == 4_096
    assert sum(len(c.switches) for c in result) == 14
    assert len(result[4].switches) == 2
    assert len(result[7].switches) == 3
    assert len(result[9].switches) == 3
    for q in (Fraction(2), Fraction(10), Fraction(30)):
        for cert in result:
            assert _certified_cell(cert, q) == _class_winners(
                landscape, cert.hydrophobic_count, q
            )
    # Independent geometry-backed scoring of each 3D finalist at q=2.
    for seq in ("H"*12, "P"*12, "HPPHPPHPPHPP",
                "HPPHPHPPHPPP", "HPHHHPHPPHPH"):
        assert landscape.partition(seq, Fraction(2)) == atlas.partition(
            seq, Fraction(2)
        )


def test_fail_closed_3d_atlas_and_continuum_budget():
    atlas = build_cubic_atlas(7)
    assert audit_cubic_atlas(atlas)
    fake = CubicContactAtlas(4, ((0, 1),), 1)
    assert not audit_cubic_atlas(fake)
    with pytest.raises(ValueError):
        build_cubic_atlas(7, max_nodes=1)
    with pytest.raises(ValueError):
        audit_cubic_atlas(atlas, max_length=6)
    with pytest.raises(ValueError):
        build_cubic_landscape(atlas, max_edge_subsets=1)
    with pytest.raises(ValueError):
        build_cubic_landscape(atlas, max_length=6)
    with pytest.raises(TypeError):
        build_cubic_landscape(build_atlas(7))
    with pytest.raises(ValueError):
        CubicContactAtlas(4, ((0, 2),), 2)
    landscape = build_cubic_landscape(atlas)
    for opts in (
        {"max_sequences": 1}, {"max_classes": 1},
        {"max_events": 1}, {"width": 0.1}, {"max_refinements": 0},
    ):
        with pytest.raises(ValueError):
            certify_complete_skyline(landscape, 3, **opts)
