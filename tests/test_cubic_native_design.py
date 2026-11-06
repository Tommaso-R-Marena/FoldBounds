"""Independently audited cubic 3D symmetry-unique ground-state discovery."""
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
import shutil
import subprocess

import pytest

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.cubic_contact_atlas import build_cubic_atlas
from foldbounds.cubic_native_design import (
    CubicNativeWitness, UniqueCubicNativeOrbit,
    find_unique_cubic_native_orbits, recover_cubic_native_witnesses,
    native_ground_probability, isolate_native_occupancy_weight,
    lattice_witness_pdb,
)

EXPECTED = {
    "HHHPHPHHPHPH": 3_299_410,
    "HHPHPHHPHPHH": 26_289_413,
    "HPHPHHPHPHHH": 8_131_232,
}


@pytest.fixture(scope="module")
def independently_checked_models(tmp_path_factory):
    planar = build_atlas(12)
    assert audit_atlas(planar)
    cubic = build_cubic_atlas(12)
    assert cubic.conformations == 6_989_025
    assert len(cubic.multiplicities) == 8_738
    compiler = shutil.which("g++")
    if compiler is None:
        pytest.skip("full independent cubic n=12 audit requires g++")
    directory = tmp_path_factory.mktemp("independent_cubic_native")
    source = Path(__file__).resolve().parents[1] / "experiments" / "independent_cubic_audit.cpp"
    executable = directory / "cubic_audit"
    subprocess.run(
        [compiler, "-O3", "-std=c++17", str(source), "-o", str(executable)],
        capture_output=True, text=True, check=True, timeout=45
    )
    result = subprocess.run(
        [str(executable), "12"], capture_output=True, text=True,
        check=True, timeout=45
    )
    independent = tuple(tuple(map(int, line.split()))
                        for line in result.stdout.splitlines())
    assert independent == cubic.multiplicities
    return cubic, planar


def test_only_three_nontrivial_symmetry_unique_native_sequences(independently_checked_models):
    cubic, planar = independently_checked_models
    certified = find_unique_cubic_native_orbits(
        cubic, planar, require_audit=False
    )
    assert len(certified) == 3
    assert {c.sequence: c.ground_contact_mask for c in certified} == EXPECTED
    assert {c.hydrophobic_count for c in certified} == {8}
    assert {c.max_contacts for c in certified} == {7}
    assert {c.ground_conformations for c in certified} == {8}
    assert {c.symmetry_orbit_size for c in certified} == {8}
    # These maps have zero planar realizations. Eight anchored cubic
    # geometries are therefore precisely one nonplanar D4 orbit, not
    # two fourfold planar orbits.
    planar_masks = {m for m, _ in planar.multiplicities}
    assert all(c.ground_contact_mask not in planar_masks for c in certified)
    assert len({min(c.sequence, c.sequence[::-1]) for c in certified}) == 2
    assert "HHPHPHHPHPHH" == "HHPHPHHPHPHH"[::-1]


def test_explicit_witness_structures_are_unique_compact_prisms(independently_checked_models):
    cubic, planar = independently_checked_models
    certified = find_unique_cubic_native_orbits(
        cubic, planar, require_audit=False
    )
    witnesses = recover_cubic_native_witnesses(cubic, certified)
    assert len(witnesses) == 3
    for witness in witnesses:
        path = witness.coordinates
        assert len(path) == len(set(path)) == 12
        assert path[:2] == ((0,0,0), (1,0,0))
        sizes = tuple(max(p[j] for p in path) - min(p[j] for p in path) + 1
                      for j in range(3))
        assert sorted(sizes) == [2, 2, 3]
        assert len(set(path)) == sizes[0] * sizes[1] * sizes[2]
        # Each unique native fold fills all 12 sites of a 2x2x3 prism.
        pdb = lattice_witness_pdb(witness)
        assert "NOT AN EXPERIMENTAL OR PREDICTED REAL PROTEIN" in pdb
        assert pdb.count("\nATOM  ") == 12
        assert pdb.count("\nCONECT") == 11
        assert "ALA" in pdb and "GLY" in pdb
        assert pdb.endswith("END\n")


def test_native_probability_rational_thresholds(independently_checked_models):
    cubic, planar = independently_checked_models
    certified = find_unique_cubic_native_orbits(
        cubic, planar, require_audit=False
    )
    for c in certified:
        s = c.sequence
        assert native_ground_probability(cubic, s, Fraction(1)) == Fraction(8, 6_989_025)
        assert native_ground_probability(cubic, s, Fraction(2)) > Fraction(0)
        assert native_ground_probability(cubic, s, Fraction(2)) < Fraction(1, 1_000)
        assert native_ground_probability(cubic, s, Fraction(100)) < Fraction(1, 2)
        assert native_ground_probability(cubic, s, Fraction(200)) > Fraction(1, 2)
        for target in (Fraction(1,2), Fraction(9,10)):
            low, high = isolate_native_occupancy_weight(
                cubic, s, target, width=Fraction(1, 10_000)
            )
            assert low <= high
            assert high - low <= Fraction(1, 10_000)
            assert native_ground_probability(cubic, s, low) <= target
            assert native_ground_probability(cubic, s, high) >= target
        half_low, half_high = isolate_native_occupancy_weight(
            cubic, s, Fraction(1,2)
        )
        if s == "HHPHPHHPHPHH":
            assert Fraction(128) < half_low < half_high < Fraction(130)
        else:
            assert Fraction(134) < half_low < half_high < Fraction(135)


def test_exact_full_conformation_mass_distribution(independently_checked_models):
    cubic, planar = independently_checked_models
    native = find_unique_cubic_native_orbits(cubic, planar, require_audit=False)
    expected = {
        "HHHPHPHHPHPH": {
            0:3956309, 1:2001292, 2:721652, 3:239472,
            4:59228, 5:10068, 6:996, 7:8
        },
        "HHPHPHHPHPHH": {
            0:3262721, 1:2379240, 2:971780, 3:297184,
            4:66920, 5:10224, 6:948, 7:8
        },
    }
    for c in native:
        spectrum = cubic.spectrum(c.sequence)
        reference = expected["HHHPHPHHPHPH"] if c.sequence != "HHPHPHHPHPHH" else expected[c.sequence]
        assert spectrum == reference
        assert sum(spectrum.values()) == cubic.conformations
        assert max(spectrum) == 7
        assert spectrum[7] == 8


def test_fail_closed_audits_resources_and_witness_integrity(independently_checked_models):
    cubic, planar = independently_checked_models
    with pytest.raises(ValueError):
        find_unique_cubic_native_orbits(cubic, planar, max_sequences=2)
    with pytest.raises(ValueError):
        find_unique_cubic_native_orbits(cubic, build_atlas(11), require_audit=False)
    with pytest.raises(ValueError):
        find_unique_cubic_native_orbits(cubic, planar, require_audit=1)
    certified = find_unique_cubic_native_orbits(cubic, planar, require_audit=False)
    with pytest.raises(ValueError):
        recover_cubic_native_witnesses(cubic, certified, max_nodes=1)
    with pytest.raises(ValueError):
        recover_cubic_native_witnesses(cubic, (certified[0], certified[0]))
    with pytest.raises(ValueError):
        isolate_native_occupancy_weight(cubic, certified[0].sequence, 0.5)
    with pytest.raises(ValueError):
        isolate_native_occupancy_weight(cubic, certified[0].sequence, Fraction(1,2), max_steps=1)
    with pytest.raises(ValueError):
        native_ground_probability(cubic, certified[0].sequence, 2.0)
    witness = recover_cubic_native_witnesses(cubic, (certified[0],))[0]
    with pytest.raises(ValueError):
        lattice_witness_pdb(replace(witness, coordinates=tuple([(0,0,0)] * 12)))
    with pytest.raises(ValueError):
        lattice_witness_pdb(witness, spacing_angstrom=3.8)
