"""Produce independently audited, synthetic CA-only PyMOL PDB witnesses.

Usage:
  python experiments/cubic_native_witness_demo.py --output-dir ./cubic_hp_witnesses

The program first checks *all 6,989,025 walks* against a separately
implemented full C++ walk/leaf-score auditor (if g++ exists) or an
unsymmetrized full Python walk audit. It audits the planar atlas too.
Only then does it report three 3D symmetry-unique ground-state designs.

PDBs are illustrative lattice curves (ALA=H, GLY=P), NOT real proteins.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
import shutil

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.cubic_contact_atlas import audit_cubic_atlas, build_cubic_atlas
from foldbounds.cubic_native_design import (
    find_unique_cubic_native_orbits, recover_cubic_native_witnesses,
    native_ground_probability, isolate_native_occupancy_weight,
    lattice_witness_pdb
)
from cubic_3d_phase_demo import independent_cpp_audit


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    cubic = build_cubic_atlas(12)
    planar = build_atlas(12)
    if not audit_atlas(planar):
        raise RuntimeError("independent full planar geometry audit failed")
    if shutil.which("g++"):
        if not independent_cpp_audit(cubic):
            raise RuntimeError("independent full cubic C++ geometry audit failed")
        print("Cubic cross-language histogram audit: PASS")
    else:
        if not audit_cubic_atlas(cubic, max_length=12, max_leaves=8_000_000):
            raise RuntimeError("independent full cubic Python geometry audit failed")
        print("Cubic full-walk Python audit: PASS")

    native = find_unique_cubic_native_orbits(
        cubic, planar, require_audit=False
    )
    expected = {
        "HHHPHPHHPHPH", "HHPHPHHPHPHH", "HPHPHHPHPHHH"
    }
    if {c.sequence for c in native} != expected:
        raise RuntimeError("unexpected native ground-orbit classification")
    witnesses = recover_cubic_native_witnesses(cubic, native)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = [
        "sequence,h_count,max_contacts,native_conformations,"
        "orbital_degeneracy,native_q_50_lower,native_q_50_upper,"
        "native_q_90_lower,native_q_90_upper"
    ]
    for witness in witnesses:
        c = witness.certificate
        half = isolate_native_occupancy_weight(
            cubic, c.sequence, Fraction(1, 2)
        )
        ninety = isolate_native_occupancy_weight(
            cubic, c.sequence, Fraction(9, 10)
        )
        file = args.output_dir / f"{c.sequence}_toy_lattice.pdb"
        file.write_text(lattice_witness_pdb(witness), encoding="utf-8")
        manifest.append(
            f"{c.sequence},{c.hydrophobic_count},{c.max_contacts},"
            f"{c.ground_conformations},{c.symmetry_orbit_size},"
            f"{half[0]},{half[1]},{ninety[0]},{ninety[1]}"
        )
        print(
            f"{c.sequence}: unique orbit size {c.symmetry_orbit_size}; "
            f"contacts {c.max_contacts}; P_native(q=2) "
            f"= {native_ground_probability(cubic, c.sequence, Fraction(2))}"
        )
        print(f"  P_native=50%: q in [{half[0]}, {half[1]}]")
        print(f"  P_native=90%: q in [{ninety[0]}, {ninety[1]}]")
        print("  PyMOL synthetic PDB:", file)
    (args.output_dir / "manifest.csv").write_text(
        "\n".join(manifest) + "\n", encoding="utf-8"
    )
    print(
        "\nAll PDBs are synthetic 3.8 A-spaced lattice CA traces "
        "(ALA=H, GLY=P), NOT actual proteins."
    )


if __name__ == "__main__":
    main()
