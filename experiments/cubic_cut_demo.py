"""Produce and independently verify non-enumerative 3D HP gap certificates.

Usage:
  python experiments/cubic_cut_demo.py --output-dir ./cubic_cut_proof

Requires C++17 compiler (g++). The external checker validates *every*
producer cut against exact integer contact-capacity and discrete
Loomis-Whitney projection bounds, then scores every surviving leaf.

These are 3D HP *toy-model* claims, not experimental protein predictions.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess

from foldbounds.cubic_cut_certificate import (
    certify_native_gap, contact_score, inspect_cut_header,
    verify_cut_certificate,
)

SEQUENCE = "HHHPPHHHHHHHPPHH"
# Explicit 2x2x4 Hamiltonian paths; first oriented bond (0,0,0)->(1,0,0).
# Independently scored: 11 and 9 nonconsecutive H-H contacts, respectively.
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


def compile_tools(output: Path) -> tuple[Path,Path]:
    compiler = shutil.which("g++")
    if not compiler:
        raise RuntimeError("C++17 g++ is required for separate proof producer and checker")
    source_dir = Path(__file__).resolve().parent
    executables = []
    for name in ("cubic_cut_producer", "cubic_cut_checker"):
        executable = output/name
        subprocess.run(
            [compiler,"-std=c++17","-O3","-Wall","-Wextra","-Werror",
             str(source_dir/f"{name}.cpp"),"-o",str(executable)],
            check=True,capture_output=True,text=True,timeout=60,
        )
        executables.append(executable)
    return tuple(executables)


def run_demo(output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    producer, checker = compile_tools(output)
    if contact_score(SEQUENCE, GROUND_WITNESS) != 11:
        raise RuntimeError("explicit ground-state witness has incorrect contact score")
    if contact_score(SEQUENCE, EXCITED_WITNESS) != 9:
        raise RuntimeError("explicit excited-state witness has incorrect contact score")

    filenames = {10:"n16_threshold10.cut",11:"n16_threshold11.cut",12:"n16_threshold12.cut"}
    paths = tuple(output/filenames[i] for i in (10,11,12))
    for threshold, file in zip((10,11,12), paths):
        subprocess.run(
            [str(producer.resolve()),SEQUENCE,str(threshold),str(file.resolve()),"geometry"],
            check=True,capture_output=True,text=True,timeout=60,
        )
    certificate = certify_native_gap(paths,checker,EXCITED_WITNESS)
    tails = {}
    for threshold, file in zip((10,11,12), paths):
        # Trust these header fields ONLY because certify_native_gap has already
        # called the independent checker successfully on all three files.
        checked = inspect_cut_header(file)
        tails[str(threshold)] = {
            "qualifying_conformations":checked.qualifying_conformations,
            "visited_nodes":checked.total_nodes,
            "checked_cut_prefixes":checked.cuts,
            "raw_certificate_bytes":file.stat().st_size,
            "sha256":sha256(file.read_bytes()).hexdigest(),
        }

    if (certificate.ground_contacts,certificate.ground_conformations,certificate.gap)!= (11,2032,2):
        raise RuntimeError("replayed native energy gap differs from expected oracle")
    population = {}
    for q in (2,100,1000,2000,5000,5200,10000):
        value=certificate.lower_native_probability(Fraction(q))
        population[str(q)]={"exact_numerator":value.numerator,
                            "exact_denominator":value.denominator,
                            "lower_bound":float(value)}
    targets={"50_percent": certificate.minimum_integer_weight(Fraction(1,2)),
             "90_percent": certificate.minimum_integer_weight(Fraction(9,10)),
             "99_percent":certificate.minimum_integer_weight(Fraction(99,100))}
    manifest = {
        "model":"3D cubic-lattice HP with first directed bond fixed +x",
        "scientific_limitation":"synthetic finite toy model, not a biological prediction",
        "sequence":SEQUENCE,
        "chain_length":len(SEQUENCE),
        "hydrophobic_count":SEQUENCE.count("H"),
        "ground_contacts":certificate.ground_contacts,
        "ground_conformations":certificate.ground_conformations,
        "exact_contact_gap":certificate.gap,
        "maximum_anchored_3d_walks":certificate.max_total_conformations,
        "coarse_bound_formula":"P_ground(q) >= 2032*q^2 / (2032*q^2 + 5^14 - 2032), for every real q>=1",
        "strong_bound_formula":"P_ground(q) >= 2032*q^11/(2032*q^11 + min( (5^14-2032)*q^9, sum_k C10[k]*q^k, sum_k C11[k]*q^k )), for all real q>=1",
        "verified_excited_upper_bins": [
            {str(k):v for k,v in envelope}
            for envelope in certificate.upper_non_ground_envelopes
        ],
        "tail_proofs":tails,
        "native_probability_lower_bounds":population,
        "integer_weight_targets":targets,
        "checker":"separately compiled C++ checker; no Python/Lean kernel execution proof",
        "ground_witness_coordinates":GROUND_WITNESS,
        "excited_witness_coordinates":EXCITED_WITNESS,
    }
    (output/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"verified_ground_count":certificate.ground_conformations,
                     "exact_contact_gap":certificate.gap,
                     "weight_for_90_percent_lower_bound":targets["90_percent"],
                     "cut_files":tails},indent=2),flush=True)
    return manifest


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args()
    run_demo(args.output_dir)
