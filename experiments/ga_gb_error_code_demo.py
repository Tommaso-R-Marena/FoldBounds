"""Real GA95/GB95 experimental NMR reporter-code audit.

Usage:
python experiments/ga_gb_error_code_demo.py --ga data/2KDL.pdb --gb data/2KDM.pdb

Requires EXACT immutable expected SHA-256 of official PDB coordinate files.
The 20 deposited NMR models are NOT independent thermodynamic samples.
"""
import argparse
import json
from dataclasses import asdict
from fractions import Fraction
from itertools import combinations
from math import isqrt
from pathlib import Path

from foldbounds.ga_gb_error_codes import (
    EXPECTED, _d2, certify_code, classify_readouts,
    load_real_models, verify_code
)


def audit_actual_models(cert,ga,gb):
    """For every deposited model, test every pattern of <=f bad channels."""
    good=0
    bad=0
    for truth,ensemble in (("GA95",ga),("GB95",gb)):
        for model in ensemble:
            clean={e.endpoints:Fraction(isqrt(_d2(model.ca[e.a],model.ca[e.b])))
                   for e in cert.reporters}
            if classify_readouts(cert,clean)!=truth:
                raise AssertionError("deposited model classified incorrectly")
            good+=1
            for k in range(cert.corrupted_channels+1):
                for corrupted in combinations(clean,k):
                    adversarial=dict(clean)
                    for pair in corrupted:
                        adversarial[pair]=Fraction(100000 if truth=="GA95" else 0)
                    if classify_readouts(cert,adversarial)!=truth:
                        raise AssertionError("worst-case adversarial bit error not corrected")
                    bad+=1
    return {"all_40_deposited_models_correct":good,
            "all_adversarial_bit_flip_trials_correct":bad,
            "not_prospective_or_equilibrium_prediction":True}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--ga",type=Path,required=True)
    parser.add_argument("--gb",type=Path,required=True)
    parser.add_argument("--out",type=Path)
    args=parser.parse_args()
    ga,gb=load_real_models(args.ga.read_bytes(),args.gb.read_bytes())
    results={}
    for channels in (5,7):
        cert=certify_code(ga,gb,channels=channels)
        if not verify_code(cert,ga,gb):
            raise AssertionError("certificate verifier rejected real input")
        observed=audit_actual_models(cert,ga,gb)
        results[str(channels)]={
            "residue_disjoint_reporters":[asdict(e) for e in cert.reporters],
            "candidate_count":cert.all_candidate_count,
            "exact_bottleneck_gap_A":cert.bottleneck_milliang/1000,
            "worst_case_pair_distance_error_A_strictly_less_than":
                float(cert.classifier_error_milliang_strict_bound/1000),
            "per_CA_site_error_A_if_both_sites_perturbed_strictly_less_than":
                float(cert.classifier_error_milliang_strict_bound/2000),
            "arbitrarily_wrong_channels_correctable":cert.corrupted_channels,
            "stronger_margin_vertex_cover":{
                "left_sites":list(cert.stronger_graph_cover_left),
                "right_sites":list(cert.stronger_graph_cover_right),
                "size":len(cert.stronger_graph_cover_left)+len(cert.stronger_graph_cover_right),
            },
            "experiment":observed,
        }
    report={
        "PDB_ensemble_sources":EXPECTED,
        "entries":["2KDL","2KDM"],
        "models_per_fold":20,
        "distance_inputs":"exact integer milli-Angstrom coordinates from official CA PDB columns",
        "guarantee_domain":"only 40 deposited NMR models and modeled bounded pair readout errors",
        "important_limitations":[
            "No new amino acid variants, independent test proteins, or live experiments",
            "Structural NMR models are not unbiased equilibrium population samples",
            "Raw protein CA distances are not automatically experimentally measured dye distances",
            "Independent fold-validation and dye errors require new experimental controls",
            "Matching optimality is certified for the explicitly defined bipartite 1..14 vs 15..56, 12..40 A residue-pair domain",
        ],
        "codes":results,
    }
    text=json.dumps(report,indent=2,sort_keys=True)
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(text+"\n")
    print(text)


if __name__=="__main__":
    main()
