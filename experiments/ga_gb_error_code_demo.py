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
    EXPECTED, _d2, _ceil_sqrt, certify_code, classify_readouts,
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



def heldout_check(cert,ga_test,gb_test):
    """Out-of-fit NMR conformers: retain negative outcomes, never hide failures."""
    clean_hits=0
    robust_hits=0
    total_robust=0
    min_margin=None
    for truth,ensemble in (("GA95",ga_test),("GB95",gb_test)):
        for model in ensemble:
            clean={e.endpoints:Fraction(isqrt(_d2(model.ca[e.a],model.ca[e.b])))
                   for e in cert.reporters}
            clean_hits+=int(classify_readouts(cert,clean)==truth)
            for r in range(cert.corrupted_channels+1):
                for corrupted in combinations(clean,r):
                    adversarial=dict(clean)
                    for pair in corrupted:
                        adversarial[pair]=Fraction(100000 if truth=="GA95" else 0)
                    robust_hits+=int(classify_readouts(cert,adversarial)==truth)
                    total_robust+=1
            for e in cert.reporters:
                squared=_d2(model.ca[e.a],model.ca[e.b])
                if e.shorter=="GA95":
                    threshold=Fraction(e.ga_range_milliang[1]+e.gb_range_milliang[0],2)
                else:
                    threshold=Fraction(e.gb_range_milliang[1]+e.ga_range_milliang[0],2)
                if truth==e.shorter:
                    margin=threshold-_ceil_sqrt(squared)
                else:
                    margin=isqrt(squared)-threshold
                min_margin=margin if min_margin is None else min(min_margin,margin)
    return {
        "heldout_models":len(ga_test)+len(gb_test),
        "heldout_clean_correct":clean_hits,
        "adversarial_pattern_trials":total_robust,
        "adversarial_pattern_correct":robust_hits,
        "smallest_certified_signed_holdout_decision_margin_A":
            float(min_margin/1000),
        "evaluation_type":"retrospective within-deposition half-model holdout; not external validation",
    }


def alternating_split_assessment(ga,gb,channels):
    """Fixed, predeclared alternating-model holdouts. No label reuse in fit."""
    results={}
    for offset in (0,1):
        train_a,train_b=ga[offset::2],gb[offset::2]
        test_a,test_b=ga[1-offset::2],gb[1-offset::2]
        try:
            fit=certify_code(train_a,train_b,channels=channels)
            verify_code(fit,train_a,train_b)
            results[f"fit_model_indices_{'1_3_5' if offset==0 else '2_4_6'}_etc"]={
                "training_models_per_fold":len(train_a),
                "fitted_bottleneck_gap_A":fit.bottleneck_milliang/1000,
                "reporters":[list(e.endpoints) for e in fit.reporters],
                **heldout_check(fit,test_a,test_b),
            }
        except ValueError as ex:
            results[f"fit_offset_{offset}"]={
                "no_certificate":str(ex),
                "heldout_test_skipped":True,
            }
    return results

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
            "fixed_half_model_holdout":alternating_split_assessment(ga,gb,channels),
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
