"""Authentic 2KDL/2KDM-only training -> mutation-site-free GA/GB NMR validation.

Run with six official PDB files. This is a retrospective structural
recognition benchmark and exhaustively audits channel erasure plus
adversarial errors 2*t+e<K, not a physical FRET assay.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from itertools import combinations
from math import isqrt
from pathlib import Path
from hashlib import sha256

from foldbounds.ga_gb_error_codes import (
    EXPECTED, load_real_models, _d2, _ceil_sqrt
)
from foldbounds.ga_gb_mutation_safe import (
    FORBIDDEN_MUTATION_SITES, TARGET_LENGTHS,
    make_safe_code, check_safe_code, decode_with_erasures,
    certified_budget, coordinate_distance_readings,
    finite_margin_to_truth
)
from experiments.ga98_gb98_transfer import (
    GA98,GB98,EXPECTED_TEST_SHA256,trusted_test_structures
)
from experiments.ga98_t25i_followup import (
    GB98_T25I, GB98_T25I_L20A, EXPECTED_FOLLOWUP_SHA256
)
from foldbounds.ga_gb_nmr import GA95,GB95,parse_pdb


def _read_and_authenticate(path:Path, expected_sha:str, name:str):
    raw=path.read_bytes()
    actual=sha256(raw).hexdigest()
    if actual!=expected_sha:
        raise ValueError(f"{name}: untrusted PDB SHA256 {actual}")
    return raw


def _read_models(files:dict[str,Path]):
    names=set(EXPECTED)|set(EXPECTED_TEST_SHA256)|set(EXPECTED_FOLLOWUP_SHA256)
    if set(files)!=names:
        raise ValueError("require all six explicitly declared RCSB entries")
    raws={}
    for k,p in files.items():
        known={**EXPECTED,**EXPECTED_TEST_SHA256,**EXPECTED_FOLLOWUP_SHA256}
        raws[k]=_read_and_authenticate(p,known[k],k)
    train_a,train_b=load_real_models(raws["2KDL"],raws["2KDM"])
    ext_a,ext_b=trusted_test_structures(raws["2LHC"],raws["2LHD"])
    held_alpha=parse_pdb(raws["2LHG"],"2LHG",GB98_T25I,expected_models=10)
    held_beta=parse_pdb(raws["2LHE"],"2LHE",GB98_T25I_L20A,expected_models=20)
    seqs=[GA95,GB95,GA98,GB98,GB98_T25I,GB98_T25I_L20A]
    for site in range(1,57):
        if site in FORBIDDEN_MUTATION_SITES:
            continue
        if len({s[site-1] for s in seqs})!=1:
            raise AssertionError(f"training/test amino acid differs at a supposedly safe site {site}")
    external={
        "2LHC":("GA95",ext_a),
        "2LHD":("GB95",ext_b),
        "2LHG":("GA95",held_alpha),
        "2LHE":("GB95",held_beta),
    }
    return train_a,train_b,external


def _wrong_class_value(edge,truth):
    """For this edge, force wrong class regardless of which fold is short."""
    return Fraction(100000 if edge.shorter==truth else 0)


def _trial_summary(code,ensemble,truth,*,max_erased=2):
    k=code.reference.target_channels
    features={edge.endpoints:edge for edge in code.reference.reporters}
    correct=0; all_trials=0; passed=0
    ood_channels=0; ood_models=0
    worst_margin=None; ties=0
    failures=[]
    per_model=[]
    for model_id,m in enumerate(ensemble,1):
        raw=coordinate_distance_readings(code,m)
        label=decode_with_erasures(code,raw)
        correct+=int(label==truth)
        ties+=int(label is None)
        margin=finite_margin_to_truth(code,m,truth)
        worst_margin=margin if worst_margin is None else min(worst_margin,margin)
        outlier_pairs=[]
        for edge in code.reference.reporters:
            d2=_d2(m.ca[edge.a],m.ca[edge.b])
            lower,upper=isqrt(d2),_ceil_sqrt(d2)
            band=(edge.ga_range_milliang if truth=="GA95"
                  else edge.gb_range_milliang)
            if lower<band[0] or upper>band[1]:
                outlier_pairs.append([edge.a,edge.b])
        ood_channels+=len(outlier_pairs)
        ood_models+=int(bool(outlier_pairs))
        trials=0; successes=0
        edges=list(code.reference.reporters)
        for e in range(min(max_erased,k)+1):
            for erased in combinations(edges,e):
                erased_pairs={x.endpoints for x in erased}
                remaining=[edge for edge in edges if edge.endpoints not in erased_pairs]
                for t in range(len(remaining)+1):
                    if not certified_budget(k,e,t):
                        continue
                    for wrong in combinations(remaining,t):
                        observations={
                            pair:value for pair,value in raw.items()
                            if pair not in erased_pairs
                        }
                        for edge in wrong:
                            observations[edge.endpoints]=_wrong_class_value(edge,truth)
                        value=decode_with_erasures(code,observations)
                        is_right=value==truth
                        successes+=int(is_right)
                        trials+=1
                        if not is_right and len(failures)<50:
                            failures.append({
                                "model_id":model_id,"missing_pairs":sorted([list(p) for p in erased_pairs]),
                                "wrong_pairs":[list(edge.endpoints) for edge in wrong],
                                "outcome":value,
                            })
        passed+=successes;all_trials+=trials
        per_model.append({
            "id":model_id,
            "clean_decision":label,
            "correct_clean":label==truth,
            "certified_true_class_min_margin_mA":[margin.numerator,margin.denominator],
            "out_of_reference_distance_envelope_pairs":outlier_pairs,
            "correct_error_erasure_trials":successes,
            "total_error_erasure_trials":trials,
        })
    return {
        "models":len(ensemble),
        "truth_fold":truth,
        "clean_correct":correct,
        "abstentions":ties,
        "error_erasure_trials_passed":passed,
        "error_erasure_trials_total":all_trials,
        "minimum_true_fold_margin_A":float(worst_margin/1000),
        "models_with_any_channel_outside_training_distance_band":ood_models,
        "individual_channels_outside_training_distance_band":ood_channels,
        "first_50_failures":failures,
        "per_model":per_model,
    }


def audit(files:dict[str,Path]):
    train_a,train_b,external=_read_models(files)
    result={
        "training_entries":["2KDL","2KDM"],
        "heldout_entries":["2LHC","2LHD","2LHG","2LHE"],
        "official_SHA256":{**EXPECTED,**EXPECTED_TEST_SHA256,**EXPECTED_FOLLOWUP_SHA256},
        "immutable_mutation_exclusion_sites":sorted(FORBIDDEN_MUTATION_SITES),
        "excluded_sites_identical_across_all_six_seqs":False,
        "remaining_reporter_site_residues_identical_across_all_six_sequences":True,
        "site_exposure_or_fluorophore_chemistry_verified":False,
        "max_erasures_tested":2,
        "sampling_warning":"deposited NMR conformers are not unbiased independent solution-state samples",
        "codes":{},
    }
    for k in TARGET_LENGTHS:
        try:
            code=make_safe_code(train_a,train_b,k)
        except ValueError as ex:
            result["codes"][str(k)]={"abstained":str(ex)}
            continue
        if not check_safe_code(code,train_a,train_b):
            raise AssertionError("mutation-safe training checker rejected candidate")
        chosen=code.reference
        evaluations={}
        for entry,(truth,models) in external.items():
            evaluations[entry]=_trial_summary(code,models,truth)
        result["codes"][str(k)]={
            "mutation_safe_pairs":[{
                "a":e.a,"b":e.b,"gap_mA":e.gap_milliang,
                "shorter":e.shorter,"GA_band_mA":list(e.ga_range_milliang),
                "GB_band_mA":list(e.gb_range_milliang),
            } for e in chosen.reporters],
            "eligible_training_candidates":code.eligible_edges,
            "optimal_training_bottleneck_gap_A":chosen.bottleneck_milliang/1000,
            "source_only_stronger_margin_vertex_cover":{
                "left":list(chosen.stronger_graph_cover_left),
                "right":list(chosen.stronger_graph_cover_right),
            },
            "correctable_binary_channel_errors_if_no_erasures":(k-1)//2,
            "max_test_erasure_sites":2,
            "heldout":{
                "clean_correct":sum(v["clean_correct"] for v in evaluations.values()),
                "total_models":sum(v["models"] for v in evaluations.values()),
                "error_erasure_trials_passed":sum(v["error_erasure_trials_passed"] for v in evaluations.values()),
                "error_erasure_trials_total":sum(v["error_erasure_trials_total"] for v in evaluations.values()),
                "models_outside_original_distance_bands":sum(v["models_with_any_channel_outside_training_distance_band"]
                                                             for v in evaluations.values()),
                "minimum_true_fold_margin_A":min(v["minimum_true_fold_margin_A"] for v in evaluations.values()),
                "by_entry":evaluations,
            },
        }
    return result


def main():
    parser=argparse.ArgumentParser()
    for entry in ("2KDL","2KDM","2LHC","2LHD","2LHG","2LHE"):
        parser.add_argument("--"+entry.lower(),type=Path,required=True)
    parser.add_argument("--out",type=Path)
    a=parser.parse_args()
    files={entry:getattr(a,entry.lower()) for entry in
           ("2KDL","2KDM","2LHC","2LHD","2LHG","2LHE")}
    output=json.dumps(audit(files),sort_keys=True,indent=2)+"\n"
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        a.out.write_text(output)
    print(output)


if __name__=="__main__":
    main()
