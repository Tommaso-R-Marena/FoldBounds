"""Mutation-heldout GA95/GB95 -> GA98/GB98 structural transfer audit.

Reference 2KDL/2KDM reporter pairs and thresholds are FROZEN in a
versioned SHA-pinned JSON before 2LHC/2LHD coordinates are inspected.

This evaluates known experimental structural folds, NOT sequence-based
fold prediction or solution-state population inference. It retains all
misclassifications and geometric tolerance violations, including
adversarially forced wrong readout decisions.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import json

from foldbounds.ga_gb_nmr import GA95, GB95, mutate, parse_pdb
from foldbounds.ga_gb_error_codes import (
    EXPECTED, CodeCertificate, Reporter, classify_readouts, _d2,
    _ceil_sqrt, certify_code, verify_code,
)

GA98 = mutate(GA95, {30: ("I", "F")})
GB98 = mutate(GB95, {20: ("A", "L")})
assert sum(a != b for a,b in zip(GA98,GB98)) == 1
assert GA98[44] == "L" and GB98[44] == "Y"


def verify_frozen_reference(ga95_data: bytes, gb95_data: bytes, reference_json: Path):
    """Rebuild original codes from authentic 2KDL/2KDM only; no target access."""
    report=json.loads(reference_json.read_text(encoding="utf-8"))
    if report["PDB_ensemble_sources"] != EXPECTED:
        raise ValueError("frozen source hashes have changed")
    if set(report["codes"]) != {"1","3","5","7","9"}:
        raise ValueError("frozen reporter family or version mismatch")
    for key,raw in (("2KDL",ga95_data),("2KDM",gb95_data)):
        if sha256(raw).hexdigest()!=EXPECTED[key]:
            raise ValueError(f"{key}: unauthorized reference coordinate bytes")
    ga=parse_pdb(ga95_data,"2KDL",GA95)
    gb=parse_pdb(gb95_data,"2KDM",GB95)
    result={}
    for k in (1,3,5,7,9):
        cert=certify_code(ga,gb,channels=k)
        if not verify_code(cert,ga,gb):
            raise ValueError("reference code failed independent verification")
        original=report["codes"][str(k)]
        reproduced_pairs=[(x.a,x.b,x.gap_milliang,x.shorter,
                           x.ga_range_milliang,x.gb_range_milliang)
                          for x in cert.reporters]
        frozen_pairs=[(x["a"],x["b"],x["gap_milliang"],x["shorter"],
                       tuple(x["ga_range_milliang"]),tuple(x["gb_range_milliang"]))
                      for x in original["residue_disjoint_reporters"]]
        if (reproduced_pairs!=frozen_pairs
                or cert.bottleneck_milliang/1000 != original["exact_bottleneck_gap_A"]
                or cert.all_candidate_count!=original["candidate_count"]):
            raise ValueError(f"frozen {k}-channel reference witness changed")
        result[k]=cert
    return result


def trusted_test_structures(ga98_data:bytes,gb98_data:bytes):
    """Check PDB identity, all 56 residue identities, full NMR ensembles."""
    if not isinstance(ga98_data,bytes) or not isinstance(gb98_data,bytes):
        raise ValueError("expected raw PDB coordinate bytes")
    ga=parse_pdb(ga98_data,"2LHC",GA98,expected_models=20)
    gb=parse_pdb(gb98_data,"2LHD",GB98,expected_models=20)
    if any(m.seq!=GA98 for m in ga) or any(m.seq!=GB98 for m in gb):
        raise ValueError("nonidentical GA98/GB98 variants")
    return ga,gb


def _geometric_class(cert:CodeCertificate, structure):
    """Actual structure-derived readouts with rigorous ±1 milli-Å brackets.

    If integer sqrt brackets straddle a classification threshold, reject
    as ambiguous rather than force a potentially incorrect direction.
    """
    raw={}
    margins={}
    for e in cert.reporters:
        d2=_d2(structure.ca[e.a],structure.ca[e.b])
        low=isqrt(d2); high=_ceil_sqrt(d2)
        if e.shorter=="GA95":
            middle=Fraction(e.ga_range_milliang[1]+e.gb_range_milliang[0],2)
            low_fold="GA95"
        else:
            middle=Fraction(e.gb_range_milliang[1]+e.ga_range_milliang[0],2)
            low_fold="GB95"
        if high<middle:
            vote=low_fold
        elif low>middle:
            vote="GB95" if low_fold=="GA95" else "GA95"
        else:
            raise ValueError("classification requires unattained sub-milliangstrom precision")
        raw[e.endpoints]=Fraction(low)
        margins[e.endpoints]=(
            middle-Fraction(high) if vote==low_fold
            else Fraction(low)-middle
        )
    direct=classify_readouts(cert,raw)
    votes={"GA95":0,"GB95":0}
    for e in cert.reporters:
        d2=_d2(structure.ca[e.a],structure.ca[e.b])
        low=isqrt(d2); high=_ceil_sqrt(d2)
        mid=Fraction(
            (e.ga_range_milliang[1]+e.gb_range_milliang[0])
            if e.shorter=="GA95"
            else (e.gb_range_milliang[1]+e.ga_range_milliang[0]),2
        )
        short=e.shorter
        v=(short if high<mid else ("GB95" if short=="GA95" else "GA95"))
        votes[v]+=1
    if direct!=("GA95" if votes["GA95"]>votes["GB95"] else "GB95"):
        raise AssertionError("independent vote reader disagrees with classifier")
    return direct,raw,margins,votes


def heldout_transfer(codes:dict[int,CodeCertificate],ga98,gb98):
    output={}
    for k,cert in sorted(codes.items()):
        by_fold={}
        clean_all=0
        robust_all=0
        trials_all=0
        min_margin=None
        min_good_votes=k
        for truth,entry,ensemble in (("GA95","2LHC",ga98),("GB95","2LHD",gb98)):
            good=0;robust=0;trials=0;per_sample=[]
            for model_id,model in enumerate(ensemble,1):
                predicted,readouts,margins,votes=_geometric_class(cert,model)
                clean=predicted==truth
                good+=int(clean)
                good_votes=votes[truth]
                min_good_votes=min(min_good_votes,good_votes)
                for e in cert.reporters:
                    distance_squared=_d2(model.ca[e.a],model.ca[e.b])
                    lower=isqrt(distance_squared)
                    upper=_ceil_sqrt(distance_squared)
                    threshold=Fraction(
                        e.ga_range_milliang[1]+e.gb_range_milliang[0]
                        if e.shorter=="GA95"
                        else e.gb_range_milliang[1]+e.ga_range_milliang[0],
                        2,
                    )
                    # Worst signed margin is based on TRUE fold,
                    # not on whichever vote the reporter happened
                    # to produce. Negative means an actual error.
                    signed=(threshold-Fraction(upper) if e.shorter==truth
                            else Fraction(lower)-threshold)
                    min_margin=signed if min_margin is None else min(min_margin,signed)
                model_robust=0;model_trials=0
                reporter_by_pair={e.endpoints:e for e in cert.reporters}
                for nwrong in range(cert.corrupted_channels+1):
                    for flipped in combinations(readouts,nwrong):
                        corrupted=dict(readouts)
                        for pair in flipped:
                            # Set the numerical distance to guarantee
                            # the opposite binary vote for THIS reporter.
                            corrupted[pair]=Fraction(
                                100000 if reporter_by_pair[pair].shorter==truth else 0
                            )
                        model_robust+=int(classify_readouts(cert,corrupted)==truth)
                        model_trials+=1
                robust+=model_robust;trials+=model_trials
                per_sample.append({
                    "model_id":model_id,"predicted_fold":predicted,
                    "correct":clean,"correct_reporter_votes":good_votes,
                    "total_reporters":k,
                    "worst_case_corrupted_patterns_passed":model_robust,
                    "worst_case_corrupted_patterns_total":model_trials,
                })
            by_fold[entry]={
                "experimental_structure_label":truth,
                "sequence":ensemble[0].seq,
                "models":len(ensemble),
                "clean_correct":good,
                "strict_corrupted_trials_correct":robust,
                "strict_corrupted_trials_total":trials,
                "per_model":per_sample,
            }
            clean_all+=good;robust_all+=robust;trials_all+=trials
        output[str(k)]={
            "source_code_minimum_reference_gap_A":cert.bottleneck_milliang/1000,
            "original_corruption_budget":cert.corrupted_channels,
            "transfer_clean_correct":clean_all,
            "transfer_clean_total":len(ga98)+len(gb98),
            "transfer_corrupted_trials_correct":robust_all,
            "transfer_corrupted_trials_total":trials_all,
            "minimum_correct_votes_per_model":min_good_votes,
            "minimum_signed_reporter_decision_margin_A":float(min_margin/1000),
            "by_entry":by_fold,
        }
    return output


def run(ga95_data:bytes,gb95_data:bytes,ga98_data:bytes,gb98_data:bytes,
        reference_json:Path):
    frozen=verify_frozen_reference(ga95_data,gb95_data,reference_json)
    ga98,gb98=trusted_test_structures(ga98_data,gb98_data)
    result={
        "experiment":"Fully frozen GA95/GB95 2009 NMR distance codes evaluated on independently deposited GA98/GB98 2012 NMR coordinates",
        "new_PDB_SHA256":{
            "2LHC":sha256(ga98_data).hexdigest(),
            "2LHD":sha256(gb98_data).hexdigest(),
        },
        "new_test_sequences":{"GA98":GA98,"GB98":GB98},
        "mutations_between_test_variants":[{
            "site":i+1,"GA98":a,"GB98":b
        } for i,(a,b) in enumerate(zip(GA98,GB98)) if a!=b],
        "models_per_test_entry":20,
        "test_structures_were_excluded_from_all_reporter_selection":True,
        "test_data_temporally_deposited_after_training_structures":True,
        "not_a_blind_prediction_of_sequences_or_population":True,
        "codes":heldout_transfer(frozen,ga98,gb98),
    }
    return result
