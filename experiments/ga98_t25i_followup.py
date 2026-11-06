"""Unretuned GA95/GB95 NMR reporter codes on second two-state mutation cascade.

Experimental structure labels come from He et al. Structure (2012):
2LHG = GB98-T25I, three-helix fold, 10 NMR models
2LHE = GB98-T25I,L20A, beta-rich fold, 20 NMR models

All pairs, bands, and thresholds are FROZEN from reference GA95/GB95.
Not a sequence-to-fold predictor, FRET measurement, or independent
thermodynamic ensemble-population assay.
"""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import isqrt

from foldbounds.ga_gb_nmr import GB95, mutate, parse_pdb
from foldbounds.ga_gb_error_codes import _d2, classify_readouts
from experiments.ga98_gb98_transfer import (
    GB98, verify_frozen_reference, _geometric_class,
)

EXPECTED_FOLLOWUP_SHA256 = {
    "2LHG": "946bb4d44567de7725aeecd334d7c1cec086737d6445d157cd94976326c84cff",
    "2LHE": "cd2fa17ac762106c4f21c464bd7709cb612e3cd96d227ee18fdccd99608e83d1",
}

GB98_T25I=mutate(GB98,{25:("T","I")})
GB98_T25I_L20A=mutate(GB98_T25I,{20:("L","A")})


def analyze_followup(ga95_raw,gb95_raw,alpha_t25i_raw,beta_t25i_l20a_raw,
                     reference_json):
    original=verify_frozen_reference(ga95_raw,gb95_raw,reference_json)
    for entry,raw in (("2LHG",alpha_t25i_raw),("2LHE",beta_t25i_l20a_raw)):
        if sha256(raw).hexdigest()!=EXPECTED_FOLLOWUP_SHA256[entry]:
            raise ValueError(f"{entry}: unrecognized external NMR coordinate bytes")
    alpha=parse_pdb(alpha_t25i_raw,"2LHG",GB98_T25I,expected_models=10)
    beta=parse_pdb(beta_t25i_l20a_raw,"2LHE",GB98_T25I_L20A,expected_models=20)
    if sum(a!=b for a,b in zip(GB98_T25I,GB98_T25I_L20A))!=1:
        raise AssertionError("followup proteins must differ at exactly one position")
    outputs={}
    for k,cert in sorted(original.items()):
        each={}
        for entry,models,truth in (
            ("2LHG",alpha,"GA95"),("2LHE",beta,"GB95")
        ):
            correct=0
            robust=0
            trials=0
            minimum_signed_margin=None
            incorrect=[]
            channel_errors=[]
            for m_id,m in enumerate(models,1):
                classification,raw,margin,votes=_geometric_class(cert,m)
                correct+=int(classification==truth)
                if classification!=truth:
                    incorrect.append(m_id)
                wrong_pairs=[list(e.endpoints) for e in cert.reporters
                             if (e.shorter if
                                 Fraction(isqrt(_d2(m.ca[e.a],m.ca[e.b]))) <
                                 (Fraction(e.ga_range_milliang[1]+e.gb_range_milliang[0],2)
                                  if e.shorter=="GA95" else
                                  Fraction(e.gb_range_milliang[1]+e.ga_range_milliang[0],2))
                                 else ("GB95" if e.shorter=="GA95" else "GA95"))!=truth]
                channel_errors.append({"model":m_id,"incorrect_reporter_pairs":wrong_pairs})
                for e in cert.reporters:
                    d2=_d2(m.ca[e.a],m.ca[e.b])
                    lower=isqrt(d2)
                    upper=lower if lower*lower==d2 else lower+1
                    threshold=(Fraction(e.ga_range_milliang[1]+e.gb_range_milliang[0],2)
                               if e.shorter=="GA95"
                               else Fraction(e.gb_range_milliang[1]+e.ga_range_milliang[0],2))
                    signed=(threshold-upper if truth==e.shorter else lower-threshold)
                    minimum_signed_margin=(signed if minimum_signed_margin is None
                                           else min(minimum_signed_margin,signed))
                by_pair={e.endpoints:e for e in cert.reporters}
                for faults in range(cert.corrupted_channels+1):
                    for corrupted in combinations(raw,faults):
                        perturbed=dict(raw)
                        for pair in corrupted:
                            perturbed[pair]=Fraction(
                                100000 if by_pair[pair].shorter==truth else 0
                            )
                        robust+=int(classify_readouts(cert,perturbed)==truth)
                        trials+=1
            each[entry]={
                "experimental_fold_class":truth,
                "sequence":models[0].seq,
                "deposited_models":len(models),
                "correct_clean":correct,
                "incorrect_model_ids":incorrect,
                "adversarial_trials_correct":robust,
                "adversarial_trials_total":trials,
                "minimum_true_fold_signed_margin_A":float(minimum_signed_margin/1000),
                "wrong_channel_diagnostics":channel_errors,
            }
        outputs[str(k)]={
            "frozen_reference_gap_A":cert.bottleneck_milliang/1000,
            "permitted_corrupted_channel_count":cert.corrupted_channels,
            "correct_clean_total":sum(v["correct_clean"] for v in each.values()),
            "total_deposited_models":sum(v["deposited_models"] for v in each.values()),
            "adversarial_trials_correct":sum(v["adversarial_trials_correct"] for v in each.values()),
            "adversarial_trials_total":sum(v["adversarial_trials_total"] for v in each.values()),
            "per_entry":each,
        }
    return {
        "experimental_depositions":["2LHG","2LHE"],
        "experimental_fold_annotations":{
            "2LHG":"3-alpha (GA-like); assigned by published 2012 NMR study",
            "2LHE":"4-beta-plus-alpha (GB-like); assigned by published 2012 NMR study",
        },
        "PDB_SHA256":{
            "2LHG":sha256(alpha_t25i_raw).hexdigest(),
            "2LHE":sha256(beta_t25i_l20a_raw).hexdigest(),
        },
        "mutations_in_second_cascade":{
            "GB98_to_2LHG":[{"site":25,"from":"T","to":"I"}],
            "2LHG_to_2LHE":[{"site":20,"from":"L","to":"A"}],
        },
        "external_models":30,
        "training_structures_unchanged_and_never_refitted":True,
        "retrospective_experimental_structural_classification_not_sequence_prediction":True,
        "codes":outputs,
    }
