"""Exact finite-conformer geometric noise radii for mutation-safe fold codes.

The thresholds are POST-HOC on 70 independently deposited NMR
coordinates. They are mathematical finite-set robustness certificates,
not training-only uncertainty estimates, label precision guarantees,
or prospective solution-phase population evidence.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path

from foldbounds.ga_gb_mutation_safe import make_safe_code,check_safe_code
from foldbounds.mutation_safe_radius import (
    risk_channels_for_budget, sorted_signed_margins,
    certify_perturbation_radius,check_perturbation_certificate
)
from experiments.mutation_safe_real_audit import _read_models


def as_ratio(v:Fraction):
    return {"numerator":v.numerator,"denominator":v.denominator,
            "Angstrom":float(v/1000)}


def main():
    p=argparse.ArgumentParser()
    for pdb in ("2KDL","2KDM","2LHC","2LHD","2LHG","2LHE"):
        p.add_argument("--"+pdb.lower(),type=Path,required=True)
    p.add_argument("--out",type=Path)
    args=p.parse_args()
    files={name:getattr(args,name.lower()) for name in
           ("2KDL","2KDM","2LHC","2LHD","2LHG","2LHE")}
    train_a,train_b,holdout=_read_models(files)
    conformers=tuple((truth,model) for entry,(truth,models) in holdout.items()
                     for model in models)
    tags=tuple((entry,i) for entry,(_,models) in holdout.items()
               for i,_ in enumerate(models,1))
    assert len(conformers)==len(tags)==70
    out={
        "claim":"exact posthoc conformer-specific quantile noise radius, not unseen-population or measured dye uncertainty",
        "source":list(files),
        "heldout_conformers":70,
        "no_test_optimizing_sites_or_thresholds":True,
        "radius_units":"integer milliangstrom rationals; Angstrom float only presentation",
        "codes":{},
    }
    for k in (1,3,5,7,9):
        code=make_safe_code(train_a,train_b,k)
        assert check_safe_code(code,train_a,train_b)
        budgets=[]
        for e,t in (
            (0,0),(1,0),(2,0),
            (0,1),(1,1),(2,1),
            (0,2),(1,2),(2,2),
        ):
            try:b=risk_channels_for_budget(k,e,t)
            except ValueError:continue
            cert=certify_perturbation_radius(
                code,conformers,erased_channels=e,
                additional_wrong_channels=t)
            assert check_perturbation_certificate(code,conformers,cert)
            actual=[sorted_signed_margins(code,m,truth)[b]
                    for truth,m in conformers]
            tight=min(actual)
            if tight>0 and cert.strict_pair_distance_radius_milliangstrom!=tight:
                raise AssertionError("geometry radius does not equal exact minimum")
            origin=tags[actual.index(tight)]
            budgets.append({
                "erased_reporters":e,
                "arbitrary_bad_reporters":t,
                "geometry_risk_channels":b,
                "strict_pair_distance_noise_radius":as_ratio(cert.strict_pair_distance_radius_milliangstrom),
                "strict_each_CA_coordinate_displacement_radius":as_ratio(cert.strict_per_CA_site_radius_milliangstrom),
                "positive_radius_certified":cert.holds_for_all_submitted_conformers,
                "tightest_entry":origin[0],
                "tightest_model_index":origin[1],
                "tightest_signed_margin":as_ratio(tight),
                "checked_against_all_70_models":True,
            })
        out["codes"][str(k)]={
            "trained_mutation_safe_pairs":[list(e.endpoints) for e in code.reference.reporters],
            "budgeted_geometric_radii":budgets,
        }
    result=json.dumps(out,sort_keys=True,indent=2)+"\n"
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        args.out.write_text(result)
    print(result)


if __name__=="__main__":
    main()
