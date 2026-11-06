"""Protein-disjoint group-conformal assurance on real experimental stability data.

Physical requirement: ΔΔG(A→B) = -ΔΔG(B→A) for the *same conditions*.
Statistical requirement: entire new protein groups exchangeable with the
calibration groups, and any model selection isolated to training groups.
No guarantee for a particular fixed protein or outside that distribution.

The measured target and published 3D-model outputs are in the S669 CSV.
This code never trains on evaluation protein groups.
"""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path

SOURCE_REPO = "Gonglab-THU/GeoStab"
SOURCE_COMMIT = "be8e6d079ea792396f3a4343a8de54273d3d3828"
SOURCE_PATH = "data/ddG/S669.csv"
SOURCE_BLOB_SHA1 = "f7a19c2ac0821371a947469231eb03d728ab884f"
MODELS = ("GeoDDG-3D", "GeoDDG-Seq", "ACDC-NN", "DDMut", "Dynamut")


class AssuranceError(ValueError):
    pass


def git_blob_sha1(raw: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()


def _number(value, field):
    try:
        if not isinstance(value, str) or not value.strip():
            raise AssuranceError(f"missing {field}")
        d = Decimal(value)
        if not d.is_finite():
            raise AssuranceError(f"nonfinite {field}")
        return d
    except InvalidOperation as exc:
        raise AssuranceError(f"invalid {field}") from exc


@dataclass(frozen=True)
class Measurement:
    name: str
    wt: str
    mutant: str
    ph: Decimal
    temperature: Decimal
    ddg: Decimal
    estimates: dict[str, Decimal]

    def predict(self, spec):
        if spec.endswith("_dir"):
            return self.estimates[spec]
        if spec.endswith("_sym"):
            name = spec[:-4]
            return (self.estimates[name+"_dir"] -
                    self.estimates[name+"_inv"])/2
        raise AssuranceError("unsupported prediction specification")


def load_s669(path, *, verify_provenance=True):
    raw = Path(path).read_bytes()
    if verify_provenance and git_blob_sha1(raw) != SOURCE_BLOB_SHA1:
        raise AssuranceError("source does not match pinned published Git blob")
    try:
        parsed = csv.DictReader(raw.decode("utf-8").splitlines())
    except UnicodeDecodeError as exc:
        raise AssuranceError("invalid UTF-8") from exc
    required = {"name", "wt_seq", "mut_seq", "ddG", "pH", "temp"}
    required.update(m+"_"+s for m in MODELS for s in ("dir","inv"))
    if parsed.fieldnames is None or not required.issubset(parsed.fieldnames):
        raise AssuranceError("missing required S669 source fields")
    records = []
    names = set()
    for row in parsed:
        if None in row:
            raise AssuranceError("malformed CSV")
        name, wt, mut = row["name"], row["wt_seq"], row["mut_seq"]
        if not name or name in names:
            raise AssuranceError("duplicate or missing variant identifier")
        names.add(name)
        if (not wt or not mut or len(wt) != len(mut)
                or sum(x != y for x,y in zip(wt,mut)) != 1):
            raise AssuranceError("non single-site mutation or invalid sequence")
        estimates = {
            m+"_"+s: _number(row[m+"_"+s],m+"_"+s)
            for m in MODELS for s in ("dir","inv")
        }
        records.append(Measurement(
            name, wt, mut, _number(row["pH"],"pH"),
            _number(row["temp"],"temp"), _number(row["ddG"],"ddG"),
            estimates
        ))
    if not records:
        raise AssuranceError("empty benchmark")
    if verify_provenance and (len(records) != 669 or len({r.wt for r in records}) != 93):
        raise AssuranceError("unexpected S669 release size")
    return records


def antisymmetric_projection(forward: Decimal, reverse: Decimal):
    """Euclidean projection onto {(t,-t)}; not a forward-error theorem."""
    t = (forward-reverse)/2
    return t, -t


def reciprocity_pythagoras(f: Decimal, r: Decimal, truth: Decimal):
    """Exact identity for a thermodynamically reversible true target.

    (f-y)^2+(r+y)^2=2*((f-r)/2-y)^2 + (f+r)^2/2.
    """
    t,_ = antisymmetric_projection(f,r)
    return (f-truth)**2+(r+truth)**2 == 2*(t-truth)**2+(f+r)**2/2


def _xorshift_shuffle(items, seed):
    if type(seed) is not int or not 0 < seed < 2**32:
        raise AssuranceError("seed must be nonzero uint32")
    state = seed
    items = sorted(items)
    for i in range(len(items)-1,0,-1):
        state ^= (state<<13)&0xffffffff
        state ^= state>>17
        state ^= (state<<5)&0xffffffff
        state &= 0xffffffff
        j = (state*(i+1))//2**32
        items[i],items[j] = items[j],items[i]
    return items


def protein_folds(records, seed=20261105):
    groups = {}
    for r in records:
        groups.setdefault(r.wt,[]).append(r)
    if len(groups) < 5:
        raise AssuranceError("too few protein groups")
    assignment={g:i%5 for i,g in enumerate(_xorshift_shuffle(groups,seed))}
    return assignment, groups


def choose_predictor(train, groups, *, enforce_reciprocity=True):
    """Select solely by unweighted protein-macro MAE in TRAIN groups."""
    specs=(m+"_sym" for m in MODELS) if enforce_reciprocity else (
        m+"_"+mode for m in MODELS for mode in ("dir","sym"))
    if not train:
        raise AssuranceError("empty training proteins")
    scores=[]
    for spec in specs:
        errors=[
            sum((abs(x.ddg-x.predict(spec)) for x in groups[g]),Decimal(0))/len(groups[g])
            for g in train
        ]
        scores.append((sum(errors,Decimal(0))/len(train),spec))
    return min(scores)[1]


def conformal_cutoff(cal, groups, spec, alpha=Fraction(1,10), *, unit="protein"):
    """Exact order statistic; None denotes an INFINITE/abstaining interval.

    For unit='protein', each score is a maximum absolute residual over ALL
    mutations of a calibration protein. Under exchangeability of whole
    protein groups, the entire new-protein panel is covered >=1-alpha
    marginally. No individual conditional guarantee and no physics proof.
    unit='mutation' is a negative control, without whole-protein validity.
    """
    if not isinstance(alpha,Fraction) or not 0 < alpha < 1:
        raise AssuranceError("alpha must be a Fraction between 0 and 1")
    if not cal or len(set(cal)) != len(cal):
        raise AssuranceError("calibration groups missing/duplicated")
    if unit=="protein":
        values = [max(abs(x.ddg-x.predict(spec)) for x in groups[g]) for g in cal]
    elif unit=="mutation":
        values = [abs(x.ddg-x.predict(spec)) for g in cal for x in groups[g]]
    else:
        raise AssuranceError("unknown calibration unit")
    values.sort()
    m=len(values)
    k=((m+1)*(alpha.denominator-alpha.numerator)+alpha.denominator-1)//alpha.denominator
    return None if k>m else values[k-1]


def evaluate(records, *, seed=20261105, alpha=Fraction(1,10),
             strategy="fixed", calibration_unit="protein"):
    """Cross-fitting: test fold i, calibration fold i+1, remaining train."""
    if strategy not in ("fixed","selected_phys","selected_any"):
        raise AssuranceError("invalid predeclared strategy")
    assignment,groups=protein_folds(records,seed)
    folds=[]
    for f in range(5):
        test=[g for g in groups if assignment[g]==f]
        cal=[g for g in groups if assignment[g]==(f+1)%5]
        train=[g for g in groups if assignment[g] not in (f,(f+1)%5)]
        if set(test)&set(cal) or set(train)&(set(test)|set(cal)):
            raise AssuranceError("data leakage between protein groups")
        spec=("GeoDDG-3D_dir" if strategy=="fixed" else
              choose_predictor(train,groups,enforce_reciprocity=(strategy=="selected_phys")))
        cutoff=conformal_cutoff(cal,groups,spec,alpha,unit=calibration_unit)
        group_hits=row_hits=sign_called=sign_right=rows=0
        for g in test:
            hits=[]
            for x in groups[g]:
                rows+=1
                pred=x.predict(spec)
                inside=cutoff is None or abs(pred-x.ddg)<=cutoff
                row_hits+=int(inside)
                hits.append(inside)
                if cutoff is not None and abs(pred)>cutoff:
                    sign_called+=1
                    sign_right+=int(x.ddg!=0 and ((pred>0)==(x.ddg>0)))
            group_hits+=int(all(hits))
        folds.append({
            "fold":f,"predictor":spec,"cal_proteins":len(cal),
            "test_proteins":len(test),"test_mutations":rows,
            "cutoff":None if cutoff is None else str(cutoff),
            "proteins_covered":group_hits,"mutations_covered":row_hits,
            "sign_decisions":sign_called,"correct_sign_decisions":sign_right,
        })
    return folds


def summarize(folds):
    proteins=sum(x["test_proteins"] for x in folds)
    rows=sum(x["test_mutations"] for x in folds)
    total_width=Decimal(0)
    for x in folds:
        if x["cutoff"] is None:
            width="Infinity"
            break
        total_width+=2*Decimal(x["cutoff"])*x["test_mutations"]
    else:
        width=str(total_width/rows)
    return {
        "proteins_covered":sum(x["proteins_covered"] for x in folds),
        "proteins_total":proteins,
        "mutations_covered":sum(x["mutations_covered"] for x in folds),
        "mutations_total":rows,
        "mean_interval_width_kcal_mol":width,
        "sign_decisions":sum(x["sign_decisions"] for x in folds),
        "correct_sign_decisions":sum(x["correct_sign_decisions"] for x in folds),
        "folds":folds,
    }
