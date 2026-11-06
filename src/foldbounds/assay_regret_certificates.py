"""Finite-sample, group-conformal regret certificates for measured protein mutants.

Scope: published S669 experimental free-energy changes and model predictions;
NOT verified molecular dynamics, atomistic physics, or a deterministic bound.
No reference outcomes may be read to choose the candidate subset, and the
predictor must be fixed independently of the calibration/test labels.

Guarantee: if entire new protein (or motif-component) groups are exchangeable,
Pr(max_panel[best_real_measurement - best_measured_in_preselected_budget]
   <= calibrated_group_quantile) >= 1-alpha, marginally over new groups.
"""
from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from io import StringIO
from math import comb
from typing import Iterable

SOURCE_REPO = "Gonglab-THU/GeoStab"
SOURCE_COMMIT = "be8e6d079ea792396f3a4343a8de54273d3d3828"
SOURCE_PATH = "data/ddG/S669.csv"
SOURCE_GIT_BLOB = "f7a19c2ac0821371a947469231eb03d728ab884f"
MODEL = "GeoDDG-3D_dir"


class AssayCertificateError(ValueError):
    """Unverifiable source, split, experiment, or assurance input."""


@dataclass(frozen=True)
class Mutation:
    name: str
    wt_seq: str
    mut_seq: str
    ph: Decimal
    temp: Decimal
    ddg: Decimal | None  # None is permitted for future UNMEASURED proposals.
    pred: Decimal


@dataclass(frozen=True)
class Regret:
    max_panel_regret: Decimal
    mean_panel_regret: Decimal
    measured_candidates: int
    available_candidates: int
    panel_count: int


@dataclass(frozen=True)
class FoldResult:
    index: int
    calibration_groups: int
    test_groups: int
    covered_groups: int
    selected_budget: int
    available_candidates: int
    cutoff: Decimal | None  # None = INFINITE; never interpret as zero.


def _decimal(text: str, name: str) -> Decimal:
    if not isinstance(text, str) or not text.strip():
        raise AssayCertificateError(f"missing {name}")
    try:
        value = Decimal(text)
    except InvalidOperation as exc:
        raise AssayCertificateError(f"invalid {name}") from exc
    if not value.is_finite():
        raise AssayCertificateError(f"nonfinite {name}")
    return value


def git_blob_sha1(data: bytes) -> str:
    return hashlib.sha1(b"blob " + str(len(data)).encode() + b"\x00" + data).hexdigest()


def parse_s669(data: bytes, *, enforce_provenance: bool = True) -> list[Mutation]:
    """Validate exact CSV release and a real single-site substitution per row."""
    if not isinstance(data, bytes):
        raise AssayCertificateError("source must be bytes")
    if enforce_provenance and git_blob_sha1(data) != SOURCE_GIT_BLOB:
        raise AssayCertificateError("published S669 Git blob provenance mismatch")
    try:
        reader = csv.DictReader(StringIO(data.decode("utf-8")))
    except UnicodeDecodeError as exc:
        raise AssayCertificateError("malformed UTF-8") from exc
    needed = {"name", "wt_seq", "mut_seq", "pH", "temp", "ddG", MODEL}
    if reader.fieldnames is None or not needed <= set(reader.fieldnames):
        raise AssayCertificateError("missing S669 fields")
    out: list[Mutation] = []
    ids: set[str] = set()
    for row in reader:
        if None in row or any(v is None for v in row.values()):
            raise AssayCertificateError("malformed source row")
        name, wt, mt = row["name"], row["wt_seq"], row["mut_seq"]
        if not name or name in ids:
            raise AssayCertificateError("missing or duplicate mutation identifier")
        ids.add(name)
        if not wt or not mt or len(wt) != len(mt) or sum(a != b for a, b in zip(wt, mt)) != 1:
            raise AssayCertificateError("source row must be one single-site substitution")
        out.append(Mutation(name, wt, mt, _decimal(row["pH"],"pH"),
                            _decimal(row["temp"],"temp"),
                            _decimal(row["ddG"],"ddG"),
                            _decimal(row[MODEL],MODEL)))
    if enforce_provenance and (len(out) != 669 or len({r.wt_seq for r in out}) != 93):
        raise AssayCertificateError("pinned dataset dimensions unexpected")
    if not out:
        raise AssayCertificateError("empty S669 input")
    return out


def assay_panels(rows: Iterable[Mutation]) -> dict[str, list[tuple[tuple[Decimal, Decimal], tuple[Mutation, ...]]]]:
    """Group same wildtype and actual stated assay conditions; discard n=1 panels.

    Different pH/temperature entries must NEVER share an optimization maximum.
    Each protein can have several panels. Entire proteins remain indivisible
    under statistical calibration/testing, regardless of panel multiplicity.
    """
    by_protein: dict[str, dict[tuple[Decimal,Decimal], list[Mutation]]] = defaultdict(lambda: defaultdict(list))
    seen: set[str] = set()
    for row in rows:
        if row.name in seen:
            raise AssayCertificateError("duplicate variant identifier")
        seen.add(row.name)
        if row.ddg is None:
            raise AssayCertificateError("retrospective assay regret requires observed ddG")
        if not all(x.is_finite() for x in (row.ph,row.temp,row.ddg,row.pred)):
            raise AssayCertificateError("nonfinite numeric source")
        by_protein[row.wt_seq][(row.ph,row.temp)].append(row)
    return {
        wt: [(cond, tuple(sorted(variants,key=lambda x:x.name)))
             for cond, variants in sorted(panels.items()) if len(variants) >= 2]
        for wt,panels in sorted(by_protein.items())
        if any(len(variants) >= 2 for variants in panels.values())
    }


def motif_components(wt_sequences: Iterable[str], *, motif_length: int = 12) -> dict[str,str]:
    """Block transitive exact k-mer overlaps (NOT proof of homology)."""
    if type(motif_length) is not int or not 6 <= motif_length <= 40:
        raise AssayCertificateError("motif_length must be an integer 6..40")
    sequences=sorted(set(wt_sequences))
    if not sequences or any(not isinstance(s,str) or not s for s in sequences):
        raise AssayCertificateError("empty or invalid wildtype sequences")
    parent=list(range(len(sequences)))

    def find(i):
        while i!=parent[i]:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i

    def union(i,j):
        a,b=find(i),find(j)
        if a!=b:
            parent[b]=a
    owner={}
    for i,seq in enumerate(sequences):
        for motif in set(seq[j:j+motif_length]
                         for j in range(max(0,len(seq)-motif_length+1))):
            if motif in owner:
                union(i,owner[motif])
            else:
                owner[motif]=i
    components=defaultdict(list)
    for i,s in enumerate(sequences):
        components[find(i)].append(s)
    return {s:min(group) for group in components.values() for s in group}


def rank_select(panel: Iterable[Mutation], budget: Fraction, *, selector: str="model") -> tuple[Mutation,...]:
    """Choose ceil(budget*n) candidates without reading experimental ddG."""
    if not isinstance(budget, Fraction) or not 0 < budget <= 1:
        raise AssayCertificateError("budget must be exact Fraction in (0,1]")
    variants=list(panel)
    if not variants:
        raise AssayCertificateError("cannot select from an empty panel")
    if len({r.name for r in variants}) != len(variants):
        raise AssayCertificateError("duplicate candidate")
    if selector=="model":
        variants.sort(key=lambda r:(-r.pred,r.name))
    elif selector=="identifier_baseline":
        variants.sort(key=lambda r:r.name)
    else:
        raise AssayCertificateError("unrecognized, potentially outcome-optimized selector")
    k=max(1,(len(variants)*budget.numerator+budget.denominator-1)//budget.denominator)
    return tuple(variants[:k])


def regret_from_panels(panels, budget: Fraction, *, selector: str="model") -> Regret:
    """Maximum decision regret across a complete protein's matched assays."""
    errors=[]; spent=available=0
    for _,panel in panels:
        if len(panel)<2 or any(x.ddg is None for x in panel):
            raise AssayCertificateError("regret needs all truth labels for each eligible panel")
        selected=rank_select(panel,budget,selector=selector)
        best=max(x.ddg for x in panel)
        chosen=max(x.ddg for x in selected)
        error=best-chosen
        if error<0:
            raise AssertionError("regret cannot be negative")
        errors.append(error)
        spent+=len(selected)
        available+=len(panel)
    if not errors:
        raise AssayCertificateError("no eligible assay panel")
    return Regret(max(errors),sum(errors,Decimal(0))/len(errors),spent,available,len(errors))


def expected_uniform_random_panel_regret(panel, budget: Fraction) -> Fraction:
    """EXACT rational expected regret over uniformly sampled k-subsets.

    For sorted y_1<=...<=y_n, maximum's rank is j with probability
    C(j-1,k-1)/C(n,k). Uses exact Fraction(Decimal), not floating
    approximation of experimental decimals or probabilities.
    """
    variants=list(panel)
    if not variants or any(x.ddg is None for x in variants):
        raise AssayCertificateError("random oracle requires all labels")
    if not isinstance(budget,Fraction) or not 0<budget<=1:
        raise AssayCertificateError("invalid budget")
    k=max(1,(len(variants)*budget.numerator+budget.denominator-1)//budget.denominator)
    values=sorted(Fraction(x.ddg) for x in variants)
    z=comb(len(values),k)
    expected=values[-1]-sum((values[j-1]*Fraction(comb(j-1,k-1),z)
                             for j in range(k,len(values)+1)),Fraction(0))
    if expected<0:
        raise AssertionError("random selection cannot exceed the true optimum")
    return expected


def protein_random_expected_regret(panels,budget:Fraction) -> Fraction:
    # Per-protein mean of within-condition panel regrets. An expectation of
    # max panel regrets would require a joint distribution of random subsets.
    if not panels:
        raise AssayCertificateError("empty panels")
    return sum((expected_uniform_random_panel_regret(panel,budget)
                for _,panel in panels),Fraction(0))/len(panels)


def conformal_rank_cutoff(calibration_scores: Iterable[Decimal],
                          alpha:Fraction=Fraction(1,10)) -> Decimal | None:
    """Rank order statistic with ∞ sentinel when finite-sample target impossible."""
    if not isinstance(alpha,Fraction) or not 0<alpha<1:
        raise AssayCertificateError("alpha must be an exact rational strictly between 0 and 1")
    a=list(calibration_scores)
    if not a or any(not isinstance(x,Decimal) or not x.is_finite() or x<0 for x in a):
        raise AssayCertificateError("nonempty nonnegative finite Decimal scores required")
    m=len(a)
    k=((m+1)*(alpha.denominator-alpha.numerator)+alpha.denominator-1)//alpha.denominator
    return None if k>m else sorted(a)[k-1]


def _shuffle(groups: Iterable[str],seed:int) -> list[str]:
    if type(seed) is not int or not 0<seed<2**32:
        raise AssayCertificateError("seed must be a nonzero 32-bit integer")
    a=sorted(groups)
    state=seed
    for i in range(len(a)-1,0,-1):
        state^=(state<<13)&0xffffffff
        state^=state>>17
        state^=(state<<5)&0xffffffff
        state&=0xffffffff
        j=(state*(i+1))//2**32
        a[i],a[j]=a[j],a[i]
    return a


def block_scores(rows:list[Mutation], budget:Fraction, *, motif_length:int|None=None,
                 selector:str="model") -> tuple[dict[str,Regret], dict[str,list[str]]]:
    """Construct complete protein or motif-cluster maximum-regret scores."""
    panels=assay_panels(rows)
    if motif_length is None:
        mapping={s:s for s in panels}
    else:
        # Construct components on ALL wildtype sequences *before* selecting
        # eligible panels, avoiding a transitive bridge through n=1 proteins.
        mapping=motif_components((r.wt_seq for r in rows),motif_length=motif_length)
    blocks=defaultdict(list)
    for protein in panels:
        blocks[mapping[protein]].append(protein)
    scores={}
    for block,members in blocks.items():
        r=[regret_from_panels(panels[s],budget,selector=selector) for s in members]
        numerator=sum((v.mean_panel_regret*v.panel_count for v in r),Decimal(0))
        n=sum(v.panel_count for v in r)
        scores[block]=Regret(max(v.max_panel_regret for v in r),numerator/n,
                             sum(v.measured_candidates for v in r),
                             sum(v.available_candidates for v in r),n)
    return scores,dict(blocks)


def evaluate(rows:list[Mutation],budget=Fraction(1,5),*,motif_length:int|None=None,
             alpha=Fraction(1,10),selector="model",seed=20261105) -> tuple[FoldResult,...]:
    """Five held-out folds; every test fold uses OTHER FOUR as calibration.

    The published model is never refitted, so no training fold is needed.
    Test outcomes NEVER enter that test fold's calibration cutoff, although
    they may serve as calibration for a different held-out fold. Full
    experimental evaluation is not an independent prospective sample.
    """
    score_map,_=block_scores(rows,budget,motif_length=motif_length,selector=selector)
    keys=sorted(score_map)
    if len(keys)<5:
        raise AssayCertificateError("too few valid groups to evaluate")
    assigned={g:i%5 for i,g in enumerate(_shuffle(keys,seed))}
    outputs=[]
    for f in range(5):
        test=[g for g in keys if assigned[g]==f]
        cal=[g for g in keys if assigned[g]!=f]
        if not test or not cal or set(test)&set(cal):
            raise AssayCertificateError("invalid calibration/test split")
        cutoff=conformal_rank_cutoff((score_map[g].max_panel_regret for g in cal),alpha)
        covered=sum(cutoff is None or score_map[g].max_panel_regret<=cutoff for g in test)
        outputs.append(FoldResult(f,len(cal),len(test),covered,
                                  sum(score_map[g].measured_candidates for g in test),
                                  sum(score_map[g].available_candidates for g in test),cutoff))
    return tuple(outputs)


def summarize(folds:tuple[FoldResult,...]) -> dict:
    n=sum(f.test_groups for f in folds)
    total=sum(f.selected_budget for f in folds)
    available=sum(f.available_candidates for f in folds)
    if not n or not total:
        raise AssayCertificateError("empty evaluation")
    if any(f.cutoff is None for f in folds):
        q=None
    else:
        q=sum((f.cutoff*f.test_groups for f in folds),Decimal(0))/n
    return {
        "groups_covered":sum(f.covered_groups for f in folds),
        "groups_total":n,
        "nominal_regret_cutoff_mean_kcal_mol":None if q is None else str(q),
        "budget_assays":total,
        "available_assays":available,
        "folds":[{"fold":f.index,"calibration":f.calibration_groups,
                  "test":f.test_groups,"covered":f.covered_groups,
                  "regret_cutoff":None if f.cutoff is None else str(f.cutoff)}
                 for f in folds],
    }


def certify_unmeasured_plan(variants: Iterable[Mutation], budget:Fraction,
                            q:Decimal|None) -> dict:
    """Create concrete forward assay list without consulting any truth labels.

    A *valid* q must have been independently obtained by group calibration;
    this helper does NOT verify the provenance of a user-supplied q.
    """
    panel=list(variants)
    if any(x.ddg is not None for x in panel):
        raise AssayCertificateError("deployment plan must not contain measured outcomes")
    if len({(x.wt_seq,x.ph,x.temp) for x in panel})!=1:
        raise AssayCertificateError("deployment plan must contain one matched-assay panel")
    if q is not None and (not isinstance(q,Decimal) or not q.is_finite() or q<0):
        raise AssayCertificateError("cutoff must be nonnegative finite Decimal or abstention")
    selected=rank_select(panel,budget)
    return {"candidate_count":len(panel),"measure_count":len(selected),
            "mutation_ids":[x.name for x in selected],
            "regret_upper_kcal_mol":None if q is None else str(q),
            "status":"abstain" if q is None else "marginal-group-assurance-only",
            "note":"Only marginal for a randomly exchangeable new whole group; not a bound for a fixed protein."}
