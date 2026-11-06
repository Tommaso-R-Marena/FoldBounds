"""Sequence-similarity-blocked stress test for real S669 stability assurance.

Two full wildtype sequences are linked when they share >=k consecutive
identical residues. Connected components are indivisible test/cal/train
units. Shared motifs may come from affinity tags rather than ancestry:
THIS IS A CONSERVATIVE MOTIF STRESS TEST, NOT TRUE HOMOLOGY INFERENCE.

A group-max conformal rank guarantee would require exchangeability of
entire motif components, not merely exchangeability of mutations or
proteins. See the experimental assurance mathematical assumptions.
"""
from __future__ import annotations

from decimal import Decimal
from fractions import Fraction

from .experimental_reciprocal_conformal import (
    AssuranceError, _xorshift_shuffle, MODELS,
    choose_predictor, conformal_cutoff,
)


def motif_components(sequences, *, k: int=12) -> dict[str,str]:
    """Return sequence -> lexicographically least representative.

    Transitive closure prevents leakage through chains A~B, B~C even when
    A and C do not share any individual long identical motif.
    """
    if type(k) is not int or not 6 <= k <= 40:
        raise AssuranceError("k must be integer 6..40")
    seqs=sorted(set(sequences))
    if not seqs or any(not isinstance(s,str) or not s for s in seqs):
        raise AssuranceError("invalid sequence collection")
    parent=list(range(len(seqs)))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]]
            i=parent[i]
        return i
    def union(i,j):
        a,b=find(i),find(j)
        if a!=b:
            parent[b]=a
    motif_owner={}
    for i,seq in enumerate(seqs):
        if len(seq)<k:
            continue
        kmers={seq[j:j+k] for j in range(len(seq)-k+1)}
        for motif in kmers:
            if motif in motif_owner:
                union(i,motif_owner[motif])
            else:
                motif_owner[motif]=i
    components={}
    for i,seq in enumerate(seqs):
        components.setdefault(find(i),[]).append(seq)
    return {seq:min(component) for component in components.values()
            for seq in component}


def motif_blocks(records, *, k=12, seed=20261105):
    representatives=motif_components((r.wt for r in records),k=k)
    blocks={}
    for r in records:
        blocks.setdefault(representatives[r.wt],[]).append(r)
    if len(blocks)<5:
        raise AssuranceError("too few motif components for 5-fold calibration")
    order=_xorshift_shuffle(blocks,seed)
    folds={g:i%5 for i,g in enumerate(order)}
    return folds,blocks,representatives


def evaluate_motif_blocked(records, *, k=12, seed=20261105,
                           alpha=Fraction(1,10), strategy="fixed",
                           calibration="component"):
    """Train/calibrate/evaluate by indivisible motif components, never mutations."""
    if strategy not in ("fixed","selected_phys"):
        raise AssuranceError("unknown predeclared predictor strategy")
    if calibration not in ("component","mutation"):
        raise AssuranceError("calibration must be component or negative-control mutation")
    split,blocks,reps=motif_blocks(records,k=k,seed=seed)
    out=[]
    for f in range(5):
        test=[g for g in blocks if split[g]==f]
        cal=[g for g in blocks if split[g]==(f+1)%5]
        train=[g for g in blocks if split[g] not in (f,(f+1)%5)]
        if set(test)&set(cal) or set(train)&(set(test)|set(cal)):
            raise AssuranceError("motif component leakage between experiment partitions")
        prediction=("GeoDDG-3D_dir" if strategy=="fixed"
                    else choose_predictor(train,blocks,enforce_reciprocity=True))
        cutoff=conformal_cutoff(cal,blocks,prediction,alpha,
                                 unit="protein" if calibration=="component" else "mutation")
        hit_components=hit_rows=sign_n=sign_right=n_rows=0
        for group in test:
            successes=[]
            for obs in blocks[group]:
                n_rows+=1
                pred=obs.predict(prediction)
                covered=(cutoff is None or abs(obs.ddg-pred)<=cutoff)
                successes.append(covered)
                hit_rows+=int(covered)
                if cutoff is not None and abs(pred)>cutoff:
                    sign_n+=1
                    sign_right+=int(obs.ddg!=0 and (pred>0)==(obs.ddg>0))
            hit_components+=int(all(successes))
        out.append({
            "fold":f, "predictor":prediction, "cal_components":len(cal),
            "test_components":len(test), "test_mutations":n_rows,
            "cutoff":None if cutoff is None else str(cutoff),
            "components_covered":hit_components, "mutations_covered":hit_rows,
            "sign_decisions":sign_n, "correct_sign_decisions":sign_right,
        })
    return out


def summarize_motif(folds):
    n_blocks=sum(f["test_components"] for f in folds)
    n_rows=sum(f["test_mutations"] for f in folds)
    if not n_blocks or not n_rows:
        raise AssuranceError("empty evaluation")
    if any(f["cutoff"] is None for f in folds):
        width="Infinity"
    else:
        width=str(sum(2*Decimal(f["cutoff"])*f["test_mutations"]
                      for f in folds)/n_rows)
    return {
        "components_covered":sum(f["components_covered"] for f in folds),
        "components_total":n_blocks,
        "mutations_covered":sum(f["mutations_covered"] for f in folds),
        "mutations_total":n_rows,
        "mean_interval_width_kcal_mol":width,
        "sign_decisions":sum(f["sign_decisions"] for f in folds),
        "correct_sign_decisions":sum(f["correct_sign_decisions"] for f in folds),
        "folds":folds,
    }
