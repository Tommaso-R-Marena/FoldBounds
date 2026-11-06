"""Independently replay every calibration claim in an anytime assay certificate.

Deliberately does not call the certificate producer's grouping, proposal,
ranking, regret, feature, or quantile functions. This is a Python independent
checker, not a Lean kernel theorem, and source data must be authenticated by
its caller (e.g. pinned parse_s669). The group-exchangeability assumption
cannot be certified by a checker.
"""
from __future__ import annotations

from decimal import Decimal
from fractions import Fraction

from .assay_regret_certificates import AssayCertificateError, Mutation
from .anytime_assay_regret import AnytimeCertificate


def _independent_blocks(rows: list[Mutation], length: int | None):
    wildtypes = sorted({r.wt_seq for r in rows})
    if length is None:
        mapping = {s:s for s in wildtypes}
    else:
        if type(length) is not int or not 6 <= length <= 40:
            raise AssayCertificateError("invalid motif component length")
        # Independent brute-force pairwise shared k-mer graph construction.
        motif_sets = [set(w[i:i+length] for i in range(len(w)-length+1))
                      for w in wildtypes]
        unseen = set(range(len(wildtypes)))
        mapping = {}
        while unseen:
            start = min(unseen)
            frontier = [start]
            unseen.remove(start)
            for i in frontier:
                for j in sorted(unseen):
                    if motif_sets[i].intersection(motif_sets[j]):
                        frontier.append(j)
                        unseen.remove(j)
            rep = min(wildtypes[i] for i in frontier)
            for i in frontier:
                mapping[wildtypes[i]]=rep
    by_wt_cond: dict[tuple[str,Decimal,Decimal],list[Mutation]] = {}
    names=set()
    for r in rows:
        if r.name in names:
            raise AssayCertificateError("repeated mutation identifier")
        names.add(r.name)
        if not all(isinstance(x,Decimal) and x.is_finite()
                   for x in (r.ph,r.temp,r.pred)):
            raise AssayCertificateError("invalid physical conditions or predictor values")
        if r.ddg is None or not isinstance(r.ddg,Decimal) or not r.ddg.is_finite():
            raise AssayCertificateError("calibration checker requires all experimental outcomes")
        by_wt_cond.setdefault((r.wt_seq,r.ph,r.temp),[]).append(r)
    blocks={}
    for (wt,ph,temp),r in by_wt_cond.items():
        if len(r)>=2:
            blocks.setdefault(mapping[wt],[]).append(sorted(r,key=lambda x:x.name))
    return blocks


def _independent_group_score(panels, budgets, weights, lam):
    scores=[]
    for p in panels:
        scores_by_prediction=sorted(p,key=lambda x:(-x.pred,x.name))
        for b,w in zip(budgets,weights):
            k=max(1,(len(p)*b.numerator+b.denominator-1)//b.denominator)
            chosen=scores_by_prediction[:k]
            rest=scores_by_prediction[k:]
            best_known=max(Fraction(r.ddg) for r in chosen)
            best_all=max(Fraction(r.ddg) for r in p)
            lost=best_all-best_known
            if lost<0:
                raise AssertionError("true regret must be nonnegative")
            if rest:
                competitor=max(Fraction(r.pred) for r in rest)
                margin=max(Fraction(0),competitor-best_known)
            else:
                margin=Fraction(0)
            scores.append((lost-lam*margin)/w)
    if not scores:
        raise AssayCertificateError("empty eligible group")
    return max(scores)


def verify_anytime_certificate(cert: AnytimeCertificate, rows: list[Mutation]) -> dict:
    """Check all group membership, scores and exact rational quantile, fail closed."""
    if not isinstance(cert,AnytimeCertificate):
        raise AssayCertificateError("certificate object required")
    B,W,lam=cert.budgets,cert.weights,cert.lambda_margin
    if (not isinstance(B,tuple) or not B or
        any(not isinstance(x,Fraction) or x<=0 or x>1 for x in B) or
        tuple(sorted(set(B)))!=B or not isinstance(W,tuple) or len(W)!=len(B) or
        any(not isinstance(x,Fraction) or x<=0 for x in W) or
        not isinstance(lam,Fraction) or lam<0):
        raise AssayCertificateError("untrusted budget grid or weights")
    if not isinstance(cert.alpha,Fraction) or not 0<cert.alpha<1:
        raise AssayCertificateError("invalid alpha")
    by_group=_independent_blocks(rows,cert.motif_length)
    group_ids=cert.calibration_group_ids
    if (not isinstance(group_ids,tuple) or not group_ids or
        tuple(sorted(set(group_ids)))!=group_ids or
        not set(group_ids)<=set(by_group)):
        raise AssayCertificateError("invalid calibration group membership")
    expected=tuple(_independent_group_score(by_group[g],B,W,lam) for g in group_ids)
    if cert.exact_group_scores != expected:
        raise AssayCertificateError("submitted calibration scores differ from independent replay")
    m=len(expected)
    rank=((m+1)*(cert.alpha.denominator-cert.alpha.numerator)+cert.alpha.denominator-1)//cert.alpha.denominator
    q=sorted(expected)[rank-1] if rank<=m else None
    if q!=cert.cutoff:
        raise AssayCertificateError("cutoff is not the exact conformal order statistic")
    return {"verified":True,"calibration_groups":m,"budget_horizons":len(B),
            "per_horizon_alpha_correction":False,
            "abstains":q is None,
            "rank":rank,
            "checker":"independent Python finite data replay; assumes exchangeable groups"}


def verified_deployment_bound(cert: AnytimeCertificate, calibration_rows: list[Mutation],
                              unmeasured_candidates: tuple[Mutation, ...],
                              budget: Fraction, measured, *,
                              assay_error_epsilon: Fraction=Fraction(0)) -> Fraction | None:
    """Safety wrapper: replay all proof receipts and reject visible group leakage.

    This is the entry point for a real proposed panel. The formula-only helper
    is intentionally separate; a user-supplied cutoff is NOT a proof receipt.
    Homology exclusion beyond exact/shared-motif strings still requires a
    proper independently curated protein-family split.
    """
    from .anytime_assay_regret import bound_from_assayed
    verify_anytime_certificate(cert, calibration_rows)
    if not unmeasured_candidates:
        raise AssayCertificateError("no deployment candidates")
    test_seq=unmeasured_candidates[0].wt_seq
    verified_groups=_independent_blocks(calibration_rows,cert.motif_length)
    cal_seqs={p[0].wt_seq for group in cert.calibration_group_ids
              for p in verified_groups[group]}
    if test_seq in cal_seqs:
        raise AssayCertificateError("same wildtype appears in calibration")
    if cert.motif_length is not None:
        k=cert.motif_length
        motifs=set(test_seq[i:i+k] for i in range(len(test_seq)-k+1))
        for seq in cal_seqs:
            if any(seq[j:j+k] in motifs for j in range(len(seq)-k+1)):
                raise AssayCertificateError("shared exact motif with calibration protein")
    return bound_from_assayed(cert,unmeasured_candidates,budget,measured,
                             assay_error_epsilon=assay_error_epsilon)
