"""Proof-carrying, multi-budget protein assay regret assurance (Python prototype).

A *single* conformal score maximizes over all predeclared assay budgets and
all condition-matched panels in a protein or label-free motif block. One
rank calibration then covers every budget simultaneously for an exchangeable
new complete block. Adaptive stopping after measured selected candidates
inherits the SAME marginal coverage event, without a second independence
assumption. This is not a conditional guarantee for a fixed protein.

Everything is exact rational arithmetic, including the quantile. No assumption
on the validity of the molecular predictor is needed for rank coverage.
"""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from typing import Iterable, Mapping

from .assay_regret_certificates import (
    AssayCertificateError, Mutation, assay_panels, motif_components,
    rank_select,
)

DEFAULT_BUDGETS = (Fraction(1, 20), Fraction(1, 10),
                   Fraction(1, 5), Fraction(1, 2))


def _frac(value: Decimal) -> Fraction:
    if not isinstance(value, Decimal) or not value.is_finite():
        raise AssayCertificateError("expected exact finite Decimal")
    return Fraction(value)


def check_grid(budgets: tuple[Fraction, ...], weights: tuple[Fraction, ...],
               lam: Fraction) -> None:
    if (not isinstance(budgets, tuple) or not budgets or
        any(not isinstance(b, Fraction) or not 0 < b <= 1 for b in budgets) or
        tuple(sorted(set(budgets))) != budgets):
        raise AssayCertificateError("budgets must be a strictly increasing nonempty tuple of Fractions in (0,1]")
    if (not isinstance(weights, tuple) or len(weights) != len(budgets) or
        any(not isinstance(w, Fraction) or w <= 0 for w in weights)):
        raise AssayCertificateError("positive exact-rational weights, one per budget, are required")
    if not isinstance(lam, Fraction) or lam < 0:
        raise AssayCertificateError("nonnegative exact-rational lambda required")


@dataclass(frozen=True)
class Horizon:
    budget: Fraction
    selected: tuple[str, ...]
    total: int
    assayed: int
    observed_margin: Fraction
    actual_regret: Fraction | None
    wildtype: str
    condition: tuple[Decimal, Decimal]


@dataclass(frozen=True)
class AnytimeCertificate:
    """Explicit calibration witness; q=None represents infinite abstention."""
    budgets: tuple[Fraction, ...]
    weights: tuple[Fraction, ...]
    lambda_margin: Fraction
    alpha: Fraction
    motif_length: int | None
    calibration_group_ids: tuple[str, ...]
    exact_group_scores: tuple[Fraction, ...]
    cutoff: Fraction | None


def evidence_for_panel(panel: tuple[Mutation, ...], budget: Fraction,
                       *, measured: Mapping[str, Decimal] | None = None,
                       strict_deployment: bool = False) -> Horizon:
    """Replay a fixed prediction-only top-k policy with selected assay data.

    In a retrospective panel, measured=None and all experimental labels are
    required; regret is recomputed. In deployment, all labels *on the input
    candidates* MUST be None and only selected assays may be supplied through
    the separate measured mapping; it is impossible to access withheld truth.
    """
    if not panel or len({p.name for p in panel}) != len(panel):
        raise AssayCertificateError("empty or duplicate panel")
    if len({(p.wt_seq, p.ph, p.temp) for p in panel}) != 1:
        raise AssayCertificateError("mixing assay condition, protein, or physical state")
    if any(not p.pred.is_finite() for p in panel):
        raise AssayCertificateError("nonfinite model prediction")
    if strict_deployment:
        if measured is None or any(p.ddg is not None for p in panel):
            raise AssayCertificateError("deployment needs truth-free candidates and separate selected assays")
    elif measured is not None:
        raise AssayCertificateError("explicit measured map only allowed in strict deployment")
    selected = rank_select(panel, budget)
    selected_ids = tuple(p.name for p in selected)
    if strict_deployment:
        if set(measured) != set(selected_ids):
            raise AssayCertificateError("assay labels cannot be missing, extra, or withdrawn from an unselected mutant")
        values = [_frac(measured[name]) for name in selected_ids]
    else:
        if any(p.ddg is None or not p.ddg.is_finite() for p in panel):
            raise AssayCertificateError("retrospective panel requires ALL observed outcomes")
        values = [_frac(p.ddg) for p in selected]
    best_selected = max(values)
    excluded = [p for p in panel if p.name not in set(selected_ids)]
    unmeasured_pred_max = max((_frac(p.pred) for p in excluded), default=None)
    margin = max(Fraction(0), unmeasured_pred_max - best_selected) if excluded else Fraction(0)
    true_regret = None
    if not strict_deployment:
        true_regret = max(_frac(p.ddg) for p in panel) - best_selected
        if true_regret < 0:
            raise AssertionError("definition of maximum regret violated")
    return Horizon(
        budget, selected_ids, len(panel), len(selected), margin, true_regret,
        panel[0].wt_seq, (panel[0].ph, panel[0].temp),
    )


def assay_group_evidence(rows: list[Mutation], budgets: tuple[Fraction, ...],
                         *, motif_length: int | None = None) -> dict[str, tuple[Horizon, ...]]:
    """Group every eligible matched-condition panel, then every motif block.

    Motif connected components are inferred on ALL wildtypes (including
    singleton-panel intermediates) to stop transitive train/test leakage.
    """
    check_grid(budgets, tuple(Fraction(1) for _ in budgets), Fraction(0))
    panels = assay_panels(rows)
    if motif_length is None:
        groups = {s: s for s in panels}
    else:
        groups = motif_components((r.wt_seq for r in rows), motif_length=motif_length)
    by_group: dict[str, list[Horizon]] = {}
    for protein, entries in panels.items():
        key = groups[protein]
        for _, p in entries:
            for budget in budgets:
                by_group.setdefault(key, []).append(evidence_for_panel(p, budget))
    return {g: tuple(ev) for g, ev in sorted(by_group.items())}


def group_score(evidence: tuple[Horizon, ...], budgets: tuple[Fraction, ...],
                weights: tuple[Fraction, ...], lam: Fraction) -> Fraction:
    check_grid(budgets, weights, lam)
    if not evidence:
        raise AssayCertificateError("cannot certify an empty experimental protein group")
    wt = dict(zip(budgets, weights))
    terms = []
    for item in evidence:
        if item.budget not in wt or item.actual_regret is None:
            raise AssayCertificateError("missing complete calibration budget panel or truth")
        terms.append((item.actual_regret - lam * item.observed_margin) / wt[item.budget])
    return max(terms)


def signed_conformal_cutoff(scores: Iterable[Fraction], alpha: Fraction) -> Fraction | None:
    """Exact signed finite-sample quantile; no unjustified clipping at zero."""
    if not isinstance(alpha, Fraction) or not 0 < alpha < 1:
        raise AssayCertificateError("alpha must be exact Fraction strictly between zero and one")
    values = list(scores)
    if not values or any(not isinstance(v, Fraction) for v in values):
        raise AssayCertificateError("missing or nonrational calibration score")
    k = ((len(values)+1)*(alpha.denominator-alpha.numerator)+alpha.denominator-1)//alpha.denominator
    return sorted(values)[k-1] if k <= len(values) else None


def propose_certificate(rows: list[Mutation], calibration_groups: Iterable[str],
                        *, budgets=DEFAULT_BUDGETS, weights=None,
                        lambda_margin=Fraction(1), alpha=Fraction(1,10),
                        motif_length=None) -> AnytimeCertificate:
    """Untrusted producer: the separate verifier must recompute these scores."""
    if weights is None:
        weights = tuple(Fraction(1) for _ in budgets)
    check_grid(budgets, weights, lambda_margin)
    evidence = assay_group_evidence(rows, budgets, motif_length=motif_length)
    groups = tuple(sorted(calibration_groups))
    if (not groups or len(groups) != len(set(groups)) or
        not set(groups) <= set(evidence)):
        raise AssayCertificateError("missing/duplicate/unknown calibration group")
    scores = tuple(group_score(evidence[g],budgets,weights,lambda_margin) for g in groups)
    return AnytimeCertificate(budgets,weights,lambda_margin,alpha,motif_length,
                               groups,scores,signed_conformal_cutoff(scores,alpha))


def robust_noise_radius(cert: AnytimeCertificate, epsilon: Fraction) -> Fraction:
    """Exact worst-case uniform measurement perturbation of one-sided scores.

    If each *calibration* candidate assay may differ from latent truth by
    epsilon, |regret_obs-regret_true|<=2epsilon and
    |margin_obs-margin_true|<=epsilon for each budget. Hence the entire
    signed score changes by <= max_b((2+lambda)*epsilon/weight_b).
    """
    check_grid(cert.budgets, cert.weights, cert.lambda_margin)
    if not isinstance(epsilon, Fraction) or epsilon < 0:
        raise AssayCertificateError("epsilon must be a nonnegative exact Fraction")
    return max((Fraction(2)+cert.lambda_margin)*epsilon/w for w in cert.weights)


def bound_from_assayed(cert: AnytimeCertificate, panel: tuple[Mutation, ...],
                       budget: Fraction, measured: Mapping[str, Decimal],
                       *, assay_error_epsilon: Fraction = Fraction(0)) -> Fraction | None:
    """Sound marginal prospectively replayable regret bound, or abstain.

    The initial calibration certificate and rank proof require independent
    verification and whole-group exchangeability. epsilon is a DECLARED
    deterministic worst-case bound on errors across ALL calibration outcomes
    and future selected assays, not estimated from S669 repeatability.
    """
    if cert.cutoff is None:
        return None
    check_grid(cert.budgets, cert.weights, cert.lambda_margin)
    if budget not in cert.budgets:
        raise AssayCertificateError("adaptive assay budget outside certified grid")
    weight = cert.weights[cert.budgets.index(budget)]
    horizon = evidence_for_panel(panel, budget, measured=measured, strict_deployment=True)
    if horizon.assayed == horizon.total:
        return Fraction(0)
    eta = robust_noise_radius(cert, assay_error_epsilon)
    # q_true <= q_observed + eta and margin_true <= margin_measured+epsilon
    value = (cert.lambda_margin * (horizon.observed_margin+assay_error_epsilon)
             + weight*(cert.cutoff+eta))
    return max(Fraction(0), value)


def retrospectively_audit(cert: AnytimeCertificate, group_evidence: tuple[Horizon, ...],
                          *, assay_error_epsilon: Fraction=Fraction(0)) -> dict:
    """Test-event audit; NEVER use a held-out score to change the certificate."""
    if not group_evidence:
        raise AssayCertificateError("empty held-out group")
    eta=robust_noise_radius(cert,assay_error_epsilon)
    covered = True
    for item in group_evidence:
        if item.actual_regret is None:
            raise AssayCertificateError("missing held-out truth")
        if cert.cutoff is None:
            continue
        w=cert.weights[cert.budgets.index(item.budget)]
        upper=max(Fraction(0),cert.lambda_margin*(item.observed_margin+assay_error_epsilon)
                  + w*(cert.cutoff+eta))
        if item.actual_regret>upper:
            covered=False
    return {"simultaneous_coverage": covered, "calibration_abstained":cert.cutoff is None}
