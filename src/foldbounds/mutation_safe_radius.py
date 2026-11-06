"""Exact, post-hoc geometry-perturbation radii for binary fold fingerprints.

For K frozen, residue-disjoint reporter pairs and a deposited NMR
conformer with true-fold signed conservative decision margins m_i:

Let q_B be the (B+1)-th SMALLEST signed margin. An arbitrary pairwise
distance perturbation |Delta d_i| < q_B can make at most B otherwise
correct reporters wrong. With up to E erasures and T additional arbitrary
wrong reporter decisions, the strict majority still identifies the
true fold if 2*(B+T)+E<K.

The maximizing integer B is floor((K-E-2T-1)/2).
The worst-case q_B over all explicitly supplied conformers is a
rigorous (but POST-HOC, finite-conformer) geometric certificate.
It is NOT a population-level, physical FRET, unseen-variant or
sequence-based folding guarantee.

If each endpoint CA coordinate is displaced by norm <=r, pairwise
distance shifts by <=2r by the reverse triangle inequality. Therefore
the per-site radius r<q_B/2 is sufficient under the same assumptions.

Do not silently discard observed negative/zero margins; in that case
there is no positive certificate at the corresponding fault budget.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from typing import Mapping, Sequence

from .ga_gb_error_codes import _d2, _ceil_sqrt
from .ga_gb_mutation_safe import MutationSafeCode, reporter_midpoint


@dataclass(frozen=True)
class PerturbationCertificate:
    channels: int
    erased_channels: int
    additional_wrong_channels: int
    geometry_risk_channels: int
    strict_pair_distance_radius_milliangstrom: Fraction
    strict_per_CA_site_radius_milliangstrom: Fraction
    conformer_count: int
    source_only_code: bool
    holds_for_all_submitted_conformers: bool


def sorted_signed_margins(safe:MutationSafeCode,model,truth:str)->tuple[Fraction,...]:
    if not isinstance(safe,MutationSafeCode) or truth not in ("GA95","GB95"):
        raise ValueError("invalid certificate or known fold identity")
    result=[]
    for e in safe.reference.reporters:
        d2=_d2(model.ca[e.a],model.ca[e.b])
        lower,upper=isqrt(d2),_ceil_sqrt(d2)
        mid=reporter_midpoint(e)
        result.append(mid-Fraction(upper) if e.shorter==truth
                      else Fraction(lower)-mid)
    if len(result)!=safe.reference.target_channels:
        raise ValueError("invalid or missing reporter candidates")
    return tuple(sorted(result))


def risk_channels_for_budget(k:int,e:int,t:int)->int:
    if any(type(x) is not int for x in (k,e,t)) or k<1 or e<0 or t<0:
        raise ValueError("invalid exact integer error/erasure budget")
    if 2*t+e>=k:
        raise ValueError("no strict majority possible under this budget")
    return (k-e-2*t-1)//2


def certify_perturbation_radius(safe:MutationSafeCode,
         conformers:Sequence[tuple[str,object]], *,
         erased_channels:int=0,additional_wrong_channels:int=0)->PerturbationCertificate:
    """Compute and verify a finite-ensemble optimal sufficient radius.

    This is an order-statistic margin certificate, not a statistical
    confidence interval and not an assumption of independent noise.
    """
    if not conformers:
        raise ValueError("must have at least one conformer")
    k=safe.reference.target_channels
    b=risk_channels_for_budget(k,erased_channels,additional_wrong_channels)
    ordered=[]
    for truth,model in conformers:
        ranked=sorted_signed_margins(safe,model,truth)
        ordered.append(ranked[b])
    q=min(ordered)
    return PerturbationCertificate(
        k,erased_channels,additional_wrong_channels,b,
        q if q>0 else Fraction(0),
        (q/2) if q>0 else Fraction(0),
        len(conformers),True,q>0,
    )


def check_perturbation_certificate(safe:MutationSafeCode,
          conformers:Sequence[tuple[str,object]],
          witness:PerturbationCertificate)->bool:
    """Reconstruct the entire margin order statistic from raw coordinates."""
    if not isinstance(witness,PerturbationCertificate):
        raise ValueError("not a valid margin witness")
    expected=certify_perturbation_radius(
        safe,conformers,erased_channels=witness.erased_channels,
        additional_wrong_channels=witness.additional_wrong_channels)
    if witness!=expected:
        raise ValueError("forged or unsupported geometric perturbation radius")
    return True


def lossless_radius_A(v:Fraction)->str:
    """Render an exact rational mÅ radius in Å without floating rounding."""
    if not isinstance(v,Fraction) or v<0:
        raise ValueError("radius must be an exact nonnegative Fraction")
    r=v/Fraction(1000)
    return f"{r.numerator}/{r.denominator}"
