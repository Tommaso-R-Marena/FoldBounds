"""Exact mutation counterfactuals by a single shared conformational ensemble.

This is an *exact finite-lattice method*, not a physical folding free-energy model.
Every geometry fingerprint encodes all contact pairs, bend count and the
mutation-INDEPENDENT compactness event. In particular, it allows exact
pairwise mutation effects without re-enumerating conformations.

For baseline w(s)>0 and ratio R_i(s)=w_i(s)/w_0(s),
  Z_i = Z_0 E_0[R_i],
  odds_i(A)/odds_0(A) = E_0[R_i | A]/E_0[R_i | not A].
For mutations at sites of the same parity on a bipartite cubic lattice,
there is never a direct contact between them and R_ij(s)=R_i(s)R_j(s).
Thus the two-mutation thermodynamic cycle is *entirely* an ensemble
covariance effect, even though no direct pair energy exists.
"""
from __future__ import annotations
from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction
from math import log
from .mutation_switch import geometric_fingerprints
from .mutation_switch_10 import geometry_fingerprints10

F = Fraction

@dataclass(frozen=True)
class CounterfactualCycle:
    length: int
    initial: str
    sites: tuple[int, int]              # one indexed
    nominal_q: F
    nominal_b: F
    baseline_event_probability: F
    single_i_event_probability: F
    single_j_event_probability: F
    double_event_probability: F
    event_conditional_covariance: F
    complement_conditional_covariance: F
    event_interaction_factor: F
    complement_interaction_factor: F
    odds_epistasis_factor: F
    direct_contact_possible: bool
    enumerated_microstates: int


def _check_sequence(seq):
    if type(seq) is not str or len(seq) not in (8, 10) or set(seq)-{'H','P'}:
        raise ValueError('Only explicit 8 or 10 residue H/P sequences supported')
    return sum(1<<i for i,s in enumerate(seq) if s=='H')


def _fingerprints(n):
    return geometric_fingerprints() if n == 8 else geometry_fingerprints10()


def reweight_cycle(seq: str, pos_i: int, pos_j: int, q: F, b: F) -> CounterfactualCycle:
    """Compute exact 2-site cycle, conditional covariances, and cycle closure.

    Four sequence partition functions are computed from ONE common reference
    microstate histogram with state-specific likelihood ratios. No erroneous
    equality R_ij=R_i*R_j is assumed for directly contacting mutated sites.
    """
    m=_check_sequence(seq)
    n=len(seq)
    if type(pos_i) is not int or type(pos_j) is not int or not 1<=pos_i<=n or not 1<=pos_j<=n or pos_i>=pos_j:
        raise ValueError('Two distinct 1-indexed mutation sites required')
    if type(q) is not F or type(b) is not F or q<=0 or b<=0:
        raise ValueError('Strictly positive exact rational contact/bend weights required')
    i,j=pos_i-1,pos_j-1
    m_i=m^(1<<i);m_j=m^(1<<j);m_ij=m_i^(1<<j)
    sums={True:[F(0) for _ in range(7)],False:[F(0) for _ in range(7)]}
    any_pair_contact=False
    states=0
    for (edges,bends,is_event),count in _fingerprints(n):
        any_pair_contact |= ((1<<i)|(1<<j)) in edges
        k0=sum((m&e)==e for e in edges)
        ki=sum((m_i&e)==e for e in edges)
        kj=sum((m_j&e)==e for e in edges)
        kij=sum((m_ij&e)==e for e in edges)
        weight=F(count)*q**k0*b**bends
        ri=q**(ki-k0)
        rj=q**(kj-k0)
        rij=q**(kij-k0)
        sums[bool(is_event)][0]+=weight
        sums[bool(is_event)][1]+=weight*ri
        sums[bool(is_event)][2]+=weight*rj
        sums[bool(is_event)][3]+=weight*rij
        sums[bool(is_event)][4]+=weight*ri*rj
        states+=count
    if states!=(13565 if n==8 else 308981):
        raise AssertionError('Lossless fingerprint grouping failed')
    e,c=sums[True],sums[False]
    if any(v<=0 for v in e[:4]+c[:4]):
        raise ValueError('Event and complement must have positive weighted masses')
    def p(idx):return e[idx]/(e[idx]+c[idx])
    # Exact conditioning; covariance is w.r.t. reference ensemble in each basin
    def cov(x):return x[4]/x[0]-(x[1]/x[0])*(x[2]/x[0])
    def fac(x):return x[3]*x[0]/(x[1]*x[2])
    # The exact finite-cycle odds interaction is odds_ij * odds_0/(odds_i * odds_j).
    odds_cycle=(e[3]/c[3])*(e[0]/c[0])/((e[1]/c[1])*(e[2]/c[2]))
    result=CounterfactualCycle(n,seq,tuple(sorted((pos_i,pos_j))),q,b,
                               p(0),p(1),p(2),p(3),cov(e),cov(c),
                               fac(e),fac(c),odds_cycle,any_pair_contact,states)
    if not any_pair_contact:
        # Algebraic theorem applied to each conditional basin:
        # E[R_i R_j] / E[R_i]E[R_j] = 1 + Cov(R_i,R_j)/(E[R_i]E[R_j])
        assert e[3]==e[4] and c[3]==c[4]
        assert result.odds_epistasis_factor==result.event_interaction_factor/result.complement_interaction_factor
    return result


def critical_energy_error_kcal(p_original_upper: F, p_mutant_lower: F,
                               kelvin: float=293.15) -> float:
    """Maximum absolute PER-CONFORMATION energy error for guaranteed switch.

    Assuming corrections exp(-delta_i/RT) with each |delta_i|<=epsilon
    and independent adversarial delta_i across all conformations. This is
    NOT a force-field accuracy estimate. Positive margin only if original
    upper<1/2<mutant lower. Common offsets cancel.
    """
    if not (F(0)<p_original_upper<F(1,2)<p_mutant_lower<F(1)):
        raise ValueError('No strict nominal minority-to-majority switch')
    if not (273.0<=kelvin<=400.0):raise ValueError('temperature in Kelvin must be plausible')
    RT=0.00198720425864083*kelvin
    return RT*0.5*min(log(float((1-p_original_upper)/p_original_upper)),
                      log(float(p_mutant_lower/(1-p_mutant_lower))))


def population_with_factor_error(p: F, lo: F, hi: F) -> tuple[F,F]:
    if not all(type(v) is F for v in (p,lo,hi)) or not F(0)<=p<=F(1) or not F(0)<lo<=hi:
        raise ValueError('exact nonnegative probability and positive bounds required')
    return (lo*p/(lo*p+hi*(1-p)), hi*p/(hi*p+lo*(1-p)))