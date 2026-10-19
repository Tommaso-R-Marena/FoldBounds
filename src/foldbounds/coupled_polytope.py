"""Sharp event-probability enclosure from coupled exact-rational mass bounds.

Given box constraints on event E and complement C and a slab constraint on
Z=E+C, minimize/maximize E/(E+C) exactly, without numerical LP libraries.
These are sharp *only relative to the supplied six bounds*; the true HP
ensemble may satisfy additional constraints not represented by the box/slab.
Tested Python reference mathematics, not Lean-kernel verified.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction


@dataclass(frozen=True)
class CoupledMassEnvelope:
    event_lower: Fraction
    event_upper: Fraction
    complement_lower: Fraction
    complement_upper: Fraction
    partition_lower: Fraction
    partition_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction


def tighten_coupled_masses(event_lower: Fraction, event_upper: Fraction,
                           complement_lower: Fraction, complement_upper: Fraction,
                           partition_lower: Fraction, partition_upper: Fraction
                           ) -> CoupledMassEnvelope:
    """Compute exact sharp P-ratios for a bounded positive-mass 2D polytope.

    Feasible region: event/complement nonnegative; event and complement each
    in their input intervals; Z=event+complement in its input interval.
    Assumes certified positive partition_lower. The optimizer's exact
    probability endpoints are attained in this region (closed compact set).

    e_min=max(Le,Lz-Uc) and c_min=max(Lc,Lz-Ue). For fixed E, P increases
    with E and decreases with C; choose the smallest possible E and then
    largest C=min(Uc,Uz-E) to minimize P. Symmetrically maximize P at the
    smallest possible C and largest E=min(Ue,Uz-C). Both selected pairs are
    feasible whenever the original mass constraints are feasible.
    """
    raw=(event_lower,event_upper,complement_lower,complement_upper,
         partition_lower,partition_upper)
    if not all(isinstance(v,Fraction) for v in raw):
        raise ValueError('all inputs must be exact Fraction values')
    le,ue,lc,uc,lz,uz=raw
    if not (Fraction(0)<=le<=ue and Fraction(0)<=lc<=uc and
            Fraction(0)<lz<=uz):
        raise ValueError('nonnegative event/complement and positive partition intervals required')
    emin=max(le,lz-uc)
    cmin=max(lc,lz-ue)
    emax=min(ue,uz-lc)
    cmax=min(uc,uz-le)
    if (emin>emax or cmin>cmax or le+lc>uz or ue+uc<lz):
        raise ValueError('mass intervals have no feasible joint realization')
    c_at_emin=min(uc,uz-emin)
    e_at_cmin=min(ue,uz-cmin)
    # Feasible witness checks eliminate any hidden zero-denominator cases.
    if (not lc<=c_at_emin or emin+c_at_emin<lz or
            not le<=e_at_cmin or cmin+e_at_cmin<lz):
        raise ValueError('infeasible coupled constraint system')
    pmin=emin/(emin+c_at_emin)
    pmax=e_at_cmin/(e_at_cmin+cmin)
    if not Fraction(0)<=pmin<=pmax<=Fraction(1):
        raise AssertionError('ratio optimization inconsistent')
    return CoupledMassEnvelope(emin,emax,cmin,cmax,
        max(lz,emin+cmin),min(uz,emax+cmax),pmin,pmax)