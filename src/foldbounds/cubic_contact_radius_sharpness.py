"""Sharp adversarial energy-error radius for a checked cubic HP inverse design.

Exact model-specific theorem, based on contact-set differences, plus explicit
adversarial contact-energy vector witnessing a GLOBAL tie at the boundary.
No full lattice enumeration beyond the independently checked high-contact cut.
Not a real-protein design result or Lean-refined executable.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .cubic_cut_certificate import contact_score
from .cubic_robust_inverse_design import (
    CheckedGroundMaps, RobustUniqueOrbit, Walk,
    _d2, _hh_mask, cubic_candidate_pairs
)


@dataclass(frozen=True)
class SharpContactRobustness:
    # All pairwise energy errors |eta_e|<2 quarter units retain exactly one
    # native D4 orbit. At radius exactly 2 this explicit vector produces a
    # globally tied, geometrically distinct competitor.
    critical_error_radius: Fraction
    adversarial_pair_errors: tuple[tuple[int,int], ...]  # bit index, η_e
    competitor_contact_mask: int
    competitor_path: Walk
    missing_target_pairs: int
    selected_energy_at_boundary: int
    competitor_energy_at_boundary: int
    high_sector_max_checked: int
    global_non_ground_energy_margin: int


def prove_sharp_contact_error_radius(
    designed: RobustUniqueOrbit,
    checked_maps: CheckedGroundMaps,
) -> SharpContactRobustness:
    """Prove exact critical robustness radius r*=2, with constructive tie.

    For every 0<=r<2 and EVERY contact perturbation vector |η_e|<=r,
    the original selected native D4 orbit uniquely maximizes energy:
      - same target map but failed distance restraint: exact advantage 2;
      - different K=M map missing a>=1 target edges and adding a others:
        advantage >=(4-2r)*a+2*(1-bonus)>0;
      - K<=M-gap: advantage >=(8-r)*(M-k)>0 for r<=2.

    At r=2 choose an observed same-restraint K=M competitor, set η=-2
    on each target-only contact and η=+2 on each competitor-only contact.
    This yields a GLOBAL tie by the preceding non-strict inequalities.
    Since such a competitor has a different contact map it cannot belong
    to the same D4 orbit. Hence r*=2 is optimal *for this Hamiltonian*.
    """
    if (not isinstance(designed,RobustUniqueOrbit)
            or not isinstance(checked_maps,CheckedGroundMaps)
            or designed.sequence!=checked_maps.sequence
            ):
        raise ValueError('invalid or incompatible checked design')
    if (designed.sequence!=checked_maps.sequence
            or designed.ground_contact_count!=checked_maps.threshold
            or designed.ground_conformations!=checked_maps.high_count
            or designed.gap<1
            or designed.target_orbit_size not in (4,8)):
        raise ValueError('incompatible checked design and ground map census')
    target=designed.designated_contact_mask
    selected=designed.target_path
    if (_hh_mask(designed.sequence,selected)!=target
            or _d2(selected,*designed.target_distance_pair)!=designed.target_squared_distance):
        raise ValueError('selected coordinate witness does not realize target design')
    complete=tuple(checked_maps.entries)
    if len(complete)!=designed.ground_conformations:
        raise ValueError('checked high-conformation census incomplete')
    matches=[]
    pair_order=cubic_candidate_pairs(len(designed.sequence))
    if len(pair_order)>63:
        raise ValueError('contact mask exceeds 63 bits')
    for mask,path in complete:
        if (_hh_mask(designed.sequence,path)!=mask or
                contact_score(designed.sequence,path)!=designed.ground_contact_count):
            raise ValueError('high conformation fails independent contact rescore')
        if (mask!=target and
                _d2(path,*designed.target_distance_pair)==designed.target_squared_distance):
            a=(target&~mask).bit_count()
            b=(mask&~target).bit_count()
            if a!=b or a<1:
                raise ArithmeticError('contact comparison is not same-score replacement')
            matches.append((a,mask,path))
    if not matches:
        raise ValueError('no same-distance competing map: sharp boundary not demonstrated')
    a,alternative,competitor=min(matches)
    perturbations=tuple(sorted(
        [(bit,-2) for bit in range(len(pair_order)) if target&(1<<bit) and not alternative&(1<<bit)]
        + [(bit,+2) for bit in range(len(pair_order)) if alternative&(1<<bit) and not target&(1<<bit)]
    ))
    errors=dict(perturbations)
    def energy(mask,path):
        fav=(mask&target).bit_count()
        rest=(mask&~target).bit_count()
        bonus=int(_d2(path,*designed.target_distance_pair)==designed.target_squared_distance)
        return (8*fav+4*rest+2*bonus
                +sum(value for bit,value in errors.items() if mask&(1<<bit)))
    selected_energy=energy(target,selected)
    challenger_energy=energy(alternative,competitor)
    if selected_energy!=challenger_energy:
        raise ArithmeticError('adversarial radius-2 construction did not yield a tie')
    high_max=max(energy(mask,path) for mask,path in complete)
    if high_max!=selected_energy:
        raise ArithmeticError('a checked high-contact state exceeds proposed tied optimum')
    # The cut certificate excludes intermediate score levels, so for EVERY
    # unseen K<=M-gap geometry the exact contact-set error comparison gives
    # this strict positive energy margin independent of its actual shape.
    low_margin=(8-2)*designed.gap
    if low_margin<=0:
        raise ArithmeticError('no bound separating lower-contact sector')
    return SharpContactRobustness(
        Fraction(2),perturbations,alternative,competitor,a,
        selected_energy,challenger_energy,high_max,low_margin
    )
