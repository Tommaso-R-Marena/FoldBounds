"""Proof-carrying exact minimum-atom-count NOE repairs and sparse radii.

For originally violated atom-pair restraints, a matching of size m forces any
feasible coordinate correction to alter >=m atoms. A candidate with exactly
m moved atoms and all restraints satisfied proves the minimum cardinality is
EXACTLY m. Each matched edge then has exactly one changed endpoint, so one
originally violated matched edge gives a radius lower bound using *rho*, not
2*rho. Integer-only independent replay checks both witnesses.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Sequence

from .noe_restraint_audit import NOERestraint,check_matching_witness
from .noe_minimax_repair import RepairCertificate,VerifiedRepair,verify_repair,SCALE


@dataclass(frozen=True)
class SparseRepairCertificate:
    repair: RepairCertificate
    disjoint_violated_ids: tuple[int, ...]
    lower_witness_id: int
    strict_sparse_lower_radius_ticks: int


@dataclass(frozen=True)
class VerifiedSparseRepair:
    checked_repair: VerifiedRepair
    minimum_number_of_atoms_moved: int
    strict_sparse_lower_radius_ticks: int
    sparse_upper_radius_ticks: int


def verify_sparse_repair(records: Sequence[NOERestraint],
                         cert: SparseRepairCertificate) -> VerifiedSparseRepair | None:
    """Independent integer check of exact support optimality + near-minimax radius."""
    try:
        if (type(cert) is not SparseRepairCertificate or
            type(cert.disjoint_violated_ids) is not tuple or
            type(cert.lower_witness_id) is not int or
            type(cert.strict_sparse_lower_radius_ticks) is not int or
            cert.strict_sparse_lower_radius_ticks < 0):
            return None
        base=verify_repair(records,cert.repair)
        if base is None:return None
        matching=cert.disjoint_violated_ids
        if (not matching or cert.lower_witness_id not in matching
            or len(cert.repair.moves)!=len(matching)
            or not check_matching_witness(records,0,matching)
            or cert.strict_sparse_lower_radius_ticks>=cert.repair.upper_radius_ticks):
            return None
        witnesses={r.identity:r for r in records}
        edge=witnesses.get(cert.lower_witness_id)
        if edge is None:return None
        s=sum((a-b)**2 for a,b in zip(edge.left_position_mA,edge.right_position_mA))
        r=cert.strict_sparse_lower_radius_ticks
        # Since m=m_matching selected moves suffice to cover m disjoint bad
        # edges, every feasible solution with <=m moves has precisely ONE
        # moved endpoint on each matched edge. Thus |d'-d|<=rho, not 2rho.
        high=SCALE*edge.upper_mA+r
        low=SCALE*edge.lower_mA-r
        impossibly_high=s*SCALE*SCALE>high*high
        impossibly_low=low>0 and s*SCALE*SCALE<low*low
        if not impossibly_high and not impossibly_low:
            return None
        return VerifiedSparseRepair(base,len(matching),r,cert.repair.upper_radius_ticks)
    except (ValueError,TypeError,KeyError,OverflowError,AttributeError):
        return None


def sparse_from_json(sample: str, source: dict)->SparseRepairCertificate:
    """Decode external candidate; *only* verify_sparse_repair grants authority."""
    from .noe_minimax_repair import certificate_from_json
    if type(source) is not dict:raise ValueError('invalid source')
    try:
        result=SparseRepairCertificate(
            certificate_from_json(sample,source),
            tuple(source['matching_ids']), source['sparse_lower_restraint_id'],
            source['strict_sparse_lower_radius_ticks'])
        if (type(result.lower_witness_id) is not int or
            any(type(k) is not int for k in result.disjoint_violated_ids)):
            raise ValueError('invalid matching IDs')
        return result
    except (KeyError,ValueError,TypeError) as e:
        raise ValueError('malformed sparse repair proposal') from e
