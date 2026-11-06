"""Proof-carrying near-global-optimal coordinate repair for selected NMR NOEs.

Integer source coordinates are in milliångströms. Integer candidate movement
vectors have resolution 1/1000 milliångström (one millionth Å). The producer
is not trusted: the independent verifier checks all observed distance
intervals, every source atom, a radius cap on all changed atoms, and an
incompatible original restraint giving a universal minimax lower bound.

Scope: the strictly selected, atom-specific NEF subset ONLY. Chemical bonds,
stereochemical feasibility, excluded NOE ambiguities and 3D fold fidelity
are not part of this mathematical optimization. This is not a protein
refinement program or a Lean kernel theorem.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isqrt
from typing import Sequence

from .noe_restraint_audit import AtomRef, NOERestraint, verify_record

SCALE = 1000  # fine ticks per milliångström; 1 tick == 1e-6 Å


@dataclass(frozen=True)
class RepairCertificate:
    sample: str
    source_count: int
    lower_restraint_id: int
    strict_lower_radius_ticks: int
    upper_radius_ticks: int
    moves: tuple[tuple[AtomRef, tuple[int, int, int]], ...]


@dataclass(frozen=True)
class VerifiedRepair:
    sample: str
    restraints: int
    atoms: int
    moved_atoms: int
    lower_witness_id: int
    # Mathematical interval is (strict_lower, upper] in ticks.
    strict_lower_radius_ticks: int
    upper_radius_ticks: int
    actual_max_movement_squared_ticks: int

    @property
    def optimality_gap_angstrom(self) -> float:
        return (self.upper_radius_ticks - self.strict_lower_radius_ticks) / 1_000_000


def _triple(x: object) -> bool:
    return type(x) is tuple and len(x) == 3 and all(type(v) is int for v in x)


def _sq(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    return sum((x - y) ** 2 for x, y in zip(a, b))


def verify_repair(
        records: Sequence[NOERestraint], cert: RepairCertificate
        ) -> VerifiedRepair | None:
    """Independently replay a global minimax bracket, fail closed.

    Let rho* be inf_{Y} max_atom ||Y_atom-X_atom|| subject to ALL listed
    source restraints. A violated pair with original distance d and
    allowed [L,U] implies rho* >= (d-U)/2 or (L-d)/2. We verify a
    STRICT lower cutoff R by exact comparison of integer squared distances.
    A complete integer candidate Y with every interval checked implies
    rho* <= claimed upper radius. No trust in the proposal optimizer.
    """
    try:
        if (type(cert) is not RepairCertificate or type(cert.sample) is not str
            or not cert.sample or type(cert.source_count) is not int
            or type(cert.lower_restraint_id) is not int
            or type(cert.strict_lower_radius_ticks) is not int
            or type(cert.upper_radius_ticks) is not int
            or not 0 <= cert.strict_lower_radius_ticks < cert.upper_radius_ticks
            or type(cert.moves) is not tuple):
            return None
        if len(records) != cert.source_count or not records:
            return None
        source: dict[AtomRef, tuple[int, int, int]] = {}
        by_id: dict[int, NOERestraint] = {}
        for r in records:
            if not verify_record(r) or r.identity in by_id:
                return None
            by_id[r.identity] = r
            for name, xyz in ((r.left, r.left_position_mA),
                              (r.right, r.right_position_mA)):
                if name in source and source[name] != xyz:
                    return None
                source[name] = xyz
        witness = by_id.get(cert.lower_restraint_id)
        if witness is None:
            return None
        s = _sq(witness.left_position_mA, witness.right_position_mA)
        r = cert.strict_lower_radius_ticks
        # No movement of size <=r could correct the selected violation.
        upper_term = SCALE * witness.upper_mA + 2*r
        lower_term = SCALE * witness.lower_mA - 2*r
        exceeds_upper = s*SCALE*SCALE > upper_term*upper_term
        below_lower = lower_term > 0 and s*SCALE*SCALE < lower_term*lower_term
        if not exceeds_upper and not below_lower:
            return None

        move: dict[AtomRef, tuple[int, int, int]] = {}
        for entry in cert.moves:
            if (type(entry) is not tuple or len(entry) != 2
                or type(entry[0]) is not AtomRef or not _triple(entry[1])
                or entry[0] in move or entry[0] not in source
                or entry[1] == (0, 0, 0)):
                return None
            move[entry[0]] = entry[1]
        upper_sq = cert.upper_radius_ticks**2
        if any(_sq(v, (0, 0, 0)) > upper_sq for v in move.values()):
            return None
        def corrected(atom: AtomRef) -> tuple[int, int, int]:
            x = source[atom]
            delta = move.get(atom, (0, 0, 0))
            return tuple(SCALE*q+d for q, d in zip(x, delta))

        for original in records:
            dist_sq = _sq(corrected(original.left), corrected(original.right))
            if not ((SCALE*original.lower_mA)**2 <= dist_sq <=
                    (SCALE*original.upper_mA)**2):
                return None
        max_sq = max((_sq(v,(0,0,0)) for v in move.values()), default=0)
        return VerifiedRepair(cert.sample, len(records), len(source), len(move),
                              witness.identity, r, cert.upper_radius_ticks,
                              max_sq)
    except (ValueError, TypeError, OverflowError, KeyError, AttributeError):
        return None


def certificate_from_json(name: str, raw: dict) -> RepairCertificate:
    """Strict conversion of an external proposal; not a proof by parsing."""
    if type(raw) is not dict or type(name) is not str:
        raise ValueError('invalid proposal')
    try:
        moves=[]
        for atom, dx, dy, dz in raw['moves']:
            if type(atom) is not str or atom.count(':') != 1:
                raise ValueError('invalid atom')
            residue, label=atom.split(':',1)
            if not residue.isdecimal() or not label.isalnum() or not label.isupper():
                raise ValueError('invalid atom name')
            moves.append((AtomRef(int(residue),label),(dx,dy,dz)))
        cert=RepairCertificate(name,raw['source_count'],raw['lower_restraint_id'],
                               raw['strict_lower_radius_ticks'],raw['upper_radius_ticks'],
                               tuple(moves))
        if any(type(z) is not int for z in (cert.source_count,cert.lower_restraint_id,
                cert.strict_lower_radius_ticks,cert.upper_radius_ticks)):
            raise ValueError('invalid numeric certificate')
        return cert
    except (KeyError, TypeError, ValueError) as e:
        raise ValueError('malformed proposal') from e
