"""Exact, fail-closed certificates for single-valued, atom-resolved NMR restraints.

Scope: source-pinned, converted NEF/PDB *model-vs-restraint* audit.  Neither
NEF NOE intensities nor ambiguity groups are independent exact distances.
Only one-row restraint IDs with explicitly mapped individual atom names are
in the supported dataset; a separate extractor enforces that inclusion rule.

Every status is computed by arbitrary-precision integer comparisons.  If each
endpoint can move independently by a Euclidean vector of norm <= eps, its
separation belongs to [max(0,d-2 eps), d+2 eps]. Both endpoints can attain
these extrema, so classifications are *sharp for each individual pair*.
Joint satisfaction of many restraints is not claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Mapping, Sequence


@dataclass(frozen=True)
class AtomRef:
    residue: int
    atom: str


@dataclass(frozen=True)
class NOERestraint:
    identity: int
    left: AtomRef
    left_position_mA: tuple[int, int, int]
    right: AtomRef
    right_position_mA: tuple[int, int, int]
    lower_mA: int
    upper_mA: int


@dataclass(frozen=True)
class NOEAudit:
    epsilon_mA: int
    evaluated: int
    certainly_satisfied: int
    certainly_violated: int
    undecided: int
    nominal_violated: int
    violation_ids: tuple[int, ...]
    disjoint_witness_ids: tuple[int, ...]

    @property
    def minimum_atoms_requiring_more_than_epsilon_movement(self) -> int:
        """The matching/vertex-cover dual lower bound, NOT a repair proposal."""
        return len(self.disjoint_witness_ids)


def _int(x: object) -> bool:
    return type(x) is int


def _triple(x: object) -> bool:
    return isinstance(x, tuple) and len(x) == 3 and all(_int(v) for v in x)


def _sq(a: tuple[int, int, int], b: tuple[int, int, int]) -> int:
    return sum((x-y)**2 for x, y in zip(a, b))


def verify_record(r: NOERestraint) -> bool:
    if type(r) is not NOERestraint or not _int(r.identity) or r.identity < 1:
        return False
    if type(r.left) is not AtomRef or type(r.right) is not AtomRef:
        return False
    for a in (r.left, r.right):
        if not _int(a.residue) or a.residue < 1 or not isinstance(a.atom, str) or not a.atom or not a.atom.isalnum() or not a.atom.isupper():
            return False
    return (_triple(r.left_position_mA) and _triple(r.right_position_mA)
            and _int(r.lower_mA) and _int(r.upper_mA)
            and 0 <= r.lower_mA <= r.upper_mA and r.left != r.right)


def record_from_row(row: Sequence[object]) -> NOERestraint:
    """Pinned schema: id, residue1, atom1, xyz1, residue2, atom2, xyz2, lo, hi."""
    if not isinstance(row, (list, tuple)) or len(row) != 13:
        raise ValueError('malformed restraint row')
    r = NOERestraint(row[0], AtomRef(row[1], row[2]), tuple(row[3:6]),
                     AtomRef(row[6], row[7]), tuple(row[8:11]), row[11], row[12])
    if not verify_record(r):
        raise ValueError('unsafe or unsupported restraint record')
    return r


def load_pinned_reference(path: str | Path) -> Mapping[str, tuple[NOERestraint, ...]]:
    """Read a *local* source-pinned extraction; do not download or trust labels."""
    j = json.loads(Path(path).read_text(encoding='utf-8'))
    if (type(j) is not dict or j.get('schema') != 'FoldBounds.NOESingleMapped.v1'
        or j.get('unit') != 'milliangstrom' or type(j.get('samples')) is not dict):
        raise ValueError('wrong reference schema')
    out = {}
    for sample, obj in j['samples'].items():
        if type(sample) is not str or type(obj) is not dict:
            raise ValueError('invalid reference grouping')
        records = tuple(record_from_row(r) for r in obj['records'])
        if not records or any(records[i].identity >= records[i+1].identity
                              for i in range(len(records)-1)):
            raise ValueError('duplicate or unordered restraint IDs')
        out[sample] = records
    return out


def classify(r: NOERestraint, epsilon_mA: int = 0) -> str:
    """Exactly one of satisfied, violated, uncertain under *all* bounded moves.

    In milliangstroms, let s = ||p-q||^2 and t=2*epsilon.
    All allowed d' are within [max(0,sqrt(s)-t), sqrt(s)+t].
    Always satisfied iff sqrt(s)-t>=lo and sqrt(s)+t<=hi.
    Always violated iff sqrt(s)-t>hi or sqrt(s)+t<lo.
    Comparisons are reduced to exact integer squares.
    """
    if not verify_record(r) or not _int(epsilon_mA) or epsilon_mA < 0:
        raise ValueError('invalid record or uncertainty radius')
    s = _sq(r.left_position_mA, r.right_position_mA)
    lo, hi, t = r.lower_mA, r.upper_mA, 2 * epsilon_mA
    if hi >= t and s >= (lo + t)**2 and s <= (hi - t)**2:
        return 'satisfied'
    if s > (hi + t)**2 or (lo > t and s < (lo - t)**2):
        return 'violated'
    return 'uncertain'


def check_matching_witness(records: Sequence[NOERestraint],
                           epsilon_mA: int, ids: Sequence[int]) -> bool:
    """Recheck pairwise disjoint, persistently violated edges independently."""
    if (not _int(epsilon_mA) or epsilon_mA < 0 or
        not isinstance(ids, (tuple, list)) or
        any(not _int(k) for k in ids) or len(ids) != len(set(ids))):
        return False
    lookup = {r.identity: r for r in records if verify_record(r)}
    if len(lookup) != len(records):
        return False
    used = set()
    for key in ids:
        r = lookup.get(key)
        if r is None or classify(r, epsilon_mA) != 'violated':
            return False
        if r.left in used or r.right in used:
            return False
        used.add(r.left)
        used.add(r.right)
    return True


def audit(records: Sequence[NOERestraint], epsilon_mA: int = 0) -> NOEAudit:
    """Conservative matching lower bound; no claim of maximum matching."""
    if not _int(epsilon_mA) or epsilon_mA < 0:
        raise ValueError('epsilon must be a nonnegative integer')
    state = {'satisfied': 0, 'violated': 0, 'uncertain': 0}
    nominal = 0
    violated = []
    seen = set()
    records = tuple(records)
    for r in records:
        if not verify_record(r) or r.identity in seen:
            raise ValueError('invalid or duplicate restraint')
        seen.add(r.identity)
        state[classify(r, epsilon_mA)] += 1
        nominal += int(classify(r, 0) == 'violated')
        if classify(r, epsilon_mA) == 'violated':
            violated.append(r)
    used_atoms, witnesses = set(), []
    for r in sorted(violated, key=lambda item: item.identity):
        if r.left not in used_atoms and r.right not in used_atoms:
            used_atoms.add(r.left)
            used_atoms.add(r.right)
            witnesses.append(r.identity)
    if not check_matching_witness(records, epsilon_mA, tuple(witnesses)):
        raise AssertionError('matching certificate failed independent check')
    if sum(state.values()) != len(records):
        raise AssertionError('restraint accounting inconsistency')
    return NOEAudit(epsilon_mA, len(records), state['satisfied'],
                    state['violated'], state['uncertain'], nominal,
                    tuple(r.identity for r in violated), tuple(witnesses))
