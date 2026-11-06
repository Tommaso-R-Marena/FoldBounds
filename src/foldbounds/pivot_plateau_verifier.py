"""Independent nonenumerative verifier for local+pivot kinetic metastability.

The verifier never calls the producer's geometry, score or move functions.
It audits EVERY one of 11(n-2) possible slots from EVERY supplied plateau
member. The full conformational space is neither generated nor traversed.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction

from .pivot_plateau_producer import PivotPlateauCertificate


class PlateauCertificateError(ValueError):
    """Malformed or incomplete proposed kinetic plateau; fail closed."""


@dataclass(frozen=True)
class PivotPlateauResult:
    n: int
    plateau_states: int
    plateau_score: int
    witness_score: int
    witness_gap: int
    proposal_slots: int
    boundary_proposal_slots: int
    boundary_penalties: tuple[tuple[int, int], ...]
    max_exit_proposals_per_state: int
    conditional_gibbs_exit_one_step: Fraction
    fixed_seed_exit_one_step_upper: Fraction
    fixed_seed_hit_by_steps_upper: Fraction
    conditional_gibbs_hit_by_steps_upper: Fraction
    all_q_fixed_seed_exit_coefficient: Fraction
    equilibrium_better_event_lower: Fraction
    equilibrium_bad_event_upper: Fraction
    total_enumerated_conformations: int


def _geometry_ok(path: tuple, n: int) -> bool:
    if not isinstance(path, tuple) or len(path) != n:
        return False
    if path[:2] != ((0,0),(1,0)):
        return False
    for point in path:
        if not isinstance(point, tuple) or len(point) != 2 or any(type(v) is not int for v in point):
            return False
    if len(frozenset(path)) != n:
        return False
    return all(abs(a[0]-b[0])+abs(a[1]-b[1]) == 1
               for a,b in zip(path, path[1:]))


def _independent_contacts(path: tuple, sequence: str) -> int:
    # Independent definition: count only positive-coordinate geometric edges,
    # excluding chain bonds. O(n) via a separate occupancy lookup, not the
    # producer's four-neighbor/ordered-index implementation.
    indices = {u:i for i,u in enumerate(path)}
    count = 0
    for i, (x,y) in enumerate(path):
        if sequence[i] != 'H':
            continue
        for v in ((x+1,y),(x,y+1)):
            j = indices.get(v, -1)
            if j >= 0 and sequence[j] == 'H' and abs(i-j) > 1:
                count += 1
    return count


def _state_proposal_targets(path: tuple):
    """Independent slot enumeration: local geometric corners + D4 action.

    Multiplicities are preserved. The local term uses the other diagonal
    corner and independent endpoint checks; the pivot term uses repeated
    quarter turns plus a reflection, not the producer's matrices.
    """
    n = len(path)
    taken = frozenset(path)
    for i in range(2,n-1):
        left, p, right = path[i-1], path[i], path[i+1]
        alt = (left[0]+right[0]-p[0], left[1]+right[1]-p[1])
        if (alt not in taken and
                abs(alt[0]-left[0])+abs(alt[1]-left[1]) == 1 and
                abs(alt[0]-right[0])+abs(alt[1]-right[1]) == 1):
            yield path[:i] + (alt,) + path[i+1:]
    x, y = path[-2]
    for dx,dy in ((0,1),(1,0),(0,-1),(-1,0)):
        new_end = (x+dx, y+dy)
        if new_end not in taken:
            yield path[:-1] + (new_end,)
    for pivot in range(1,n-1):
        x0,y0 = path[pivot]
        protected = set(path[:pivot+1])
        tail = tuple((x-x0,y-y0) for x,y in path[pivot+1:])
        for reflect in (False,True):
            for turns in range(4):
                if not reflect and turns == 0:
                    continue
                transformed = []
                for x,y in tail:
                    if reflect:
                        y = -y
                    for _ in range(turns):
                        x,y = -y,x
                    transformed.append((x0+x,y0+y))
                end = tuple(transformed)
                if end == path[pivot+1:] or any(v in protected for v in end):
                    continue
                yield path[:pivot+1] + end


def verify_pivot_plateau(
    cert: PivotPlateauCertificate, q: Fraction, steps: int, *,
    max_states: int = 64, max_beads: int = 1600,
) -> PivotPlateauResult:
    """Certify strong-ergodicity-compatible kinetic/equilibrium separation.

    For fixed seed and target event A={k>k_plateau}, with at most one native-
    improving witness required:
      P_seed(hit A by T) <= min(1, T max_{u in S} hazard(u)).
    This holds from any seed in S, without equilibration, irreducibility,
    mixing, or a full lattice-ensemble enumeration.

    Each boundary proposal has a loss d>=1, giving hazard(u) =
      sum_{valid exit slots from u} q**(-d) / [11(n-2)].
    Therefore for ALL REAL q>=1, hazard(u) <= B_u/[11(n-2) q].

    With a valid better witness of contact score k+gap and <=3**(n-2)
    canonical rooted SAWs, the global Gibbs measure obeys
      P_eq(k>k0) >= q**gap/(q**gap + 3**(n-2)).

    Exact rational checks: no producer algorithms, no untrusted score/cuts.
    """
    if not isinstance(cert, PivotPlateauCertificate):
        raise PlateauCertificateError("certificate type mismatch")
    if not isinstance(q, Fraction) or q < 1:
        raise PlateauCertificateError("q must be an exact Fraction >=1")
    if type(steps) is not int or steps < 0:
        raise PlateauCertificateError("time horizon must be a nonnegative integer")
    if type(max_states) is not int or max_states < 1 or type(max_beads) is not int or max_beads < 4:
        raise PlateauCertificateError("invalid verifier budget")
    sequence = cert.sequence
    if not isinstance(sequence, str) or len(sequence) < 4 or len(sequence) > max_beads or set(sequence)-{'H','P'}:
        raise PlateauCertificateError("invalid sequence or max_beads exceeded")
    n = len(sequence)
    if not isinstance(cert.plateau, tuple) or len(cert.plateau) == 0 or len(cert.plateau) > max_states:
        raise PlateauCertificateError("empty or over-budget plateau")
    try:
        S = set(cert.plateau)
    except (ValueError,TypeError) as exc:
        raise PlateauCertificateError("malformed plateau paths") from exc
    if len(S) != len(cert.plateau):
        raise PlateauCertificateError("duplicate plateau member")
    if not all(_geometry_ok(p,n) for p in S) or not _geometry_ok(cert.better_witness,n):
        raise PlateauCertificateError("noncanonical, self-intersecting or invalid path")
    if not _geometry_ok(cert.seed,n) or cert.seed not in S:
        raise PlateauCertificateError("invalid starting conformation")
    k0 = _independent_contacts(cert.seed,sequence)
    if any(_independent_contacts(p,sequence) != k0 for p in S):
        raise PlateauCertificateError("plateau contains unequal-energy states")
    kw = _independent_contacts(cert.better_witness,sequence)
    gap = kw-k0
    if gap <= 0:
        raise PlateauCertificateError("better-state witness does not improve contacts")
    if cert.better_witness in S:
        raise PlateauCertificateError("better witness lies in plateau")
    score_cache = {}
    rows = {}
    adjacency = {}
    total_losses = Counter()
    all_boundary = 0
    max_exit_slots = 0
    for p in S:
        losses = Counter()
        eq_neighbors = set()
        for nxt in _state_proposal_targets(p):
            if nxt in score_cache:
                k = score_cache[nxt]
            else:
                k = _independent_contacts(nxt, sequence)
                score_cache[nxt] = k
            if nxt in S:
                if k != k0:
                    raise PlateauCertificateError("unequal-energy internal transition")
                eq_neighbors.add(nxt)
            else:
                if k >= k0:
                    raise PlateauCertificateError("omitted neutral/uphill move: cut not closed")
                losses[k0-k] += 1
        rows[p] = losses
        adjacency[p] = eq_neighbors
        total_losses.update(losses)
        n_out = sum(losses.values())
        all_boundary += n_out
        max_exit_slots = max(max_exit_slots,n_out)
    # The submitted set is exactly the seed's zero-barrier connected component.
    known = {cert.seed}
    todo = [cert.seed]
    for p in todo:
        for nxt in adjacency[p]:
            if nxt not in known:
                known.add(nxt)
                todo.append(nxt)
    if known != S:
        raise PlateauCertificateError("extraneous/disconnected plateau members")
    if not all_boundary:
        # A globally irreducible pivot chain cannot have a nontrivial
        # closed nonnative class; still mathematically safe if no boundary.
        raise PlateauCertificateError("not a metastable cut: no exit proposals")
    M = 11*(n-2)
    row_exit = {p:sum((count*q**(-d) for d,count in losses.items()),Fraction())
                for p, losses in rows.items()}
    hazard = max(row_exit.values())/M
    global_q_ceiling = Fraction(max_exit_slots, M) / q
    if hazard > global_q_ceiling:
        raise AssertionError("all-q barrier certificate inconsistent")
    gibbs_hazard = sum(row_exit.values(),Fraction()) / (M*len(S))
    # One witness lower bounds the target partition mass; at most 3**(n-2)
    # self-avoiding continuations of the first fixed root bond exist.
    upper_count = 3**(n-2)
    q_gap = q**gap
    pi_good_lo = q_gap/(q_gap + upper_count)
    seed_hit_hi = min(Fraction(1), steps*hazard)
    gibbs_hit_hi = min(Fraction(1), steps*gibbs_hazard)
    return PivotPlateauResult(n,len(S),k0,kw,gap,M,all_boundary,
       tuple(sorted(total_losses.items())),max_exit_slots,
       gibbs_hazard,hazard,seed_hit_hi,gibbs_hit_hi,
       Fraction(max_exit_slots,M),pi_good_lo,1-pi_good_lo,0)
