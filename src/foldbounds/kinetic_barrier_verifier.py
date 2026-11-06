"""Independent complete-state replay for rooted 2D HP kinetic cut certificates.

The move rule is: choose uniformly one of n-2 movable beads (indices 2..n-1)
and one of the four neighboring sites of its preceding bead. Reject invalid
relocations; accept valid moves x->y with min(1, q**(k(y)-k(x))).
This symmetric proposal obeys detailed balance even if the state graph is
not irreducible. All rational probabilities are exact, for q >= 1.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .kinetic_barriers import BasinCertificate

NEIGHBORS = ((0, 1), (0, -1), (1, 0), (-1, 0))


class KineticCertificateError(ValueError):
    """Proposed basin or its required completeness conditions are invalid."""


@dataclass(frozen=True)
class KineticCutResult:
    total_states: int
    native_score: int
    native_states: int
    native_equilibrium_probability: Fraction
    basin_states: int
    basin_top_score: int
    basin_top_states: int
    seed_score: int
    boundary_edges: int
    boundary_top_score: int | None
    barrier_order: int | None
    proposal_slots: int
    basin_mass: Fraction
    boundary_flux: Fraction
    stationary_exit_per_step: Fraction
    all_q_coefficient: Fraction
    all_q_exit_per_step_bound: Fraction
    exit_by_steps_upper: Fraction
    deterministic_seed_exit_by_steps_upper: Fraction
    deterministic_all_q_exit_per_step_bound: Fraction
    optimal_cut_verified_at_q: bool
    closed_class: bool


def _walks_independent(n: int, cap: int) -> list[tuple]:
    # Alternative iterative depth-first enumeration and direction ordering.
    root = ((0, 0), (1, 0))
    pending = [(root, frozenset(root))]
    result = []
    while pending:
        p, occupied = pending.pop()
        if len(p) == n:
            result.append(p)
            if len(result) > cap:
                raise KineticCertificateError("walk budget exceeded; fail closed")
            continue
        x, y = p[-1]
        for dx, dy in NEIGHBORS:
            new = (x + dx, y + dy)
            if new not in occupied:
                pending.append((p + (new,), occupied | {new}))
    return result


def _score_independent(path: tuple, sequence: str) -> int:
    k = 0
    for i in range(len(path)):
        if sequence[i] != "H":
            continue
        for j in range(i + 3, len(path)):
            if sequence[j] == "H":
                if abs(path[i][0] - path[j][0]) + abs(path[i][1] - path[j][1]) == 1:
                    k += 1
    return k


def _moves_independent(path: tuple) -> set[tuple]:
    out = set()
    used = set(path)
    n = len(path)
    for i in range(2, n - 1):
        a, b, old = path[i-1], path[i+1], path[i]
        # On a square lattice, only the opposite corner of a 1x1 square
        # can be a second common neighbor of both adjacent chain beads.
        other = (a[0] + b[0] - old[0], a[1] + b[1] - old[1])
        if other != old and other not in used and (
            abs(a[0]-other[0])+abs(a[1]-other[1]) == 1 and
            abs(b[0]-other[0])+abs(b[1]-other[1]) == 1
        ):
            out.add(path[:i] + (other,) + path[i+1:])
    x, y = path[-2]
    for dx, dy in NEIGHBORS:
        end = (x + dx, y + dy)
        if end != path[-1] and end not in used:
            out.add(path[:-1] + (end,))
    return out


def verify_kinetic_basin(
    certificate: BasinCertificate,
    q: Fraction,
    steps: int,
    *,
    max_states: int = 40000,
) -> KineticCutResult:
    """Recompute the entire finite ensemble, basin closure, and cut weights.

    For X_0 ~ pi(.|S), pi(x) proportional to q**k(x), and any t >= 0:
        Pr(hit S^c by t) <= min(1, t * Q(S,S^c) / pi(S)).
    Since pi(.|S) <= pi(.)/pi(S) coordinatewise, this follows from
    stationarity and a union bound, without an independence assumption.
    No native state lies in S, hence the same bound holds for native hits.

    For all real q>=1, additionally Q(S,S^c)/pi(S) <=
        B/(4(n-2)*a) * q**(-d),
    where B is the number of directed cut edges, a is the count of maximum
    score states in S, d=max_score(S)-max_cut_min_score.
    """
    if not isinstance(certificate, BasinCertificate):
        raise KineticCertificateError("expected a BasinCertificate")
    seq = certificate.sequence
    if not isinstance(seq, str) or len(seq) < 4 or set(seq) - {"H", "P"}:
        raise KineticCertificateError("invalid HP sequence")
    if not isinstance(q, Fraction) or q < 1:
        raise KineticCertificateError("exact rational q >= 1 required")
    if type(steps) is not int or steps < 0:
        raise KineticCertificateError("steps must be a nonnegative integer")
    if type(max_states) is not int or max_states < 1:
        raise KineticCertificateError("invalid budget")
    if certificate.mode not in ("superlevel", "arbitrary", "exact-mincut"):
        raise KineticCertificateError("unknown basin certificate mode")
    if type(certificate.threshold) is not int or certificate.threshold < (1 if certificate.mode == "superlevel" else 0):
        raise KineticCertificateError("invalid integer threshold")
    if certificate.mode == "exact-mincut":
        if certificate.flow_q != q or not isinstance(certificate.flow_q, Fraction):
            raise KineticCertificateError("maxflow witness is only valid at its exact q")
    elif certificate.flow_q is not None or certificate.flow_edges:
        raise KineticCertificateError("unexpected flow metadata for unoptimized cut")

    paths = _walks_independent(len(seq), max_states)
    universe = set(paths)
    members = certificate.members
    if not isinstance(members, tuple) or not members:
        raise KineticCertificateError("empty or malformed basin")
    try:
        S = set(members)
        duplicate = len(S) != len(members)
        valid_seed = certificate.seed in S
    except (TypeError, ValueError) as exc:
        raise KineticCertificateError("malformed path data") from exc
    if duplicate:
        raise KineticCertificateError("duplicate basin members")
    if not S <= universe or not valid_seed:
        raise KineticCertificateError("unknown path or seed outside submitted basin")
    scores = {p: _score_independent(p, seq) for p in paths}
    if certificate.mode == "superlevel" and any(scores[p] < certificate.threshold for p in S):
        raise KineticCertificateError("basin includes a below-threshold state")
    score_max = max(scores.values())
    top = max(scores[p] for p in S)
    if top >= score_max:
        raise KineticCertificateError("basin contains global native ground states")
    # Check connectedness and completeness independently, over the entire cut.
    reached = {certificate.seed}
    queue = [certificate.seed]
    cut: list[tuple] = []
    for p in S:
        for other in _moves_independent(p):
            if other not in universe:
                raise KineticCertificateError("move outside enumerated state universe")
            if p not in _moves_independent(other):
                raise KineticCertificateError("asymmetric move graph")
            if other not in S:
                if certificate.mode == "superlevel" and scores[other] >= certificate.threshold:
                    raise KineticCertificateError("omitted high-score basin neighbor")
                cut.append((p, other))
    for p in queue:
        for other in _moves_independent(p):
            if other in S and other not in reached:
                reached.add(other)
                queue.append(other)
    if reached != S:
        raise KineticCertificateError("disconnected/extra basin members")

    Z = sum((q ** s for s in scores.values()), Fraction(0))
    native_count = sum(s == score_max for s in scores.values())
    native_p = native_count * q**score_max / Z
    ZS = sum((q ** scores[p] for p in S), Fraction(0))
    slots = 4 * (len(seq) - 2)
    flux = sum((q ** min(scores[a], scores[b]) for a, b in cut), Fraction(0))
    leak = flux / (slots * ZS)
    max_cut = max((min(scores[a], scores[b]) for a, b in cut), default=None)
    top_count = sum(scores[p] == top for p in S)
    seed_score = scores[certificate.seed]
    coefficient = Fraction(len(cut), slots * top_count)
    barrier = top - max_cut if max_cut is not None else None
    all_q = coefficient * (q**(-barrier)) if barrier is not None else Fraction(0)
    deterministic_bound = min(Fraction(1), steps * flux / (slots * q**seed_score))
    seed_all_q = (Fraction(len(cut), slots) * q**(max_cut-seed_score)
                  if max_cut is not None else Fraction(0))
    if leak > all_q or leak > 1 or flux/(slots*q**seed_score) > seed_all_q:
        raise AssertionError("mathematical cut inequality violated")

    optimal_verified = False
    if certificate.mode == "exact-mincut":
        # Exact integral feasible flow plus matching cut certifies optimality.
        # The verifier does not trust the producer's max-flow algorithm.
        scale = q.denominator ** score_max
        balance: dict[tuple, int] = {}
        used_edges: set[frozenset] = set()
        if not isinstance(certificate.flow_edges, tuple):
            raise KineticCertificateError("malformed flow witness")
        for item in certificate.flow_edges:
            try:
                u, v, amount = item
                if u not in universe or v not in universe or u == v:
                    raise KineticCertificateError("flow references unknown states")
                if type(amount) is not int or amount == 0:
                    raise KineticCertificateError("flow amount must be nonzero integer")
                key = frozenset((u, v))
            except (TypeError, ValueError) as exc:
                raise KineticCertificateError("malformed flow edge") from exc
            if key in used_edges or v not in _moves_independent(u):
                raise KineticCertificateError("duplicate flow edge or nonlocal transition")
            used_edges.add(key)
            cap = scale * q**min(scores[u],scores[v])
            if cap.denominator != 1 or abs(amount) > cap:
                raise KineticCertificateError("flow exceeds exact edge capacity")
            balance[u] = balance.get(u, 0) + amount
            balance[v] = balance.get(v, 0) - amount
        for p, val in balance.items():
            if p == certificate.seed:
                if val < 0:
                    raise KineticCertificateError("negative source outflow")
            elif scores[p] == score_max:
                if val > 0:
                    raise KineticCertificateError("native sink has net outflow")
            elif val != 0:
                raise KineticCertificateError("flow conservation violated")
        total_flow = balance.get(certificate.seed, 0)
        native_inflow = -sum(val for p,val in balance.items() if scores[p] == score_max)
        cut_capacity = scale * flux
        if cut_capacity.denominator != 1 or total_flow != native_inflow or total_flow != cut_capacity:
            raise KineticCertificateError("maxflow / cut duality witness mismatched")
        optimal_verified = True

    return KineticCutResult(
        len(paths), score_max, native_count, native_p,
        len(S), top, top_count, seed_score, len(cut), max_cut, barrier, slots,
        ZS, flux, leak, coefficient, all_q,
        min(Fraction(1), steps*leak), deterministic_bound, seed_all_q,
        optimal_verified, not cut,
    )
