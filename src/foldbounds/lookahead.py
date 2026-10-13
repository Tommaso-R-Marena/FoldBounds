"""Hierarchical exact-rational bounds for the anchored 2D lattice HP model.

Research prototype: all guarantees concern this finite toy model. This Python
code and its implementation are NOT yet formally verified by Lean.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache

from .matching import matching_score_cap
from .lattice import (
    Certificate, CertificateError, DIRECTIONS, Node, Path, ROOT,
    VerifiedBounds, _check_q, _check_sequence, _contacts, _neighbors,
    _possible_hh_pairs, _validate_path, exact_spectrum,
)


@dataclass(frozen=True)
class Mass:
    e_lo: Fraction = Fraction(0)
    e_hi: Fraction = Fraction(0)
    n_lo: Fraction = Fraction(0)
    n_hi: Fraction = Fraction(0)
    z_hi: Fraction = Fraction(0)


def add(items: list[Mass]) -> Mass:
    return Mass(*(sum((getattr(x, field) for x in items), Fraction(0))
                  for field in ('e_lo', 'e_hi', 'n_lo', 'n_hi', 'z_hi')))


def score_cap(sequence: str, path: Path, pairs: tuple[tuple[int, int], ...]) -> tuple[int, int]:
    """Necessary contact constraints: parity, reachability and H coordination."""
    m, n = len(path), len(sequence)
    current = _contacts(path, pairs)
    occupied = set(path)
    tip = path[-1]
    candidates = 0
    for i, j in pairs:
        if j < m:
            continue  # fixed pairs already included in current
        if i < m:
            x, y = path[i]
            if all((x + dx, y + dy) in occupied for dx, dy in DIRECTIONS):
                continue
            # Chain step j must be within (j-m+1) moves of current tip;
            # to contact i it must lie in the distance-one neighborhood.
            if abs(x-tip[0]) + abs(y-tip[1]) > j-m+2:
                continue
        candidates += 1
    capacity = 0
    for i, residue in enumerate(sequence):
        if residue != 'H':
            continue
        budget = 3 if i == 0 or i == n-1 else 2
        if i < m:
            seen = sum(1 for a,b in pairs if b < m and (a == i or b == i)
                       and abs(path[a][0]-path[b][0])+abs(path[a][1]-path[b][1]) == 1)
            x,y = path[i]
            free = sum((x+dx,y+dy) not in occupied for dx,dy in DIRECTIONS)
            budget = min(max(0, budget-seen), free)
        capacity += budget
    return current, current + min(candidates, capacity // 2)



@lru_cache(maxsize=64)
def walk_count_upper(remaining: int) -> int:
    """Bound count of valid completions by universal square-lattice SAWs.

    Prefix obstacles only remove extensions. The first step has at most
    three orientations (the previous chain site is occupied). For any one
    orientation, the first k steps form one of exactly c_k anchored SAWs.
    Any remaining steps have at most three choices. Compute c_k directly,
    cap k at 8 to keep this cheap; cache the exact integer result.
    """
    if type(remaining) is not int or remaining < 0:
        raise ValueError('remaining must be a nonnegative integer')
    if remaining == 0:
        return 1
    k = min(remaining, 8)
    anchored = sum(exact_spectrum('P'*(k+1)).values())
    return 3 * anchored * 3 ** (remaining-k)

def ray_witnesses(path: Path, length: int, pairs: tuple[tuple[int, int], ...],
                  q: Fraction, threshold: int) -> tuple[Fraction, Fraction]:
    """One explicitly constructed full-length conformation per event class."""
    best_e = best_n = Fraction(0)
    x,y = path[-1]
    for dx,dy in DIRECTIONS:
        future = tuple((x+k*dx,y+k*dy) for k in range(1,length-len(path)+1))
        full = path + future
        if len(set(full)) != length:
            continue
        weight = q ** _contacts(full, pairs)
        if _contacts(full, pairs) >= threshold:
            best_e = max(best_e, weight)
        else:
            best_n = max(best_n, weight)
    return best_e, best_n


def subtree(sequence: str, path: Path, pairs: tuple[tuple[int,int], ...],
            q: Fraction, threshold: int, horizon: int, use_matching: bool = False) -> Mass:
    n = len(sequence)
    if len(path) == n:
        score = _contacts(path, pairs)
        w = q ** score
        return Mass(w,w,0,0,w) if score >= threshold else Mass(0,0,w,w,w)
    score_lo, score_hi = (matching_score_cap(sequence, path, pairs) if use_matching
                          else score_cap(sequence, path, pairs))
    cap = walk_count_upper(n-len(path)) * q ** score_hi
    e_up = cap if score_hi >= threshold else Fraction(0)
    n_up = cap if score_lo < threshold else Fraction(0)
    e_wit,n_wit = ray_witnesses(path,n,pairs,q,threshold)
    base = Mass(e_wit,e_up,n_wit,n_up,cap)
    if horizon == 0:
        return base
    children = add([subtree(sequence,child,pairs,q,threshold,horizon-1,use_matching)
                    for child in _neighbors(path)])
    # Intersect two sound enclosures, guaranteeing a nested bound hierarchy.
    return Mass(max(base.e_lo,children.e_lo),min(base.e_hi,children.e_hi),
                max(base.n_lo,children.n_lo),min(base.n_hi,children.n_hi),
                min(base.z_hi,children.z_hi,base.e_hi+base.n_hi))


def verify_lookahead(cert: Certificate, *, lookahead: int = 0,
                     use_matching: bool = False) -> VerifiedBounds:
    """Independently check a complete frontier; return rational bounds.

    The frontier supplies paths, never trusted mass/energy numbers. Each
    bounded subtree uses a sound contact upper cap, explicit ray witnesses,
    and optionally lookahead levels of exact branching. Runtime can be
    exponential in lookahead; this is not a polynomial-time guarantee.
    """
    _check_sequence(cert.sequence)
    q = _check_q(cert.q)
    if type(cert.event_at_least) is not int:
        raise CertificateError('event threshold must be integer')
    if type(use_matching) is not bool:
        raise ValueError('use_matching must be a boolean')
    if type(lookahead) is not int or lookahead < 0:
        raise ValueError('lookahead must be nonnegative integer')
    n = len(cert.sequence)
    pairs = _possible_hh_pairs(cert.sequence)
    nodes: dict[Path,str] = {}
    prefixes: set[Path] = set()
    for node in cert.frontier:
        if not isinstance(node,Node) or node.kind not in ('bound','exact','dead'):
            raise CertificateError('invalid node')
        _validate_path(node.path,n)
        if node.kind == 'exact' and len(node.path) != n:
            raise CertificateError('exact node must be terminal')
        if node.kind == 'dead' and (len(node.path) == n or _neighbors(node.path)):
            raise CertificateError('dead node has a valid continuation')
        if node.path in nodes:
            raise CertificateError('duplicate node')
        nodes[node.path] = node.kind
        for k in range(2,len(node.path)+1):
            prefixes.add(node.path[:k])
    if ROOT not in prefixes:
        raise CertificateError('missing canonical root')
    seen: set[Path] = set()
    masses: list[Mass] = []
    exact_count = bound_count = 0

    def visit(path: Path) -> None:
        nonlocal exact_count,bound_count
        if path in nodes:
            seen.add(path)
            if nodes[path] == 'dead':
                return  # Zero completions: checked by the validator above.
            masses.append(subtree(cert.sequence,path,pairs,q,cert.event_at_least,lookahead,use_matching))
            if len(path) == n:
                exact_count += 1
            else:
                bound_count += 1
            return
        if len(path) == n:
            raise CertificateError('missing complete conformation')
        for child in _neighbors(path):
            if child not in prefixes:
                raise CertificateError('missing nonempty continuation subtree')
            visit(child)
    visit(ROOT)
    if seen != set(nodes):
        raise CertificateError('shadowed nodes')
    total = add(masses)
    elo,ehi,nlo,nhi = total.e_lo,total.e_hi,total.n_lo,total.n_hi
    if not (0 <= elo <= ehi and 0 <= nlo <= nhi):
        raise AssertionError('invalid mass bounds')
    # The globally straight anchored walk proves that Z > 0.
    p_lo = Fraction(1) if nhi == 0 else (elo/(elo+nhi) if elo else Fraction(0))
    p_hi = Fraction(0) if ehi == 0 else (ehi/(ehi+nlo) if nlo else Fraction(1))
    if not 0 <= p_lo <= p_hi <= 1:
        raise AssertionError('invalid probability bounds')
    return VerifiedBounds(elo+nlo,min(ehi+nhi,total.z_hi),elo,ehi,nlo,nhi,p_lo,p_hi,
                          exact_count,bound_count)
