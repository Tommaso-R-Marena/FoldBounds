"""Exact disjoint union of monotone-quadrant and nonmonotone <=2-turn suffixes.

Every candidate path is recomputed from cardinal directions, checked for full
self avoidance, and scored by exact HP contacts. The monotone restricted family
is counted by the independent integer DP; extra paths are filtered to exclude
its precisely characterized overlap. This is Python-tested, NOT a Lean proof.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations,product
from .lattice import (DIRECTIONS,Path,_check_q,_check_sequence,
                      _validate_path,_contacts,_possible_hh_pairs)
from .monotone import monotone_family_mass

_QUADS=(frozenset(((1,0),(0,1))),frozenset(((1,0),(0,-1))),
        frozenset(((-1,0),(0,1))),frozenset(((-1,0),(0,-1))))

@dataclass(frozen=True)
class UnionLower:
    event: Fraction
    complement: Fraction
    distinct_paths: int
    nonmonotone_turn_paths: int

def monotone_turn_union(sequence: str,prefix: Path,q: Fraction,
                        event_at_least: int,max_turns: int=2) -> UnionLower:
    """Exactly count monotone union (low-turn \\ monotone), both event classes.

    For at most two turns, every direction word has <=3 constant-direction
    blocks. Unique positive block lengths and unequal consecutive directions
    give exactly one representation. This is polynomial in suffix length for
    the fixed turn budget; prefix-core exploration can still be exponential.
    """
    _check_sequence(sequence)
    _check_q(q)
    _validate_path(prefix,len(sequence))
    if type(event_at_least) is not int or type(max_turns) is not int or not 0<=max_turns<=2:
        raise ValueError('invalid threshold or maximum turns')
    base=monotone_family_mass(sequence,prefix,q,event_at_least)
    n,m=len(sequence),len(prefix)
    r=n-m
    if r<2 or max_turns==0:
        return UnionLower(base.event,base.complement,base.distinct_paths,0)
    pairs=_possible_hh_pairs(sequence)
    event=complement=Fraction(0)
    extra=0
    fixed=set(prefix)
    for changes in range(1,min(max_turns,r-1)+1):
        for cuts in combinations(range(1,r),changes):
            lens=tuple(b-a for a,b in zip((0,)+cuts,cuts+(r,)))
            for blocks in product(DIRECTIONS,repeat=changes+1):
                if any(a==b for a,b in zip(blocks,blocks[1:])):
                    continue
                if any(set(blocks)<=cone for cone in _QUADS):
                    continue  # Already included EXACTLY ONCE in monotone DP.
                x,y=prefix[-1]
                seen=set(fixed)
                full=list(prefix)
                valid=True
                for (dx,dy),length in zip(blocks,lens):
                    for _ in range(length):
                        x+=dx
                        y+=dy
                        if (x,y) in seen:
                            valid=False
                            break
                        seen.add((x,y))
                        full.append((x,y))
                    if not valid:
                        break
                if not valid:
                    continue
                score=_contacts(tuple(full),pairs)
                weight=q**score
                if score>=event_at_least:
                    event+=weight
                else:
                    complement+=weight
                extra+=1
    return UnionLower(base.event+event,base.complement+complement,
                      base.distinct_paths+extra,extra)
