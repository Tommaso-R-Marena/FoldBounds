"""Constructive lower mass from disjoint *families* of lattice completions.

An untrusted finite 2D HP model-relative reference implementation. Four
monotone quadrant classes are counted by exact integer dynamic programming
without enumerating their individual conformations. Their pairwise overlaps
are precisely the four cardinal rays (for positive suffix length), which are
subtracted once. All arithmetic is exact; no physics or Lean-verified
executable-equivalence claim is made.
"""
from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from fractions import Fraction

from .lattice import Path, _check_q, _check_sequence, _contacts, _possible_hh_pairs, _validate_path


@dataclass(frozen=True)
class MonotoneLower:
    event: Fraction
    complement: Fraction
    distinct_paths: int

    @property
    def partition(self) -> Fraction:
        return self.event+self.complement


def monotone_family_mass(sequence: str, prefix: Path, q: Fraction,
                         event_at_least: int) -> MonotoneLower:
    """Count weighted union of 4 monotone-quadrant suffix families exactly.

    A quadrant permits only two fixed signs of coordinate increments:
    (E,N), (E,S), (W,N), or (W,S). Each step strictly increases a signed
    coordinate potential, so paths within a quadrant never self-intersect
    after the prefix. Fixed prefix sites are explicitly forbidden. No
    nonconsecutive suffix-to-suffix contact can occur: vertices separated
    by >=2 steps differ by >=2 in the monotone potential, and adjacent
    grid sites differ by exactly 1.

    Each straight cardinal ray belongs to exactly two quadrants; other
    positive-length quadrant pairs have empty intersection, and triple
    intersections are empty. Thus sum(quadrants) - sum(rays) is the exact
    weighted mass of a disjoint union of *actual* valid completions.

    Runtime is polynomial in suffix length and number of prefix contacts
    for the restricted families, not in the total number of conformations.
    """
    _check_sequence(sequence)
    _check_q(q)
    if type(event_at_least) is not int:
        raise ValueError("event_at_least must be an integer")
    n,m=len(sequence),len(prefix)
    _validate_path(prefix,n)
    pairs=_possible_hh_pairs(sequence)
    fixed=_contacts(prefix,pairs)
    if m==n:
        weight=q**fixed
        return MonotoneLower(weight if fixed>=event_at_least else Fraction(0),
                             weight if fixed<event_at_least else Fraction(0),1)

    occupied=set(prefix)
    tip=prefix[-1]
    cross_by_index={
        j: tuple(prefix[i] for i,k in pairs if k==j and i<m)
        for j in range(m,n)
    }

    def histogram(moves: tuple[tuple[int,int], ...]) -> dict[int,int]:
        # DP state is (current site, number of prefix-suffix H-H contacts).
        # Nonconsecutive suffix-suffix contacts are provably absent.
        states={(tip,0):1}
        for j in range(m,n):
            nxt=defaultdict(int)
            fixed_h_positions=cross_by_index[j]
            for ((x,y),score),count in states.items():
                for dx,dy in moves:
                    new=(x+dx,y+dy)
                    if new in occupied:
                        continue
                    delta=sum(abs(new[0]-px)+abs(new[1]-py)==1
                              for px,py in fixed_h_positions)
                    nxt[(new,score+delta)]+=count
            states=nxt
            if not states:
                break
        totals=defaultdict(int)
        for (_,score),multiplicity in states.items():
            totals[score]+=multiplicity
        return totals

    east=(1,0); west=(-1,0); north=(0,1); south=(0,-1)
    quadrants=((east,north),(east,south),(west,north),(west,south))
    rays=(east,west,north,south)
    union=defaultdict(int)
    for moves in quadrants:
        for score,count in histogram(moves).items():
            union[score]+=count
    for direction in rays:
        for score,count in histogram((direction,)).items():
            union[score]-=count
    if any(count<0 for count in union.values()):
        raise AssertionError("monotone quadrant inclusion-exclusion failed")
    event=complement=Fraction(0)
    for cross,multiplicity in union.items():
        weight=multiplicity*q**(fixed+cross)
        if fixed+cross>=event_at_least:
            event+=weight
        else:
            complement+=weight
    return MonotoneLower(event,complement,sum(union.values()))
