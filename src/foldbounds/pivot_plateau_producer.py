"""UNTRUSTED proposal of nonenumerative pivot+corner-move energy plateaus.

A full 2D polymer conformational ensemble is not enumerated. Instead a closed
superlevel component is proposed by exploring its own members. Every proposal
is separately replayed by pivot_plateau_verifier.py.
"""
from __future__ import annotations

from dataclasses import dataclass

Point = tuple[int, int]
Path = tuple[Point, ...]
# All nonidentity elements of the dihedral group D4 of the square.
DIHEDRAL = (
    (0, -1, 1, 0), (-1, 0, 0, -1), (0, 1, -1, 0),
    (1, 0, 0, -1), (-1, 0, 0, 1), (0, 1, 1, 0), (0, -1, -1, 0),
)
COMPASS = ((1, 0), (-1, 0), (0, 1), (0, -1))


class PlateauProposalFailure(ValueError):
    """Refuse to issue an incomplete, higher-energy, or oversized proposal."""


@dataclass(frozen=True)
class PivotPlateauCertificate:
    sequence: str
    seed: Path
    plateau: tuple[Path, ...]
    better_witness: Path


def rectangular_hamiltonian_snake(width: int, height: int) -> Path:
    """Canonical-root Hamiltonian path of a fully occupied width*height box.

    The end or penultimate bead lies in the box's interior for width>=4,
    height>=3, so it is immobile under the local corner/end relocation rule.
    """
    if type(width) is not int or type(height) is not int or width < 4 or height < 3:
        raise ValueError("rectangle requires integer width>=4, height>=3")
    path = [(x, 0) for x in range(width)]
    path.extend((width - 1, y) for y in range(1, height))
    path.extend((x, height-1) for x in range(width-2, -1, -1))
    for y in range(height-2, 0, -1):
        row = list(range(width-1))
        if (height-2-y) % 2:
            row.reverse()
        path.extend((x,y) for x in row)
    assert len(path) == width*height
    return tuple(path)


def hh_score(path: Path, sequence: str) -> int:
    positions = {pos:i for i,pos in enumerate(path)}
    total = 0
    for i, (x,y) in enumerate(path):
        if sequence[i] != 'H':
            continue
        for dx,dy in COMPASS:
            j = positions.get((x+dx,y+dy), -1)
            if j > i+1 and sequence[j] == 'H':
                total += 1
    return total


def proposed_distinct_moves(path: Path):
    """Yield (proposal slot, state), one per valid *slot*, not unique state.

    M = 11*(n-2) uniform slots: local i=2..n-1 with four compass
    choices, and suffix pivots j=1..n-2 with seven D4 nonidentity maps.
    Rejected, identity, and self moves remain self-loops and are omitted.
    """
    n = len(path)
    occupied = set(path)
    for i in range(2, n):
        x, y = path[i-1]
        for a, (dx, dy) in enumerate(COMPASS):
            point = (x+dx, y+dy)
            if point in occupied:
                continue
            if i+1 < n:
                nx, ny = path[i+1]
                if abs(point[0]-nx)+abs(point[1]-ny) != 1:
                    continue
            yield ('local', i, a), path[:i] + (point,) + path[i+1:]
    for pivot in range(1, n-1):
        x0, y0 = path[pivot]
        original = path[pivot+1:]
        protected = set(path[:pivot+1])
        for t, (a,b,c,d) in enumerate(DIHEDRAL):
            moved = tuple((x0+a*(x-x0)+b*(y-y0), y0+c*(x-x0)+d*(y-y0))
                          for x,y in original)
            # The isometry preserves distinctness within the suffix.
            if moved == original or any(p in protected for p in moved):
                continue
            yield ('pivot', pivot, t), path[:pivot+1] + moved


def propose_rectangular_plateau(
    height: int, *, width: int = 4, max_plateau_states: int = 128
) -> PivotPlateauCertificate:
    """Explore only the equal-score closure of a nonoptimal Hamiltonian box.

    For the constructive high-score witness, use a width-8, height/2 box;
    width=4 and even height>=10 guarantee n matches and its score is larger.
    The producer is *not* the assurance boundary; see the independent checker.
    """
    if type(height) is not int or type(width) is not int or width != 4 or height < 10 or height % 2:
        raise ValueError("family requires width=4 and even height>=10")
    if type(max_plateau_states) is not int or max_plateau_states < 1:
        raise ValueError("invalid plateau-state budget")
    seed = rectangular_hamiltonian_snake(width, height)
    sequence = 'H' * len(seed)
    witness = rectangular_hamiltonian_snake(8, height//2)
    k = hh_score(seed, sequence)
    if hh_score(witness, sequence) <= k:
        raise PlateauProposalFailure("constructive witness fails improvement test")
    reached = {seed}
    queue = [seed]
    for state in queue:
        for _, target in proposed_distinct_moves(state):
            score = hh_score(target, sequence)
            if score > k:
                raise PlateauProposalFailure("found an unprotected uphill proposal")
            if score == k and target not in reached:
                reached.add(target)
                queue.append(target)
                if len(reached) > max_plateau_states:
                    raise PlateauProposalFailure("plateau budget exceeded; no partial proof")
    return PivotPlateauCertificate(sequence, seed, tuple(sorted(reached)), witness)
