"""Finite verified short-block counts yield long-chain SAW completion caps.

For any fixed legal 2D SAW prefix and r remaining sites, split the r
steps into blocks k1,...,kj. Every full extension injectively encodes
a tuple of free rooted self-avoiding block walks, so its number is at most
A(k1)*...*A(kj), where each A(k) is *exactly enumerated* for k <= 11.
These are software/model-relative combinatorial certificates, not
Lean-kernel verification or real-protein experimental statements.
"""
from __future__ import annotations
from dataclasses import dataclass
from functools import lru_cache

ROOT=((0,0),(1,0))
DIRECTIONS=((1,0),(0,1),(-1,0),(0,-1))


@dataclass(frozen=True)
class SAWCountCertificate:
    remaining: int
    block_sizes: tuple[int, ...]
    count_upper: int
    block_max: int


@lru_cache(maxsize=24)
def exact_rooted_block_counts(block_max: int = 8) -> tuple[int, ...]:
    """Exactly count anchored two-dimensional SAWs for block lengths 0..max."""
    if type(block_max) is not int or not 1 <= block_max <= 11:
        raise ValueError("block_max must be an integer from 1 to 11")
    paths = [ROOT]
    counts = [1]
    for _ in range(block_max):
        next_paths = []
        for path in paths:
            x,y = path[-1]
            occupied = set(path)
            for dx,dy in DIRECTIONS:
                site=(x+dx,y+dy)
                if site not in occupied:
                    next_paths.append(path+(site,))
        paths=next_paths
        counts.append(len(paths))
    return tuple(counts)


def verify_saw_count_cap(cert: SAWCountCertificate) -> bool:
    """Re-enumerate base blocks and independently validate the product."""
    if not isinstance(cert, SAWCountCertificate):
        raise ValueError("incorrect SAW count certificate type")
    if (type(cert.remaining) is not int or cert.remaining < 0 or
            type(cert.block_max) is not int or not 1 <= cert.block_max <= 11 or
            type(cert.count_upper) is not int or cert.count_upper < 1 or
            not isinstance(cert.block_sizes, tuple) or
            any(type(k) is not int or not 1 <= k <= cert.block_max
                for k in cert.block_sizes)):
        raise ValueError("malformed SAW count certificate")
    if sum(cert.block_sizes) != cert.remaining:
        raise ValueError("step-block partition has wrong length")
    counts = exact_rooted_block_counts(cert.block_max)
    product = 1
    for k in cert.block_sizes:
        product *= counts[k]
    if product != cert.count_upper:
        raise ValueError("SAW cap does not equal product of verified block counts")
    return True


@lru_cache(maxsize=256)
def optimal_saw_count_cap(remaining: int, block_max: int = 8) -> SAWCountCertificate:
    """Minimize exact short-block-count products by integer dynamic programming."""
    if type(remaining) is not int or remaining < 0:
        raise ValueError("remaining must be nonnegative integer")
    counts=exact_rooted_block_counts(block_max)
    table: list[tuple[int, tuple[int, ...]]] = [(1,())]
    for r in range(1,remaining+1):
        table.append(min((table[r-k][0]*counts[k],table[r-k][1]+(k,))
                         for k in range(1,min(r,block_max)+1)))
    count, blocks = table[-1]
    proof=SAWCountCertificate(remaining,blocks,count,block_max)
    verify_saw_count_cap(proof)
    assert count <= 3**remaining
    return proof
