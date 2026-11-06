"""Exact finite-model ground-state population and ambiguity certificates.

Applies only to anchored all-H simple-cubic SAWs of lengths 2^3 or 3^3.
The exact Hamiltonian counts and short-block SAW numbers are independently
replayable via the C++ oracles packaged with this research module. Their
mathematical assumptions are documented; no Lean executable refinement or
experimental protein-folding prediction is claimed.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction

from .filled_box_contact_rigidity import contact_map_confidence_upper

# Canonical first +x Hamiltonian paths for side 2 and 3; checked by two C++
# complete-search implementations (the second audits every contact mask).
CUBE_GROUND_COUNTS = {2: 24, 3: 826768}

# Number of unrestricted rooted cubic self-avoiding suffixes of k steps,
# after anchoring the first +x bond; separately enumerated in C++.
CUBIC_FREE_BLOCK_COUNTS = (
    1, 5, 25, 121, 589, 2821, 13565, 64661, 308981, 1468313, 6989025,
)


@dataclass(frozen=True)
class CubicBlockUpper:
    extra_steps: int
    blocks: tuple[int, ...]
    bound: int


@dataclass(frozen=True)
class CubicGroundPopulation:
    side: int
    q: Fraction
    minimum_contacts: int
    correction_lower: Fraction
    correction_upper: Fraction
    ground_paths: int
    ground_contact_count: int
    total_walk_upper: CubicBlockUpper
    event_probability_lower: Fraction
    single_fold_probability_upper: Fraction
    one_ground_contact_map_probability_upper: Fraction
    ground_contact_map_count: int


def _optimal_block_product(r: int) -> tuple[int, tuple[int, ...]]:
    counts = CUBIC_FREE_BLOCK_COUNTS
    dp = [1]
    ways: list[tuple[int, ...]] = [()]
    for remaining in range(1, r + 1):
        candidate = min(
            (dp[remaining - k] * counts[k], ways[remaining - k] + (k,))
            for k in range(1, min(remaining, len(counts) - 1) + 1)
        )
        dp.append(candidate[0])
        ways.append(candidate[1])
    return dp[r], ways[r]


def optimal_cubic_walk_cap(extra_steps: int) -> CubicBlockUpper:
    """Exact-minimal product of independently enumerated local walk caps.

    A legal extension divided into contiguous blocks injects into the
    Cartesian product of its local blocks after forgetting older obstacles.
    Thus the product is a rigorous *count upper bound*, not an asymptotic fit.
    """
    if type(extra_steps) is not int or extra_steps < 0 or extra_steps > 10000:
        raise ValueError("extra_steps must be an integer in [0,10000]")
    value, blocks = _optimal_block_product(extra_steps)
    return CubicBlockUpper(extra_steps, blocks, value)


def verify_cubic_walk_cap(cert: CubicBlockUpper) -> bool:
    if not isinstance(cert, CubicBlockUpper):
        return False
    if type(cert.extra_steps) is not int or not 0 <= cert.extra_steps <= 10000:
        return False
    if (not isinstance(cert.blocks, tuple) or
        any(type(k) is not int or not 1 <= k < len(CUBIC_FREE_BLOCK_COUNTS)
            for k in cert.blocks)):
        return False
    if sum(cert.blocks) != cert.extra_steps:
        return False
    product = 1
    for k in cert.blocks:
        product *= CUBIC_FREE_BLOCK_COUNTS[k]
    if type(cert.bound) is not int or cert.bound != product:
        return False
    best, _ = _optimal_block_product(cert.extra_steps)
    return cert.bound == best


def _rational(x: int | Fraction, name: str) -> Fraction:
    if isinstance(x, bool) or type(x) not in (int, Fraction):
        raise ValueError(f"{name} must be a rational Fraction or integer")
    return Fraction(x)


def cubic_ground_population(side: int, q: int | Fraction,
                            minimum_contacts: int, *,
                            correction_lower: int | Fraction = 1,
                            correction_upper: int | Fraction = 1
                            ) -> CubicGroundPopulation:
    """Simultaneous high-contact event and structure/contact-map no-go.

    Let C be the EXACT all-H cubic ground-state count, K the number of
    nonbackbone contacts for every Hamiltonian cube, and W an upper bound on
    ALL cubic walks of length n. For event A_t={k>=t}, t<=K, every non-event
    state weighs <=u*q^(t-1), while C ground states each weigh >=l*q^K:

       P(A_t) >= l*C*q^K / (l*C*q^K + u*W*q^(t-1)).

    For any single ground fold, the C equal-energy ground paths imply
    P(fold)<=u/(u+l*(C-1)). For one complete filled-cube contact map,
    exactly 8 first-East geometries have that map, giving
    P(map)<=8*u/(8*u+l*(C-8)); neither claim uses high-q assumptions.
    """
    if type(side) is not int or side not in CUBE_GROUND_COUNTS:
        raise ValueError("audited side must be 2 or 3")
    if type(minimum_contacts) is not int:
        raise ValueError("event threshold must be an exact integer")
    weight, l, u = (_rational(q, "q"),
                    _rational(correction_lower, "correction_lower"),
                    _rational(correction_upper, "correction_upper"))
    if weight < 1 or not 0 < l <= u:
        raise ValueError("require exact q>=1 and 0<correction_lower<=correction_upper")
    n = side**3
    K = 2*n - 3*side**2 + 1
    if not 1 <= minimum_contacts <= K:
        raise ValueError("the event threshold must be in [1,K]")
    C = CUBE_GROUND_COUNTS[side]
    count = optimal_cubic_walk_cap(n - 2)
    if not verify_cubic_walk_cap(count):
        raise AssertionError("count certificate failed independent replay")
    gap = K - minimum_contacts + 1
    ratio = l*C*weight**gap / (l*C*weight**gap + u*count.bound)
    micro_upper = contact_map_confidence_upper(1, C,
                                               correction_lower=l,
                                               correction_upper=u)
    map_upper = contact_map_confidence_upper(8, C,
                                             correction_lower=l,
                                             correction_upper=u)
    return CubicGroundPopulation(
        side, weight, minimum_contacts, l, u, C, K, count,
        ratio, micro_upper, map_upper, C//8,
    )


def verify_cubic_ground_population(cert: CubicGroundPopulation) -> bool:
    """Fail-closed exact-rational replay, independent of float estimates."""
    if not isinstance(cert, CubicGroundPopulation):
        return False
    try:
        other = cubic_ground_population(cert.side, cert.q, cert.minimum_contacts,
                                         correction_lower=cert.correction_lower,
                                         correction_upper=cert.correction_upper)
    except (ValueError, AssertionError, TypeError, OverflowError):
        return False
    return cert == other
