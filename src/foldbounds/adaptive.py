"""Anytime conformational-certificate refinement, with untrusted search policies.

Each selected frontier is independently re-verified. Greedy policies affect
only work allocation; their predictions are NOT part of the trust boundary.
All results remain finite 2D HP toy-model enclosures, not experimental data.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from .lattice import Certificate, Path, VerifiedBounds, refine_frontier
from .lookahead import verify_lookahead, subtree
from .lattice import _possible_hh_pairs

Policy = Literal['breadth', 'mass', 'probability']


@dataclass(frozen=True)
class RefinementSnapshot:
    step: int
    chosen: Path | None
    bounded_nodes: int
    exact_paths: int
    partition_lower: Fraction
    partition_upper: Fraction
    probability_lower: Fraction
    probability_upper: Fraction

    @property
    def probability_width(self) -> Fraction:
        return self.probability_upper-self.probability_lower


@dataclass(frozen=True)
class RefinementRun:
    certificate: Certificate
    bounds: VerifiedBounds
    history: tuple[RefinementSnapshot, ...]
    policy: Policy


def adaptive_refine(cert: Certificate, *, steps: int, lookahead: int = 0,
                    use_matching: bool = True,
                    policy: Policy = 'probability') -> RefinementRun:
    """Refine at most steps nodes of a sound complete prefix frontier.

    Breadth selects the shortest canonical prefix. Mass selects greatest
    uncertainty in event+complement weight; probability evaluates every
    legal one-step refinement and chooses smallest resulting event-probability
    interval, breaking ties with partition-width and canonical path order.
    This greedy scoring can be expensive: no runtime dominance is promised.
    """
    if type(steps) is not int or steps < 0:
        raise ValueError('steps must be a nonnegative integer')
    if policy not in ('breadth','mass','probability'):
        raise ValueError('invalid policy')
    current=verify_lookahead(cert,lookahead=lookahead,use_matching=use_matching)

    def snapshot(step: int, chosen: Path | None) -> RefinementSnapshot:
        return RefinementSnapshot(step,chosen,current.bounded_subtrees,
                                  current.exact_paths,current.partition_lower,
                                  current.partition_upper,current.probability_lower,
                                  current.probability_upper)

    history=[snapshot(0,None)]
    pairs=_possible_hh_pairs(cert.sequence)
    for step in range(1,steps+1):
        available=[node.path for node in cert.frontier if node.kind=='bound']
        if not available:
            break
        if policy=='breadth':
            chosen=min(available,key=lambda p:(len(p),p))
            next_cert=refine_frontier(cert,chosen)
            next_bounds=verify_lookahead(next_cert,lookahead=lookahead,use_matching=use_matching)
        elif policy=='mass':
            def uncertainty(path: Path) -> Fraction:
                m=subtree(cert.sequence,path,pairs,cert.q,cert.event_at_least,
                          lookahead,use_matching)
                return (m.e_hi-m.e_lo)+(m.n_hi-m.n_lo)
            chosen=min(available,key=lambda p:(-uncertainty(p),len(p),p))
            next_cert=refine_frontier(cert,chosen)
            next_bounds=verify_lookahead(next_cert,lookahead=lookahead,use_matching=use_matching)
        else:
            options=[]
            for path in available:
                proposed=refine_frontier(cert,path)
                bound=verify_lookahead(proposed,lookahead=lookahead,use_matching=use_matching)
                key=(bound.probability_upper-bound.probability_lower,
                     bound.partition_upper-bound.partition_lower,len(path),path)
                options.append((key,path,proposed,bound))
            _,chosen,next_cert,next_bounds=min(options,key=lambda item:item[0])
        if not (current.partition_lower<=next_bounds.partition_lower and
                next_bounds.partition_upper<=current.partition_upper and
                current.probability_lower<=next_bounds.probability_lower and
                next_bounds.probability_upper<=current.probability_upper):
            raise AssertionError('refinement violated the verified nesting invariant')
        cert,current=next_cert,next_bounds
        history.append(snapshot(step,chosen))
    return RefinementRun(cert,current,tuple(history),policy)
