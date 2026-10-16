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
    audit_frontiers: tuple[Certificate, ...]
    lookahead: int
    use_matching: bool
    use_transfer: bool
    witness_nodes: int
    use_monotone: bool
    transfer_memory: int
    use_turn_family: bool = False
    terminal_box: bool = False



def _intersect(a: VerifiedBounds, b: VerifiedBounds) -> VerifiedBounds:
    """Intersect independently verified sound rational mass enclosures."""
    elo,ehi=max(a.event_lower,b.event_lower),min(a.event_upper,b.event_upper)
    nlo,nhi=max(a.complement_lower,b.complement_lower),min(a.complement_upper,b.complement_upper)
    zlo=max(a.partition_lower,b.partition_lower,elo+nlo)
    zhi=min(a.partition_upper,b.partition_upper,ehi+nhi)
    p_lo=Fraction(1) if nhi==0 else (elo/(elo+nhi) if elo else Fraction(0))
    p_hi=Fraction(0) if ehi==0 else (ehi/(ehi+nlo) if nlo else Fraction(1))
    p_lo=max(a.probability_lower,b.probability_lower,p_lo)
    p_hi=min(a.probability_upper,b.probability_upper,p_hi)
    if not (0<=elo<=ehi and 0<=nlo<=nhi and 0<=zlo<=zhi
            and 0<=p_lo<=p_hi<=1):
        raise AssertionError('intersection of sound certificates is inconsistent')
    return VerifiedBounds(zlo,zhi,elo,ehi,nlo,nhi,p_lo,p_hi,
                          b.exact_paths,b.bounded_subtrees)


def verify_refinement_run(run: RefinementRun) -> VerifiedBounds:
    """Replay all frontier checks; recompute the cumulative enclosure.

    The audit trail, not merely the latest frontier, supports the strongest
    combined bounds. No values recorded in snapshots or run.bounds are trusted.
    """
    if (not run.audit_frontiers or run.certificate!=run.audit_frontiers[-1]
            or len(run.history)!=len(run.audit_frontiers)
            or run.history[0].step!=0 or run.history[0].chosen is not None):
        raise ValueError('invalid refinement audit trail')
    reference=run.audit_frontiers[0]
    identity=(reference.sequence,reference.q,reference.event_at_least)
    result=None
    for i,certificate in enumerate(run.audit_frontiers):
        if (certificate.sequence,certificate.q,certificate.event_at_least)!=identity:
            raise ValueError('audit trail mixes distinct thermodynamic problems')
        if i:
            chosen=run.history[i].chosen
            if (chosen is None or run.history[i].step!=i or
                refine_frontier(run.audit_frontiers[i-1],chosen)!=certificate):
                raise ValueError('audit trail includes a non-refinement transition')
        raw=verify_lookahead(certificate,lookahead=run.lookahead,
                            use_matching=run.use_matching,
                            use_transfer=run.use_transfer,
                            witness_nodes=run.witness_nodes,
                            use_monotone=run.use_monotone,
                            transfer_memory=run.transfer_memory,
                            use_turn_family=run.use_turn_family,
                            terminal_box=run.terminal_box)
        result=raw if result is None else _intersect(result,raw)
        h=run.history[i]
        if (h.partition_lower!=result.partition_lower or
            h.partition_upper!=result.partition_upper or
            h.probability_lower!=result.probability_lower or
            h.probability_upper!=result.probability_upper or
            h.exact_paths!=result.exact_paths or
            h.bounded_nodes!=result.bounded_subtrees):
            raise ValueError('stored snapshot disagrees with replayed certificate')
    if run.bounds!=result:
        raise ValueError('refinement result does not match checked audit trail')
    return result

def adaptive_refine(cert: Certificate, *, steps: int, lookahead: int = 0,
                    use_matching: bool = True, use_transfer: bool = False,
                    witness_nodes: int = 0, use_monotone: bool = False,
                    transfer_memory: int = 2,
                    use_turn_family: bool = False,
                    terminal_box: bool = False,
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
    current=verify_lookahead(cert,lookahead=lookahead,use_matching=use_matching,use_transfer=use_transfer,
                             witness_nodes=witness_nodes,use_monotone=use_monotone,
                             transfer_memory=transfer_memory,
                             use_turn_family=use_turn_family,
                             terminal_box=terminal_box)

    def snapshot(step: int, chosen: Path | None) -> RefinementSnapshot:
        return RefinementSnapshot(step,chosen,current.bounded_subtrees,
                                  current.exact_paths,current.partition_lower,
                                  current.partition_upper,current.probability_lower,
                                  current.probability_upper)

    history=[snapshot(0,None)]
    frontiers=[cert]
    pairs=_possible_hh_pairs(cert.sequence)
    for step in range(1,steps+1):
        available=[node.path for node in cert.frontier if node.kind=='bound']
        if not available:
            break
        if policy=='breadth':
            chosen=min(available,key=lambda p:(len(p),p))
            next_cert=refine_frontier(cert,chosen)
            next_bounds=verify_lookahead(next_cert,lookahead=lookahead,use_matching=use_matching,use_transfer=use_transfer,
                             witness_nodes=witness_nodes,use_monotone=use_monotone,
                             transfer_memory=transfer_memory,
                             use_turn_family=use_turn_family,
                             terminal_box=terminal_box)
        elif policy=='mass':
            def uncertainty(path: Path) -> Fraction:
                m=subtree(cert.sequence,path,pairs,cert.q,cert.event_at_least,
                          lookahead,use_matching,True,use_transfer,witness_nodes,use_monotone,transfer_memory,use_turn_family,terminal_box)
                return (m.e_hi-m.e_lo)+(m.n_hi-m.n_lo)
            chosen=min(available,key=lambda p:(-uncertainty(p),len(p),p))
            next_cert=refine_frontier(cert,chosen)
            next_bounds=verify_lookahead(next_cert,lookahead=lookahead,use_matching=use_matching,use_transfer=use_transfer,
                             witness_nodes=witness_nodes,use_monotone=use_monotone,
                             transfer_memory=transfer_memory,
                             use_turn_family=use_turn_family,
                             terminal_box=terminal_box)
        else:
            options=[]
            for path in available:
                proposed=refine_frontier(cert,path)
                bound=verify_lookahead(proposed,lookahead=lookahead,use_matching=use_matching,use_transfer=use_transfer,
                             witness_nodes=witness_nodes,use_monotone=use_monotone,
                             transfer_memory=transfer_memory,
                             use_turn_family=use_turn_family,
                             terminal_box=terminal_box)
                combined=_intersect(current,bound)
                key=(combined.probability_upper-combined.probability_lower,
                     combined.partition_upper-combined.partition_lower,len(path),path)
                options.append((key,path,proposed,bound))
            _,chosen,next_cert,next_bounds=min(options,key=lambda item:item[0])
        # A deeper subtree bound may be wider than its old parent bound.
        # Retain both independently verified certificates via intersection.
        cert,current=next_cert,_intersect(current,next_bounds)
        frontiers.append(cert)
        history.append(snapshot(step,chosen))
    return RefinementRun(cert,current,tuple(history),policy,tuple(frontiers),
                         lookahead,use_matching,use_transfer,witness_nodes,use_monotone,transfer_memory,use_turn_family,terminal_box)
