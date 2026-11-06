"""Mutation-invariant, proof-carrying NMR geometry codes with erasures.

Training only: SHA-256-authenticated 2KDL GA95 and 2KDM GB95 deposited
C-alpha coordinates. Reject every candidate label site changed anywhere
in the predeclared GA95/GB95/GA98/GB98/T25I/L20A mutation cascade.

For K residue-disjoint reporter pairs, if at most e sites cause channel
erasures and at most t remaining channels are arbitrarily wrong, then
ordinary non-tied majority is correct whenever 2*t+e<K and all other
reporter votes are correct. We do not pretend to know which remaining
readouts are erroneous or whether their physical measurement noise
satisfies the original model envelope.

The K-matching has a checkable strict-margin vertex-cover obstruction.
All contact-distance bands and thresholds remain exact rational or
integer milli-Angstrom quantities. This is retrospective experimental
coordinate classification, NOT a protein sequence-to-fold predictor
or a calibrated FRET experiment.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from math import isqrt
from typing import Mapping

from .ga_gb_error_codes import (
    EXPECTED, CodeCertificate, Reporter, extract_reporters, matching,
    cover_certificate, _d2, _ceil_sqrt
)
from .ga_gb_nmr import GA95, GB95

FORBIDDEN_MUTATION_SITES=frozenset({20,25,30,45})
TARGET_LENGTHS=(1,3,5,7,9)


@dataclass(frozen=True)
class MutationSafeCode:
    """A fully checkable reference matching certificate on safe sites."""
    reference: CodeCertificate
    forbidden_sites: frozenset[int]
    eligible_edges: int


def _safe_edges(training_a,training_b):
    base=extract_reporters(training_a,training_b,left_end=14)
    return tuple(edge for edge in base
                 if edge.a not in FORBIDDEN_MUTATION_SITES
                 and edge.b not in FORBIDDEN_MUTATION_SITES)


def make_safe_code(a,b,channels:int)->MutationSafeCode:
    if type(channels) is not int or channels<1 or channels%2!=1:
        raise ValueError("channels must be a positive odd integer")
    edges=_safe_edges(a,b)
    optimum=None
    for threshold in sorted({edge.gap_milliang for edge in edges},reverse=True):
        witnessed=matching(edges,threshold)
        if len(witnessed)>=channels:
            optimum=(threshold,tuple(sorted(witnessed[:channels],
                                            key=lambda x:x.endpoints)))
            break
    if optimum is None:
        raise ValueError(f"mutation-site-free {channels}-channel code is unavailable")
    gap,selected=optimum
    left,right=cover_certificate(edges,gap)
    c=CodeCertificate(
        selected,channels,(channels-1)//2,gap,left,right,
        len(edges),dict(EXPECTED),14,len(a),Fraction(gap,2)
    )
    result=MutationSafeCode(c,FORBIDDEN_MUTATION_SITES,len(edges))
    check_safe_code(result,a,b)
    return result


def check_safe_code(safe:MutationSafeCode,a,b)->bool:
    """Untrusted certificate checker; does not call make_safe_code."""
    if not isinstance(safe,MutationSafeCode):
        raise ValueError("unknown matching certificate")
    c=safe.reference
    edges=_safe_edges(a,b)
    by_pair={e.endpoints:e for e in edges}
    k=c.target_channels
    if (safe.forbidden_sites!=FORBIDDEN_MUTATION_SITES
            or c.coordinate_hashes!=EXPECTED or type(k) is not int
            or k<1 or k%2!=1 or c.left_end!=14
            or c.source_models_per_fold!=len(a) or len(a)!=len(b)):
        raise ValueError("model identity or site policy mismatch")
    if c.corrupted_channels!=(k-1)//2:
        raise ValueError("invalid code correction capacity")
    if (len(c.reporters)!=k
            or len({p for e in c.reporters for p in e.endpoints})!=2*k
            or any(p in FORBIDDEN_MUTATION_SITES
                   for e in c.reporters for p in e.endpoints)):
        raise ValueError("mutated/duplicate label sites")
    if (not c.reporters
            or any(by_pair.get(e.endpoints)!=e for e in c.reporters)
            or min(e.gap_milliang for e in c.reporters)!=c.bottleneck_milliang
            or c.classifier_error_milliang_strict_bound!=Fraction(c.bottleneck_milliang,2)
            or c.all_candidate_count!=len(edges) or safe.eligible_edges!=len(edges)):
        raise ValueError("wrong or forged reporter distance bands")
    left,right=c.stronger_graph_cover_left,c.stronger_graph_cover_right
    if (len(set(left))!=len(left) or len(set(right))!=len(right)
            or len(left)+len(right)>=k):
        raise ValueError("invalid strict-optimality cover size")
    if any(e.a not in left and e.b not in right
           for e in edges if e.gap_milliang>c.bottleneck_milliang):
        raise ValueError("uncovered stronger candidate edge")
    return True


def reporter_midpoint(edge:Reporter)->Fraction:
    if edge.shorter=="GA95":
        return Fraction(edge.ga_range_milliang[1]+edge.gb_range_milliang[0],2)
    if edge.shorter=="GB95":
        return Fraction(edge.gb_range_milliang[1]+edge.ga_range_milliang[0],2)
    raise ValueError("unknown shorter fold")


def channel_vote(edge:Reporter, value:Fraction)->str|None:
    if not isinstance(value,Fraction) or value<0:
        raise ValueError("measurement must be a nonnegative exact Fraction")
    mid=reporter_midpoint(edge)
    if value==mid:
        return None
    return edge.shorter if value<mid else ("GA95" if edge.shorter=="GB95" else "GB95")


def decode_with_erasures(safe:MutationSafeCode,
                        observations:Mapping[tuple[int,int],Fraction|None],
                        *, missing_sites:frozenset[int]=frozenset())->str|None:
    """No guess when channels tie or are ambiguous; invalid keys fail closed.

    Missing sites indicate an unusable tag at that residue, making its
    reporter an erasure; other supplied None values are channel erasures.
    Missing reporter entries are also erasures. Any explicit numeric
    reading on a missing site is ignored, not trusted.
    """
    if not isinstance(safe,MutationSafeCode) or not isinstance(observations,Mapping):
        raise ValueError("invalid code or observation mapping")
    if not isinstance(missing_sites,frozenset) or any(type(i) is not int for i in missing_sites):
        raise ValueError("missing sites must be an exact integer frozenset")
    by_pair={e.endpoints:e for e in safe.reference.reporters}
    if set(observations)-set(by_pair):
        raise ValueError("unrecognized reporter sites")
    correct={"GA95":0,"GB95":0}
    for pair,e in by_pair.items():
        if set(pair)&missing_sites:
            continue
        reading=observations.get(pair)
        if reading is None:
            continue
        vote=channel_vote(e,reading)
        if vote is not None:
            correct[vote]+=1
    if correct["GA95"]==correct["GB95"]:
        return None
    return "GA95" if correct["GA95"]>correct["GB95"] else "GB95"


def certified_budget(channels:int,erased_sites:int,wrong_readouts:int)->bool:
    """Exact soundness criterion: errors counted twice, erasures once."""
    if any(type(v) is not int for v in (channels,erased_sites,wrong_readouts)):
        raise ValueError("integer budgets only")
    if channels<1 or erased_sites<0 or wrong_readouts<0:
        raise ValueError("negative budgets prohibited")
    if erased_sites+wrong_readouts>channels:
        return False
    return 2*wrong_readouts+erased_sites<channels


def coordinate_distance_readings(code:MutationSafeCode,model):
    """Conservative PDB C-alpha readout, abstaining at sqrt ambiguity."""
    values={}
    for e in code.reference.reporters:
        square=_d2(model.ca[e.a],model.ca[e.b])
        lower,upper=isqrt(square),_ceil_sqrt(square)
        threshold=reporter_midpoint(e)
        if lower<=threshold<=upper:
            values[e.endpoints]=None
        else:
            values[e.endpoints]=Fraction(lower)
    return values


def finite_margin_to_truth(code:MutationSafeCode,model,truth:str)->Fraction:
    """Exact conservative, truth-labeled minimum margin; can be negative."""
    if truth not in ("GA95","GB95"):
        raise ValueError("unknown experimental fold label")
    out=None
    for e in code.reference.reporters:
        d2=_d2(model.ca[e.a],model.ca[e.b])
        lo,hi=isqrt(d2),_ceil_sqrt(d2)
        mid=reporter_midpoint(e)
        signed=mid-hi if e.shorter==truth else lo-mid
        out=signed if out is None else min(out,signed)
    if out is None:raise ValueError("empty reporter code")
    return out
