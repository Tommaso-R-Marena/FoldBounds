"""Error-correcting structural fingerprints for authenticated GA95/GB95 NMR models.

A residue-disjoint N-/C-terminal bipartite K-matching maximizes the WORST
pair-distance gap across every deposited 2KDL/2KDM NMR conformer. All
distance comparisons use exact integer milli-Angstrom squared lengths,
with conservative integer square-root enclosing intervals.

The constructive K-matching is a classification witness. Its optimality
has a separate independently checkable certificate: a vertex cover of
size <K covering ALL edges with strictly greater margin (Konig's
theorem), so no strictly stronger K-matching can exist. No floating LP.

If K=2f+1 channels are used, a majority vote remains correct after up
to f arbitrary erroneous binary channels, provided every uncorrupted
readout has pair-distance error strictly < gap/2 and is drawn from the
DECLARED deposited NMR conformational envelope. This does NOT validate
dye placement, new mutants, thermal populations or biological response.
"""
from __future__ import annotations
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from typing import Mapping

from .ga_gb_nmr import GA95, GB95, MUTATIONS, parse_pdb

EXPECTED = {
    "2KDL": "fb75cc55eb8b89e0cec5e2be69901b348e60916375b9c410a0573edfca1db591",
    "2KDM": "1c49e5faa25b443285432d935c621151294ed3b636b0106e6b1a7f7ccf2d00ec",
}

def _milli(point):
    out = []
    for x in point:
        v = Decimal(str(x)) * 1000
        if v != v.to_integral_value():
            raise ValueError("PDB coordinates must occupy exact 0.001 Angstrom grid")
        out.append(int(v))
    return tuple(out)

def _d2(a,b):
    return sum((x-y)**2 for x,y in zip(_milli(a),_milli(b)))

def _ceil_sqrt(n):
    floor = isqrt(n)
    return floor if floor*floor == n else floor+1

@dataclass(frozen=True)
class Reporter:
    a: int
    b: int
    gap_milliang: int
    shorter: str
    ga_range_milliang: tuple[int,int]
    gb_range_milliang: tuple[int,int]

    @property
    def endpoints(self):
        return (self.a,self.b)

def load_real_models(ga:bytes,gb:bytes):
    for entry,raw in (("2KDL",ga),("2KDM",gb)):
        if sha256(raw).hexdigest()!=EXPECTED[entry]:
            raise ValueError(f"{entry}: unrecognized coordinate bytes; fail closed")
    a=parse_pdb(ga,"2KDL",GA95)
    b=parse_pdb(gb,"2KDM",GB95)
    if len(a)!=20 or len(b)!=20:
        raise ValueError("require 20 NMR conformers of each fold")
    return a,b

def extract_reporters(a,b, *, left_end=14, min_sep=4,
                      distance_min_ang=12, distance_max_ang=40):
    """Conservative integer bands enclose all C-alpha Euclidean distances."""
    if type(left_end) is not int or not 3<=left_end<=50:
        raise ValueError("invalid N-terminal division")
    if type(min_sep) is not int or not 1<=min_sep<=56:
        raise ValueError("invalid sequence separation")
    if not (type(distance_min_ang) is int and type(distance_max_ang) is int
            and 0<distance_min_ang<distance_max_ang):
        raise ValueError("assay range must be in positive exact integer Angstroms")
    if len(a)!=20 or len(b)!=20 or any(m.seq!=GA95 for m in a) or any(m.seq!=GB95 for m in b):
        raise ValueError("real sequence identity and 20-model constraint violated")
    result=[]
    for i in range(1,left_end+1):
        if i in MUTATIONS:continue
        for j in range(max(left_end+1,i+min_sep),57):
            if j in MUTATIONS:continue
            da=[_d2(m.ca[i],m.ca[j]) for m in a]
            db=[_d2(m.ca[i],m.ca[j]) for m in b]
            a_lo,a_hi=isqrt(min(da)),_ceil_sqrt(max(da))
            b_lo,b_hi=isqrt(min(db)),_ceil_sqrt(max(db))
            if min(a_lo,b_lo)<distance_min_ang*1000 or max(a_hi,b_hi)>distance_max_ang*1000:
                continue
            if a_hi<b_lo:
                gap=b_lo-a_hi;shorter="GA95"
            elif b_hi<a_lo:
                gap=a_lo-b_hi;shorter="GB95"
            else:
                continue
            result.append(Reporter(i,j,gap,shorter,(a_lo,a_hi),(b_lo,b_hi)))
    return tuple(sorted(result,key=lambda e:(-e.gap_milliang,e.a,e.b)))

def matching(edges,threshold:int):
    """Finite exact bipartite maximum matching using augmenting paths."""
    eligible=[e for e in edges if e.gap_milliang>=threshold]
    adj={}
    for e in eligible:adj.setdefault(e.a,[]).append(e.b)
    by_pair={(e.a,e.b):e for e in eligible}
    matched_right={}
    def augment(a,seen):
        for b in adj.get(a,()):
            if b in seen:continue
            seen.add(b)
            if b not in matched_right or augment(matched_right[b],seen):
                matched_right[b]=a
                return True
        return False
    for a in sorted(adj):
        augment(a,set())
    return tuple(sorted((by_pair[(a,b)] for b,a in matched_right.items()),
                        key=lambda e:(e.a,e.b)))

def cover_certificate(edges,strict_margin):
    """Small endpoint cover proving NO stronger K-matching exists."""
    qualified=[e for e in edges if e.gap_milliang>strict_margin]
    adj={}
    for e in qualified:adj.setdefault(e.a,set()).add(e.b)
    m=matching(qualified,strict_margin+1)
    a_to_b={e.a:e.b for e in m}
    b_to_a={e.b:e.a for e in m}
    all_a=set(adj)
    reach_a=set(all_a-set(a_to_b))
    reach_b=set()
    queue=list(sorted(reach_a))
    while queue:
        a=queue.pop()
        for b in adj.get(a,()):
            if a_to_b.get(a)==b or b in reach_b:continue
            reach_b.add(b)
            if b in b_to_a and b_to_a[b] not in reach_a:
                reach_a.add(b_to_a[b]);queue.append(b_to_a[b])
    left=tuple(sorted(all_a-reach_a))
    right=tuple(sorted(reach_b))
    if any(e.a not in left and e.b not in right for e in qualified):
        raise AssertionError("incomplete vertex cover")
    if len(left)+len(right)!=len(m):
        raise AssertionError("matching-cover equality failed")
    return left,right

@dataclass(frozen=True)
class CodeCertificate:
    reporters:tuple[Reporter,...]
    target_channels:int
    corrupted_channels:int
    bottleneck_milliang:int
    stronger_graph_cover_left:tuple[int,...]
    stronger_graph_cover_right:tuple[int,...]
    all_candidate_count:int
    coordinate_hashes:Mapping[str,str]
    left_end:int
    source_models_per_fold:int
    classifier_error_milliang_strict_bound:Fraction

def certify_code(a,b, *, channels=5,left_end=14):
    if type(channels) is not int or channels<1 or channels%2!=1:
        raise ValueError("positive odd channel count required")
    edges=extract_reporters(a,b,left_end=left_end)
    optimal=None
    for margin in sorted({e.gap_milliang for e in edges},reverse=True):
        chosen=matching(edges,margin)
        if len(chosen)>=channels:
            optimal=(margin,tuple(sorted(chosen[:channels],key=lambda e:(e.a,e.b))))
            break
    if optimal is None:
        raise ValueError(f"{channels} disjoint positive-gap reporters unavailable")
    margin,chosen=optimal
    left,right=cover_certificate(edges,margin)
    if len(left)+len(right)>=channels:
        raise AssertionError("upper-cover optimality gate failed")
    if len({p for e in chosen for p in e.endpoints})!=2*channels:
        raise AssertionError("reporter endpoints reused")
    if min(e.gap_milliang for e in chosen)!=margin:
        raise AssertionError("claimed bottleneck not attained")
    return CodeCertificate(chosen,channels,(channels-1)//2,margin,left,right,
                           len(edges),dict(EXPECTED),left_end,20,Fraction(margin,2))

def verify_code(cert:CodeCertificate,a,b):
    """Recompute every candidate and check the witness and cover, not solver output."""
    if not isinstance(cert,CodeCertificate):raise ValueError("bad certificate")
    edges=extract_reporters(a,b,left_end=cert.left_end)
    by_pair={(e.a,e.b):e for e in edges}
    if cert.source_models_per_fold!=20 or cert.coordinate_hashes!=EXPECTED:
        raise ValueError("untrusted structure provenance")
    if len(cert.reporters)!=cert.target_channels or cert.target_channels%2!=1:
        raise ValueError("bad correction code length")
    if cert.corrupted_channels!=(cert.target_channels-1)//2:
        raise ValueError("bad fault tolerance claim")
    if len({p for e in cert.reporters for p in e.endpoints})!=2*cert.target_channels:
        raise ValueError("overlapping residue reporter sites")
    if any(by_pair.get(e.endpoints)!=e for e in cert.reporters):
        raise ValueError("unsupported reporter bands")
    if min(e.gap_milliang for e in cert.reporters)!=cert.bottleneck_milliang:
        raise ValueError("invalid minimum interval gap")
    if cert.classifier_error_milliang_strict_bound!=Fraction(cert.bottleneck_milliang,2):
        raise ValueError("invalid deterministic readout error allowance")
    left=set(cert.stronger_graph_cover_left);right=set(cert.stronger_graph_cover_right)
    if (len(left)+len(right)>=cert.target_channels or
            len(left)!=len(cert.stronger_graph_cover_left) or
            len(right)!=len(cert.stronger_graph_cover_right)):
        raise ValueError("unsound optimality cover size")
    if any(e.a not in left and e.b not in right for e in edges
           if e.gap_milliang>cert.bottleneck_milliang):
        raise ValueError("stronger reporter uncovered")
    if cert.all_candidate_count!=len(edges):
        raise ValueError("candidate count mismatch")
    return True
