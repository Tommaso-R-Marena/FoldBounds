"""Source-only NOE upper-metric closure with independently replayed witnesses.

If d(u,v)<=U on every included restraint edge, triangle inequality proves
d(i,j)<=sum(path U). No deposited coordinates are used to infer a new pair.
Every result is conditional on simultaneous truth of the source upper limits.
"""
from __future__ import annotations
from dataclasses import dataclass
from heapq import heappop,heappush
from typing import Sequence
from .noe_restraint_audit import AtomRef,NOERestraint,verify_record


@dataclass(frozen=True)
class ProximityWitness:
    atoms: tuple[AtomRef, ...]
    restraint_ids: tuple[int, ...]
    upper_bound_mA: int

    @property
    def pair(self) -> tuple[AtomRef,AtomRef]:
        return (self.atoms[0],self.atoms[-1])

    @property
    def sequence_separation(self) -> int:
        return abs(self.atoms[0].residue-self.atoms[-1].residue)


def _graph(records:Sequence[NOERestraint]):
    adj,refs,observed={}, {}, set()
    for r in records:
        if not verify_record(r) or r.identity in refs:
            raise ValueError("invalid/duplicate source restraint")
        refs[r.identity]=r
        e=frozenset((r.left,r.right))
        if len(e)!=2 or e in observed:
            raise ValueError("duplicate atom-pair input")
        observed.add(e)
        adj.setdefault(r.left,[]).append((r.right,r.upper_mA,r.identity))
        adj.setdefault(r.right,[]).append((r.left,r.upper_mA,r.identity))
    return adj,refs,observed


def verify_proximity_witness(records:Sequence[NOERestraint],
                             witness:ProximityWitness, *, cutoff_mA:int=8000,
                             require_unmeasured:bool=True)->bool:
    """Independent path checker; does not trust Dijkstra or coordinates."""
    if (type(witness) is not ProximityWitness or type(cutoff_mA) is not int
        or cutoff_mA<0 or type(require_unmeasured) is not bool):
        return False
    try:
        _,refs,observed=_graph(records)
        atoms,ids=witness.atoms,witness.restraint_ids
        if (not isinstance(atoms,tuple) or not isinstance(ids,tuple)
            or len(atoms)<2 or len(ids)!=len(atoms)-1
            or atoms[0]==atoms[-1] or len(set(atoms))!=len(atoms)
            or any(type(a) is not AtomRef for a in atoms)
            or any(type(k) is not int for k in ids)):
            return False
        if require_unmeasured and frozenset((atoms[0],atoms[-1])) in observed:
            return False
        mass=0
        for left,right,key in zip(atoms,atoms[1:],ids):
            r=refs.get(key)
            if r is None or frozenset((r.left,r.right))!=frozenset((left,right)):
                return False
            mass+=r.upper_mA
        return (type(witness.upper_bound_mA) is int and
                mass==witness.upper_bound_mA and mass<=cutoff_mA)
    except (ValueError,TypeError,IndexError,AttributeError):
        return False


def discover_proximities(records:Sequence[NOERestraint], *,
                         cutoff_mA:int=8000)->tuple[ProximityWitness,...]:
    """Find new implied pairs and replay every source-only proof path."""
    if type(cutoff_mA) is not int or cutoff_mA<0:
        raise ValueError("invalid cutoff")
    records=tuple(records)
    if not records:
        return ()
    adj,_,observed=_graph(records)
    atoms=sorted(adj,key=lambda a:(a.residue,a.atom))
    index={a:i for i,a in enumerate(atoms)}
    result=[]
    for source_index,source in enumerate(atoms):
        dist={source:0};prev={}
        frontier=[(0,source_index,source)]
        while frontier:
            weight,_,at=heappop(frontier)
            if weight!=dist[at] or weight>cutoff_mA:
                continue
            for neighbor,upper,rid in adj[at]:
                d=weight+upper
                if d>cutoff_mA or (neighbor in dist and d>=dist[neighbor]):
                    continue
                dist[neighbor]=d
                prev[neighbor]=(at,rid)
                heappush(frontier,(d,index[neighbor],neighbor))
        for target in atoms[source_index+1:]:
            if target not in dist or frozenset((source,target)) in observed:
                continue
            path=[target];links=[];cur=target
            while cur!=source:
                before,rid=prev[cur]
                links.append(rid);path.append(before);cur=before
            path.reverse();links.reverse()
            witness=ProximityWitness(tuple(path),tuple(links),dist[target])
            if not verify_proximity_witness(records,witness,cutoff_mA=cutoff_mA):
                raise AssertionError("producer path rejected by checker")
            result.append(witness)
    return tuple(result)
