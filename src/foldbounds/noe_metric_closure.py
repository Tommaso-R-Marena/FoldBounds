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


@dataclass(frozen=True)
class UpperOnlySharpnessCertificate:
    """A 1D embedding satisfying every included upper-distance restraint.

    Together with a path upper bound of U, an endpoint separation equal to U
    proves U is the *best possible* universal bound from upper-only data.
    Lower-distance limits, stereochemistry and excluded restraints are not used.
    """
    coordinate_potential_mA: tuple[tuple[AtomRef,int], ...]


def propose_upper_only_sharpness(records:Sequence[NOERestraint],
                                 witness:ProximityWitness)->UpperOnlySharpnessCertificate:
    """Construct a source-distance potential; independently check the result."""
    records=tuple(records)
    if not verify_proximity_witness(records,witness,
                                    cutoff_mA=witness.upper_bound_mA):
        raise ValueError("invalid source path witness")
    adj,_,_=_graph(records)
    nodes=sorted(adj,key=lambda p:(p.residue,p.atom))
    source=witness.atoms[0]
    dist={source:0}
    index={p:i for i,p in enumerate(nodes)}
    queue=[(0,index[source],source)]
    while queue:
        d,_,u=heappop(queue)
        if dist[u]!=d:continue
        for v,upper,_ in adj[u]:
            nd=d+upper
            if v not in dist or nd<dist[v]:
                dist[v]=nd
                heappush(queue,(nd,index[v],v))
    candidate=UpperOnlySharpnessCertificate(tuple((p,dist.get(p,0)) for p in nodes))
    if not verify_upper_only_sharpness(records,witness,candidate):
        raise AssertionError("shortest-path dual potential rejected")
    return candidate


def verify_upper_only_sharpness(records:Sequence[NOERestraint],
                                witness:ProximityWitness,
                                cert:UpperOnlySharpnessCertificate)->bool:
    """Exact primal path + dual Lipschitz potential optimality check.

    Path of weight U proves all Euclidean embeddings satisfy d(s,t)<=U.
    Integer potential phi with phi(s)=0, phi(t)=U and
    |phi(u)-phi(v)|<=U_uv for every observed upper constraint explicitly
    realizes all upper restraints by 3D coordinates (phi(v),0,0). Hence no
    smaller universal bound is possible using these upper limits alone.
    """
    if type(cert) is not UpperOnlySharpnessCertificate:
        return False
    try:
        records=tuple(records)
        if not verify_proximity_witness(records,witness,
                                        cutoff_mA=witness.upper_bound_mA):
            return False
        adj,_,_=_graph(records)
        values=cert.coordinate_potential_mA
        if not isinstance(values,tuple) or len(values)!=len(adj):
            return False
        phi={}
        for item in values:
            if not isinstance(item,tuple) or len(item)!=2:
                return False
            atom,d=item
            if type(atom) is not AtomRef or type(d) is not int or atom in phi:
                return False
            phi[atom]=d
        if set(phi)!=set(adj):
            return False
        if (phi[witness.atoms[0]]!=0 or
            phi[witness.atoms[-1]]!=witness.upper_bound_mA):
            return False
        return all(abs(phi[r.left]-phi[r.right])<=r.upper_mA
                   for r in records)
    except (ValueError,TypeError,IndexError,AttributeError,KeyError):
        return False
