"""Frozen 5-fold source-only NMR upper-distance restraint holdout.

This evaluates whether metric closure of *training upper limits alone* can
imply coarse <=8 Å proximity for withheld atom pairs. The complete lower+
upper restraint problem is NOT assessed by this benchmark. Distance labels,
PDB coordinates, and withheld upper limits are never read by the proposer.
The heldout source is a processed NEF file, not a second independent assay.
"""
from __future__ import annotations
from dataclasses import dataclass
from heapq import heappush, heappop
from typing import Sequence

from .noe_restraint_audit import AtomRef, NOERestraint, verify_record
from .noe_metric_closure import ProximityWitness, verify_proximity_witness


@dataclass(frozen=True)
class HoldoutFold:
    fold: int
    train: int
    heldout: int
    connected: int
    proximity_le_8a: int
    inferred_upper_stronger_than_heldout: int
    entire_interval_guaranteed_by_upper_only: int


def frozen_fold(sample: str, restraint_id: int, num_folds: int = 5) -> int:
    if (type(sample) is not str or not sample or type(restraint_id) is not int
        or restraint_id <= 0 or type(num_folds) is not int or num_folds < 2):
        raise ValueError('invalid split')
    h = 2166136261
    for c in f'{sample}:{restraint_id}':
        h = ((h ^ ord(c)) * 16777619) & 0xffffffff
    return h % num_folds


def _propose(train: Sequence[NOERestraint], source: AtomRef,
             target: AtomRef) -> ProximityWitness | None:
    """Untrusted Dijkstra proposal, verified by an independent path replay."""
    graph: dict[AtomRef, list[tuple[AtomRef,int,int]]] = {}
    for r in train:
        graph.setdefault(r.left, []).append((r.right,r.upper_mA,r.identity))
        graph.setdefault(r.right, []).append((r.left,r.upper_mA,r.identity))
    if source not in graph or target not in graph:
        return None
    dist = {source:0}; previous={};queue=[(0,source.residue,source.atom,source)]
    while queue:
        d,_,_,u=heappop(queue)
        if d!=dist[u]:continue
        if u==target:break
        for v,cost,identity in graph[u]:
            nd=d+cost
            if v not in dist or nd<dist[v]:
                dist[v]=nd;previous[v]=(u,identity)
                heappush(queue,(nd,v.residue,v.atom,v))
    if target not in dist:return None
    atoms=[target];ids=[];now=target
    while now!=source:
        p,rid=previous[now];ids.append(rid);atoms.append(p);now=p
    atoms.reverse();ids.reverse()
    proposal=ProximityWitness(tuple(atoms),tuple(ids),dist[target])
    # In training, this is an *unmeasured* atom-pair query.
    if not verify_proximity_witness(train,proposal,
                 cutoff_mA=proposal.upper_bound_mA,require_unmeasured=True):
        raise AssertionError('invalid proposed source-only metric path')
    return proposal


def five_fold_audit(sample: str, records: Sequence[NOERestraint]) -> tuple[HoldoutFold,...]:
    source=tuple(records)
    if not source or any(not verify_record(r) for r in source):
        raise ValueError('invalid records')
    if len({r.identity for r in source})!=len(source):
        raise ValueError('repeated restraint IDs')
    result=[]
    for fold in range(5):
        training=tuple(r for r in source if frozen_fold(sample,r.identity)!=fold)
        heldout=tuple(r for r in source if frozen_fold(sample,r.identity)==fold)
        connected=near=strict=whole=0
        for r in heldout:
            w=_propose(training,r.left,r.right)
            if w is None:continue
            connected+=1
            if w.upper_bound_mA<=8000:near+=1
            if w.upper_bound_mA<r.upper_mA:strict+=1
            # The source-only upper model permits a collapsed 1D embedding
            # of all atoms, so cannot prove any strictly positive lower limit.
            if r.lower_mA==0 and w.upper_bound_mA<=r.upper_mA:whole+=1
        result.append(HoldoutFold(fold,len(training),len(heldout),connected,near,strict,whole))
    return tuple(result)
