"""Untrusted reproducible proposal generator; exact checker is authoritative.

Floating-point alternating pair projection supplies CANDIDATE coordinate
repairs only. The independent integer-only verifier is mandatory and can
reject any incorrect, inconsistent, or overoptimistic proposed result.
"""
from __future__ import annotations
from math import sqrt, isqrt
from pathlib import Path
import json
import sys

from foldbounds.noe_restraint_audit import NOERestraint,load_pinned_reference,verify_record
from foldbounds.noe_minimax_repair import (
    RepairCertificate,SCALE,verify_repair)


def _lower_violation(r: NOERestraint,radius_ticks:int)->bool:
    s=sum((x-y)**2 for x,y in zip(r.left_position_mA,r.right_position_mA))
    above=(r.upper_mA*SCALE+2*radius_ticks)**2 < s*SCALE*SCALE
    below_limit=r.lower_mA*SCALE-2*radius_ticks
    below=below_limit>0 and s*SCALE*SCALE<below_limit*below_limit
    return above or below


def propose_repair(sample: str,records:tuple[NOERestraint,...],*,
                   max_sweeps:int=120)->RepairCertificate:
    """Untrusted floating proposer followed by mandatory exact replay."""
    if type(sample) is not str or not sample or not records or type(max_sweeps) is not int or max_sweeps<=0:
        raise ValueError('invalid sample or budget')
    nodes={};ids=set()
    for r in records:
        if not verify_record(r) or r.identity in ids:
            raise ValueError('invalid or duplicate source record')
        ids.add(r.identity)
        for key,pos in ((r.left,r.left_position_mA),(r.right,r.right_position_mA)):
            if key in nodes and nodes[key]!=pos:
                raise ValueError('conflicting source coordinates')
            nodes[key]=pos
    lower_candidates=[r for r in records if _lower_violation(r,0)]
    if not lower_candidates:
        raise ValueError('all source restraints already satisfied: this demo needs a strict lower witness')
    for margin in (0.001,0.003,0.01,0.05):
        positions={k:list(map(float,x)) for k,x in nodes.items()}
        for _ in range(max_sweeps):
            n=0
            for r in records:
                left,right=positions[r.left],positions[r.right]
                v=[a-b for a,b in zip(left,right)]
                d=sqrt(sum(x*x for x in v))
                if d==0:
                    raise ValueError('degenerate pair: a geometric proposal unavailable')
                target=(r.lower_mA+margin if d<r.lower_mA else
                        r.upper_mA-margin if d>r.upper_mA else None)
                if target is None:continue
                eta=(target-d)/(2*d)
                for j in range(3):
                    left[j]+=v[j]*eta
                    right[j]-=v[j]*eta
                n+=1
            if not n:break
        moves=[]
        for atom,source in nodes.items():
            delta=tuple(round(positions[atom][j]*SCALE)-source[j]*SCALE for j in range(3))
            if delta!=(0,0,0):moves.append((atom,delta))
        max_sq=max((sum(x*x for x in m) for _,m in moves),default=0)
        upper=isqrt(max_sq)+(isqrt(max_sq)**2<max_sq)
        if upper<=0:continue
        best=(0,None)
        for r in lower_candidates:
            lo,hi=0,upper
            while lo+1<hi:
                mid=(lo+hi)//2
                if _lower_violation(r,mid):lo=mid
                else:hi=mid
            if _lower_violation(r,lo) and lo>best[0]:best=(lo,r.identity)
        if best[1] is None or upper<=best[0]:continue
        cert=RepairCertificate(sample,len(records),best[1],best[0],upper,tuple(moves))
        if verify_repair(records,cert) is not None:
            return cert
    raise RuntimeError('ABSTAIN: no fully checked proposal within the attempt budget')


if __name__=='__main__':
    root=Path(__file__).resolve().parents[1]
    ref=load_pinned_reference(root/'data/pinned_unambiguous_noe_2png_6nbn.json')
    for name,records in ref.items():
        cert=propose_repair(name,records)
        check=verify_repair(records,cert)
        assert check is not None
        print(name,check)
