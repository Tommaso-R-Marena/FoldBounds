"""Exact interval-box proof of conditional NMR / frozen-local-geometry no-go.

The independent checker exhaustively partitions continuous real-position
domains for each hydrogen. A box is cut only if interval geometry contradicts
a necessary shell/slab constraint, or its objective upper bound is already
below a declared cap. Every root is proved to enclose ALL admissible
displacements by a separate four-anchor cofactor inequality. Failure to
finish within the node budget is inconclusive, never accepted as a proof.
All arithmetic is exact Python int and Fraction; tick scale 1e-6 Å.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isqrt
from .noe_scaffold_dual import _dot, _sub, _ceil, _validate

SCALE=1000


@dataclass(frozen=True)
class EndpointBound:
    target: str
    cap_mA2: int
    root_mA: tuple
    quartet: tuple


@dataclass(frozen=True)
class BranchBoundProof:
    delta_mA: int
    h: EndpointBound
    hb: EndpointBound
    max_nodes_per_endpoint: int = 10000


@dataclass(frozen=True)
class BranchBoundResult:
    delta_mA: int
    objective_upper_mA2: int
    necessary_lower_mA2: int
    nodes_h: int
    nodes_hb: int
    eliminated_geometry_h: int
    eliminated_geometry_hb: int


def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])


def _root_valid(old,anchors,delta,spec):
    """Independent, exact cofactor proof that root contains ALL possibilities."""
    points=[anchors[k] for k in spec.quartet]
    row=[_sub(p,points[0]) for p in points[1:]]
    cof=(_cross(row[1],row[2]),_cross(row[2],row[0]),_cross(row[0],row[1]))
    det=_dot(row[0],cof[0])
    if det==0:return False
    e=[2*_ceil(_dot(_sub(old,p),_sub(old,p)))*delta+delta*delta for p in points]
    for j in range(3):
        bound=Fraction(sum(abs(cof[i][j])*(e[0]+e[i+1]) for i in range(3)),2*abs(det))
        lo,hi=spec.root_mA[j]
        if lo>hi or Fraction(-lo)<bound or Fraction(hi)<bound:
            return False
    return True


def _replay_one(old,objective,anchors,delta,spec,budget):
    d=delta*SCALE
    shells=[];safety={}
    for name,anchor in anchors.items():
        diff=tuple((old[i]-anchor[i])*SCALE for i in range(3))
        s=_dot(diff,diff)
        lower=max(0,isqrt(s)-d)**2
        upper=(_ceil(s)+d)**2
        shells.append((diff,lower,upper))
        safety[name]=2*_ceil(s)*d+d*d
    slabs=[]
    for i,j in combinations(tuple(anchors),2):
        normal=tuple((anchors[i][k]-anchors[j][k])*SCALE for k in range(3))
        slabs.append((normal,safety[i]+safety[j]))
    def optimistic(box):
        return sum(objective[j]*(box[j][1] if objective[j]>=0 else box[j][0])
                   for j in range(3))
    def contradicted(box):
        for diff,lo,hi in shells:
            low=0;high=0
            for j in range(3):
                a=diff[j]+box[j][0]
                b=diff[j]+box[j][1]
                low += 0 if a<=0<=b else min(a*a,b*b)
                high += max(a*a,b*b)
            if low>hi or high<lo:return True
        for normal,limit in slabs:
            minp=sum(normal[j]*(box[j][0] if normal[j]>=0 else box[j][1])
                     for j in range(3))
            maxp=sum(normal[j]*(box[j][1] if normal[j]>=0 else box[j][0])
                     for j in range(3))
            if 2*minp>limit or 2*maxp< -limit:return True
        return False
    root=tuple(tuple(z*SCALE for z in side) for side in spec.root_mA)
    stack=[(root,0)]
    nodes=0;geometry=0
    while stack:
        box,depth=stack.pop()
        nodes+=1
        if nodes>budget:return None
        if optimistic(box)<=spec.cap_mA2*SCALE:continue
        if contradicted(box):
            geometry+=1
            continue
        spread=[abs(objective[j])*(box[j][1]-box[j][0]) for j in range(3)]
        axis=max(range(3),key=lambda j:spread[j])
        low,high=box[axis]
        # Do not claim a proof if the numeric box grid is too coarse.
        if high-low<=1 or depth>100:return None
        middle=(low+high)//2
        a=list(box);a[axis]=(low,middle)
        b=list(box);b[axis]=(middle,high)
        children=[tuple(a),tuple(b)]
        children.sort(key=optimistic)
        for child in children:stack.append((child,depth+1))
    return nodes,geometry


def verify_branch_bound_no_go(raw,proof):
    """Complete exact certificate: roots, all boxes, and global directional gap."""
    try:
        if (not _validate(raw) or type(proof) is not BranchBoundProof or
            type(proof.delta_mA) is not int or proof.delta_mA<0 or
            type(proof.max_nodes_per_endpoint) is not int or
            not 1<=proof.max_nodes_per_endpoint<=1000000 or
            type(proof.h) is not EndpointBound or type(proof.hb) is not EndpointBound):
            return None
        both=(proof.h,proof.hb)
        for s in both:
            if (type(s.target) is not str or type(s.cap_mA2) is not int or
                s.cap_mA2<0 or type(s.root_mA) is not tuple or len(s.root_mA)!=3 or
                type(s.quartet) is not tuple or len(s.quartet)!=4 or
                len(set(s.quartet))!=4 or
                any(type(q) is not tuple or len(q)!=2 or
                    any(type(z) is not int for z in q) for q in s.root_mA)):
                return None
        if proof.h.target!='H' or proof.hb.target!='HB':return None
        anchors=raw['fixed_heavy_atoms_val69']
        h=raw['target_hydrogens_val69']
        if any(name not in anchors for s in both for name in s.quartet):
            return None
        direction=_sub(h['HB'],h['H'])
        norm2=_dot(direction,direction)
        required=norm2-raw['processed_restraint_upper_mA']*_ceil(norm2)
        total=proof.h.cap_mA2+proof.hb.cap_mA2
        if required<=0 or total>=required:return None
        for s in both:
            if not _root_valid(h[s.target],anchors,proof.delta_mA,s):
                return None
        left=_replay_one(h['H'],direction,anchors,proof.delta_mA,proof.h,
                         proof.max_nodes_per_endpoint)
        right=_replay_one(h['HB'],tuple(-x for x in direction),anchors,
                          proof.delta_mA,proof.hb,proof.max_nodes_per_endpoint)
        if left is None or right is None:return None
        return BranchBoundResult(proof.delta_mA,total,required,left[0],right[0],
                                 left[1],right[1])
    except (ValueError,TypeError,KeyError,OverflowError,AttributeError,ZeroDivisionError):
        return None


def default_val69_branch_bound():
    return BranchBoundProof(
       199,
       EndpointBound('H',1835000,
                     ((-430,430),(-790,790),(-1700,1700)),
                     ('N','CA','CG1','CG2')),
       EndpointBound('HB',800000,
                     ((-1300,1300),(-630,630),(-1750,1750)),
                     ('CB','CG1','CG2','C')),
       10000)
