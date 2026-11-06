"""Exact sparse-distance reconstruction and bounded-noise contact assurance.

Coordinate models use integer milliangstroms. Three noncollinear anchors
plus three distances and an orientation sign per remaining site reconstruct
exactly. Four noncoplanar fixed anchors plus four distances per remaining
site give an exact *linear* worst-case error bound. Established Euclidean
trilateration, not a protein-folding predictor or a new mathematics claim.
"""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import isqrt

def _sub(a,b):
    return tuple(x-y for x,y in zip(a,b))
def _dot(a,b):
    return sum(x*y for x,y in zip(a,b))
def _cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def _sq(a,b):
    return _dot(_sub(a,b),_sub(a,b))
def _check_point(p):
    return isinstance(p,tuple) and len(p)==3 and all(type(x) is int for x in p)
def point_digest(p):
    return sha256("|".join(f"{x},{y},{z}" for x,y,z in p).encode("ascii")).hexdigest()

def select_anchors(points):
    p=tuple(points)
    if len(p)<4 or not all(_check_point(x) for x in p):
        raise ValueError("need four integer three-dimensional sites")
    a=min(range(len(p)),key=lambda i:(p[i],i))
    b=max((i for i in range(len(p)) if i!=a),
          key=lambda i:(_sq(p[i],p[a]),-i))
    u=_sub(p[b],p[a])
    c=max((i for i in range(len(p)) if i not in (a,b)),
          key=lambda i:(_dot(_cross(u,_sub(p[i],p[a])),
                             _cross(u,_sub(p[i],p[a]))),-i))
    n=_cross(u,_sub(p[c],p[a]))
    d=max((i for i in range(len(p)) if i not in (a,b,c)),
          key=lambda i:(abs(_dot(n,_sub(p[i],p[a]))),-i))
    if _dot(n,_sub(p[d],p[a]))==0:
        raise ValueError("no affinely independent four anchors")
    return (a,b,c,d)

def _basis4(a):
    v=tuple(_sub(x,a[0]) for x in a[1:])
    cof=(_cross(v[1],v[2]),_cross(v[2],v[0]),_cross(v[0],v[1]))
    det=_dot(v[0],cof[0])
    if det==0:raise ValueError("coplanar anchors")
    norm=max(sum(abs(c[j]) for c in cof) for j in range(3))
    return v,cof,det,norm

@dataclass(frozen=True)
class SparseDistanceCertificate:
    count:int
    anchor_indices:tuple[int,int,int,int]
    anchors:tuple
    squared_to_anchors:tuple
    recovered_sha256:str

def make_certificate(points,anchors=None):
    p=tuple(points)
    if len(p)<4 or not all(_check_point(x) for x in p):
        raise ValueError("bad input coordinates")
    if anchors is None:anchors=select_anchors(p)
    if (not isinstance(anchors,tuple) or len(anchors)!=4 or
        len(set(anchors))!=4 or any(type(i) is not int or not 0<=i<len(p) for i in anchors)):
        raise ValueError("invalid anchor indices")
    aa=tuple(p[i] for i in anchors)
    _basis4(aa)
    selected=set(anchors)
    rows=tuple(None if i in selected else tuple(_sq(x,a) for a in aa)
               for i,x in enumerate(p))
    return SparseDistanceCertificate(len(p),anchors,aa,rows,point_digest(p))

def reconstruct_exact(c):
    if type(c) is not SparseDistanceCertificate:
        raise ValueError("wrong certificate type")
    ids,aa,n=c.anchor_indices,c.anchors,c.count
    if (type(n) is not int or n<4 or not isinstance(ids,tuple) or len(ids)!=4 or
        len(set(ids))!=4 or any(type(i) is not int or not 0<=i<n for i in ids) or
        not isinstance(aa,tuple) or len(aa)!=4 or not all(_check_point(a) for a in aa) or
        not isinstance(c.squared_to_anchors,tuple) or len(c.squared_to_anchors)!=n or
        not isinstance(c.recovered_sha256,str) or len(c.recovered_sha256)!=64):
        raise ValueError("malformed certificate")
    v,cof,det,_=_basis4(aa)
    positions=[None]*n
    for i,a in zip(ids,aa):
        if c.squared_to_anchors[i] is not None:raise ValueError("anchor has distance row")
        positions[i]=a
    for i,row in enumerate(c.squared_to_anchors):
        if i in ids:continue
        if not isinstance(row,tuple) or len(row)!=4 or any(type(s) is not int or s<0 for s in row):
            raise ValueError("bad distances")
        rhs=tuple(_dot(v[k],v[k])+row[0]-row[k+1] for k in range(3))
        nums=tuple(sum(cof[k][j]*rhs[k] for k in range(3)) for j in range(3))
        if any(x%(2*det) for x in nums):raise ValueError("no integer reconstruction")
        q=tuple(aa[0][j]+nums[j]//(2*det) for j in range(3))
        if any(_sq(q,a)!=d for a,d in zip(aa,row)):raise ValueError("inconsistent distances")
        positions[i]=q
    ans=tuple(positions)
    if point_digest(ans)!=c.recovered_sha256:raise ValueError("source digest mismatch")
    return ans

def verify_certificate(c):
    try:reconstruct_exact(c);return True
    except (ValueError,TypeError,IndexError,OverflowError,AttributeError):return False

@dataclass(frozen=True)
class ContactRobustness:
    all_pairs:int
    guaranteed_contacts:int
    guaranteed_noncontacts:int
    undecided_pairs:int
    nominal_contacts:int
    max_coordinate_error_mA:Fraction
    max_pair_distance_error_mA:Fraction
    anchor_inverse_infty:Fraction

def robust_contacts(c,*,cutoff_mA=8000,epsilon_mA=10,min_sequence_separation=4):
    if (any(type(v) is not int for v in (cutoff_mA,epsilon_mA,min_sequence_separation))
        or cutoff_mA<=0 or epsilon_mA<0 or min_sequence_separation<1):
        raise ValueError("invalid thresholds")
    p=reconstruct_exact(c)
    _,_,det,norm=_basis4(c.anchors)
    denom=abs(det);hs=[0]*len(p);selected=set(c.anchor_indices)
    for i,row in enumerate(c.squared_to_anchors):
        if i in selected:continue
        dist=max(isqrt(s)+(isqrt(s)**2<s) for s in row)
        eta=2*dist*epsilon_mA+epsilon_mA**2
        hs[i]=eta*norm
    yes=no=maybe=nominal=total=0
    rr=cutoff_mA*denom
    for i in range(len(p)):
        for j in range(i+min_sequence_separation,len(p)):
            total+=1
            s=_sq(p[i],p[j]);nominal+=int(s<=cutoff_mA**2)
            err=2*(hs[i]+hs[j])
            if rr>=err and s*denom**2<=(rr-err)**2:yes+=1
            elif s*denom**2>(rr+err)**2:no+=1
            else:maybe+=1
    h=max(hs)
    return ContactRobustness(total,yes,no,maybe,nominal,
                             Fraction(h,denom),Fraction(4*h,denom),
                             Fraction(norm,denom))

@dataclass(frozen=True)
class SignedThreeAnchorCertificate:
    count:int
    anchor_indices:tuple[int,int,int]
    anchors:tuple
    squared_to_anchors:tuple
    orientation_signs:tuple
    recovered_sha256:str

def make_signed_certificate(points,anchors=None):
    p=tuple(points)
    if len(p)<4 or not all(_check_point(x) for x in p):
        raise ValueError("bad coordinates")
    if anchors is None:anchors=select_anchors(p)[:3]
    if (not isinstance(anchors,tuple) or len(anchors)!=3 or len(set(anchors))!=3 or
        any(type(i) is not int or not 0<=i<len(p) for i in anchors)):
        raise ValueError("bad anchor indices")
    a=tuple(p[i] for i in anchors)
    N=_cross(_sub(a[1],a[0]),_sub(a[2],a[0]))
    if _dot(N,N)==0:raise ValueError("collinear anchors")
    excluded=set(anchors)
    rows=tuple(None if i in excluded else tuple(_sq(v,w) for w in a)
               for i,v in enumerate(p))
    signs=[]
    for i,v in enumerate(p):
        V=0 if i in excluded else _dot(N,_sub(v,a[0]))
        signs.append((V>0)-(V<0))
    return SignedThreeAnchorCertificate(len(p),anchors,a,rows,tuple(signs),point_digest(p))

def reconstruct_signed_exact(c):
    if type(c) is not SignedThreeAnchorCertificate:raise ValueError("wrong certificate type")
    n,ids,a=c.count,c.anchor_indices,c.anchors
    if (type(n) is not int or n<4 or not isinstance(ids,tuple) or len(ids)!=3 or
        len(set(ids))!=3 or any(type(i) is not int or not 0<=i<n for i in ids) or
        not isinstance(a,tuple) or len(a)!=3 or not all(_check_point(x) for x in a) or
        not isinstance(c.squared_to_anchors,tuple) or len(c.squared_to_anchors)!=n or
        not isinstance(c.orientation_signs,tuple) or len(c.orientation_signs)!=n or
        not isinstance(c.recovered_sha256,str) or len(c.recovered_sha256)!=64):
        raise ValueError("malformed signed certificate")
    u,v=_sub(a[1],a[0]),_sub(a[2],a[0])
    N=_cross(u,v);area=_dot(N,N)
    if not area:raise ValueError("degenerate anchors")
    g11,g22,g12=_dot(u,u),_dot(v,v),_dot(u,v)
    cof=(_cross(v,N),_cross(N,u),N)
    out=[None]*n
    for i,pos in zip(ids,a):
        if c.squared_to_anchors[i] is not None or c.orientation_signs[i]!=0:
            raise ValueError("bad anchor row")
        out[i]=pos
    for i,row in enumerate(c.squared_to_anchors):
        if i in ids:continue
        sign=c.orientation_signs[i]
        if (not isinstance(row,tuple) or len(row)!=3 or
            any(type(s) is not int or s<0 for s in row) or
            type(sign) is not int or sign not in (-1,0,1)):
            raise ValueError("bad signed distances")
        p1=g11+row[0]-row[1];p2=g22+row[0]-row[2]
        if p1%2 or p2%2:raise ValueError("squared-distance parity")
        t1,t2=p1//2,p2//2
        volume_sq=area*row[0]-(g22*t1*t1+g11*t2*t2-2*g12*t1*t2)
        if volume_sq<0:raise ValueError("negative volume square")
        V=isqrt(volume_sq)
        if V*V!=volume_sq or (V==0 and sign!=0) or (V>0 and sign==0):
            raise ValueError("not a square or invalid chirality")
        V*=sign
        nums=tuple(cof[0][j]*t1+cof[1][j]*t2+cof[2][j]*V for j in range(3))
        if any(x%area for x in nums):raise ValueError("nonintegral point")
        q=tuple(a[0][j]+nums[j]//area for j in range(3))
        if any(_sq(q,t)!=s for t,s in zip(a,row)) or _dot(N,_sub(q,a[0]))!=V:
            raise ValueError("inconsistent distances and sign")
        out[i]=q
    ans=tuple(out)
    if point_digest(ans)!=c.recovered_sha256:raise ValueError("wrong source digest")
    return ans

def verify_signed_certificate(c):
    try:reconstruct_signed_exact(c);return True
    except (ValueError,TypeError,IndexError,OverflowError,AttributeError):return False
