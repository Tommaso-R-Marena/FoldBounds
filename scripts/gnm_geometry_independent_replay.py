#!/usr/bin/env python3
"""Independent exact-integer replay of RNase A geometric contact no-goods.

This is a DIFFERENT implementation from the producer: no SciPy, no geometry
producer imports, no inverse generator, no NMR target. One pinned source ZIP,
all binary64 positions lifted to a common integer denominator, and rational
Farkas multipliers checked by exact arithmetic. Counts accepted abstract
binary assignments with a separately implemented component DFS.
"""
from __future__ import annotations
import argparse,hashlib,io,itertools,json,math,zipfile
from collections import defaultdict
from fractions import Fraction as Fr
from pathlib import Path
import numpy as np
SOURCE='ec077437f4bd4b17f3c4527891bada306a77aab1953794ec5b08d76ca8c86cff'
FORMAT='FOLDBOUNDS_REFERENCE_ONLY_GEOMETRIC_NOGOOD_V1'

def checked_source(archive_path):
    blob=Path(archive_path).read_bytes()
    if hashlib.sha256(blob).hexdigest()!=SOURCE:raise ValueError('source ZIP hash changed')
    with zipfile.ZipFile(io.BytesIO(blob)) as z:
        if z.testzip():raise ValueError('bad ZIP CRC')
        members=set(z.namelist());manifest=json.loads(z.read('SHA256SUMS.json'))
        if members-{'SHA256SUMS.json'}!=set(manifest):raise ValueError('source manifest coverage')
        for file,dig in manifest.items():
            if hashlib.sha256(z.read(file)).hexdigest()!=dig:raise ValueError('source manifest mismatch')
        npz=z.read('FROZEN_REFERENCE_ONLY.npz')
        with np.load(io.BytesIO(npz),allow_pickle=False) as data:
            points=np.asarray(data['reference_coords'],dtype=np.float64)
            if points.shape!=(124,3) or not np.isfinite(points).all():raise ValueError('source geometry')
            raw=[[float(x).as_integer_ratio() for x in row] for row in points]
    den=math.lcm(20,*(q for row in raw for _,q in row))
    xyz=[[p*(den//q) for p,q in row] for row in raw]
    return xyz,den,hashlib.sha256(npz).hexdigest()

def optional_edge_indices(xyz,den):
    # Completely integer geometry. Radius=10A; each coordinate perturbed by ±.05A.
    two=den//10;threshold=100*den*den
    uncertain=[];sure=0;nominal=0;possible=0
    for i in range(124):
        for j in range(i+1,124):
            a=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
            dmin=sum(max(0,u-two)**2 for u in a)
            dmax=sum((u+two)**2 for u in a)
            dnom=sum(u*u for u in a)
            on_all=dmax<=threshold
            may_on=dmin<=threshold
            sure+=on_all;possible+=may_on;nominal+=dnom<=threshold
            if may_on and not on_all:uncertain.append((i,j))
    if (sure,nominal,possible,len(uncertain))!=(972,1020,1080,108):
        raise ValueError('unexpected contact graph')
    return uncertain

def exact_linear_halfspace(edge,on,xyz,den):
    i,j=edge;r=[Fr(xyz[i][k]-xyz[j][k],den) for k in range(3)]
    sq=sum(x*x for x in r)
    sign=1 if on else -1
    vec={}
    for k in range(3):
        vec[(i,k)]=2*sign*r[k]
        vec[(j,k)]=-2*sign*r[k]
    # Quadratic remainder bound from |delta_ik-delta_jk| <=1/10.
    b=Fr(100)-sq if on else sq-Fr(100)+Fr(3,100)
    return vec,b

def strict_fraction(text):
    if type(text) is not str or len(text)>1024:raise ValueError('invalid rational literal')
    q=Fr(text)
    if str(q)!=text:raise ValueError('noncanonical rational literal')
    return q

def certificate_contradiction(source_coords,den,edges,states,weights):
    eqns=[exact_linear_halfspace(edge,bool(v),source_coords,den)
          for edge,v in zip(edges,states)]
    if len(eqns)!=len(weights) or not weights or not any(v>0 for v in weights) or any(v<0 for v in weights):
        raise ValueError('invalid Farkas multipliers')
    keys=set().union(*(eq[0].keys() for eq in eqns))
    # For every coordinate vector in [-1/20,+1/20]^3, this expression
    # cannot fall below minus (1/20) times its L1 coefficient norm.
    lower=-Fr(1,20)*sum(abs(sum(w*eqns[i][0].get(k,0) for i,w in enumerate(weights))) for k in keys)
    upper=sum(w*eqns[i][1] for i,w in enumerate(weights))
    return lower-upper

def independent_count(total_edges,forbidden):
    # Build components via DFS, not union-find. For each component brute-force
    # bit combinations; max size 7 in this particular source certificate.
    adjacency={i:set() for i in range(total_edges)}
    for edgeids,_ in forbidden:
        for i in edgeids:
            adjacency[i].update(edgeids)
    seen=set();components=[];count=1
    for start in range(total_edges):
        if start in seen:continue
        stack=[start];nodes=set()
        while stack:
            u=stack.pop()
            if u in nodes:continue
            nodes.add(u);stack.extend(adjacency[u]-nodes)
        seen.update(nodes);nodes=sorted(nodes)
        clauses=[(ids,states) for ids,states in forbidden if ids[0] in nodes]
        if len(nodes)>20:raise ValueError('unbounded component - refuse naive enumeration')
        good=0
        for bits in itertools.product((0,1),repeat=len(nodes)):
            assignment=dict(zip(nodes,bits))
            if all(any(assignment[e]!=v for e,v in zip(edgeids,states))
                   for edgeids,states in clauses):good+=1
        components.append({'edge_ids':nodes,'forbidden_clauses':len(clauses),
                           'valid_assignments':good,'total_assignments':2**len(nodes)})
        count*=good
    return count,sorted(components,key=lambda c:(-len(c['edge_ids']),c['edge_ids'][0]))

def replay(stage1,certificate):
    data=json.loads(Path(certificate).read_text());xyz,den,npzsha=checked_source(stage1)
    if (data.get('format')!=FORMAT or data.get('source_stage1_sha256')!=SOURCE or
        data.get('source_npz_sha256')!=npzsha or data.get('box_epsilon_A')!='1/20' or
        data.get('contact_cutoff_squared_A')!='100'):
        raise ValueError('metadata changed')
    edges=optional_edge_indices(xyz,den)
    if data.get('n_optional_edges')!=108 or data.get('optional_edges')!=[list(e) for e in edges]:
        raise ValueError('optional contacts incomplete')
    clauses=data.get('forbidden_patterns')
    if not isinstance(clauses,list) or len(clauses)!=20:raise ValueError('incomplete no-good clauses')
    forbidden=[];sizes=[];previous=None
    for cl in clauses:
        ids=cl.get('edge_ids');states=cl.get('forbidden_on_bits')
        if (not isinstance(ids,list) or len(ids) not in (2,3) or
            any(type(x) is not int for x in ids) or ids!=sorted(set(ids)) or
            ids[0]<0 or ids[-1]>=len(edges)):
            raise ValueError('invalid clause edge identifiers')
        if (not isinstance(states,list) or len(states)!=len(ids) or
            any(type(v) is not int or v not in (0,1) for v in states)):
            raise ValueError('invalid no-good truth assignment')
        key=(len(ids),ids,states)
        # JSON lists cannot be tuple-compared against same? fine.
        if previous is not None and key<=previous:raise ValueError('unsorted or duplicate clauses')
        previous=key
        weights=[strict_fraction(v) for v in cl.get('lambda',[])]
        margin=certificate_contradiction(xyz,den,[edges[i] for i in ids],states,weights)
        if margin<=0 or str(margin)!=cl.get('gap'):raise ValueError('Farkas witness not sound')
        forbidden.append((ids,states));sizes.append(len(ids))
    if (sizes.count(2),sizes.count(3))!=(8,12) or (
        (data.get('n_pairs_impossible'),data.get('n_triples_impossible'))!=(8,12)):
        raise ValueError('wrong pattern counts')
    # Every forbidden 3-pattern must be additional, not merely duplicate of a pair.
    pairs=[(tuple(ids),tuple(st)) for ids,st in forbidden if len(ids)==2]
    for ids,st in forbidden:
        if len(ids)==3:
            if any((tuple(ids[k] for k in (i,j)),(st[i],st[j])) in pairs
                   for i,j in ((0,1),(0,2),(1,2))):
                raise ValueError('redundant 3-edge nogood')
    count,components=independent_count(108,forbidden)
    if (data.get('components')!=components or
        data.get('independent_binary_assignments')!=str(2**108) or
        data.get('not_yet_excluded_assignments')!=str(count) or
        data.get('excluded_assignments')!=str(2**108-count) or
        data.get('excluded_fraction')!=str(Fr(2**108-count,2**108))):
        raise ValueError('forged count or component decomposition')
    return {'status':'PASS','source_only':True,'n_nogood_certificates':len(clauses),
            'exact_rational_margins':True,'at_least_excluded_fraction':str(Fr(2**108-count,2**108)),
            'excluded_percent':float(100*(1-Fr(count,2**108))),
            'allowed_binary_assignments_upper':count,
            'certificate_sha256':hashlib.sha256(Path(certificate).read_bytes()).hexdigest()}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage1');parser.add_argument('certificate');args=parser.parse_args()
    print(json.dumps(replay(args.stage1,args.certificate),indent=2))