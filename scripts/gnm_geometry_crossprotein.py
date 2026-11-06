#!/usr/bin/env python3
"""Source-only feasibility certificates for the two earlier frozen proteins.

Only CspB/Crambin_FROZEN_BEFORE_TARGET.npz is accessed. The NMR evaluation
members inside the same archive are deliberately never read.
"""
import argparse,hashlib,io,itertools,json,zipfile
from fractions import Fraction as F
from pathlib import Path
import numpy as np
from gnm_geometry_producer import (EPS,FORMAT,pair_witness,higher_order_witness,
    verify_witness,exact_component_count)

ZIP_SHA='666a5d734fbd61be0f182c954474841bb552b71b3a1e9156d6c21b8c2f095f95'
EXPECTED={'CspB':(67,509,543,561,52),'Crambin':(46,318,334,354,36)}

def frozen_coords(path,name):
    if name not in EXPECTED:raise ValueError('unsupported protein')
    if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=ZIP_SHA:
        raise ValueError('wrong frozen multi-protein archive')
    with zipfile.ZipFile(path) as archive:
        data=archive.read(f'{name}_FROZEN_BEFORE_TARGET.npz')
    with np.load(io.BytesIO(data),allow_pickle=False) as npz:
        coords=np.asarray(npz['reference_coords'],dtype=np.float64)
    if coords.shape!=(EXPECTED[name][0],3) or not np.isfinite(coords).all():
        raise ValueError('invalid reference')
    return coords,hashlib.sha256(data).hexdigest()

def edges_from_reference(coords):
    xyz=[[F.from_float(float(v)) for v in row] for row in coords]
    certain=[];nominal=[];possible=[]
    for i,j in itertools.combinations(range(len(coords)),2):
        r=[abs(xyz[i][k]-xyz[j][k]) for k in range(3)]
        dmin=sum(max(F(0),v-2*EPS)**2 for v in r)
        dmax=sum((v+2*EPS)**2 for v in r)
        dnom=sum(v*v for v in r)
        if dmax<=100:certain.append((i,j))
        if dnom<=100:nominal.append((i,j))
        if dmin<=100:possible.append((i,j))
    return xyz,certain,nominal,possible

def produce(path,name,outpath):
    coords,npzsha=frozen_coords(path,name)
    xyz,certain,nominal,possible=edges_from_reference(coords)
    opt=sorted(set(possible)-set(certain));expected=EXPECTED[name]
    if (len(certain),len(nominal),len(possible),len(opt))!=expected[1:]:
        raise ValueError('frozen source contact topology changed')
    clauses=[];pair_used=set()
    for i,j in itertools.combinations(range(len(opt)),2):
        if not (set(opt[i])&set(opt[j])):continue
        for flags in itertools.product((False,True),repeat=2):
            cert=pair_witness(opt[i],opt[j],*flags,xyz)
            if cert:
                weights,_=max(cert,key=lambda entry:entry[1])
                clauses.append((([i,j],list(map(int,flags))),weights))
                pair_used.add(((i,j),flags))
    pairs=len(clauses)
    for ids in itertools.combinations(range(len(opt)),3):
        if len(set(opt[ids[0]])|set(opt[ids[1]])|set(opt[ids[2]]))>4:continue
        for flags in itertools.product((False,True),repeat=3):
            if any(((ids[a],ids[b]),(flags[a],flags[b])) in pair_used
                   for a,b in ((0,1),(0,2),(1,2))):continue
            weights=higher_order_witness([opt[i] for i in ids],flags,xyz)
            if weights:clauses.append(((list(ids),list(map(int,flags))),weights))
    encoded=[]
    for (ids,flags),weights in clauses:
        encoded.append({'edge_ids':ids,'forbidden_on_bits':flags,
                        'lambda':[str(w) for w in weights],
                        'gap':str(verify_witness([opt[i] for i in ids],flags,weights,xyz))})
    encoded.sort(key=lambda c:(len(c['edge_ids']),c['edge_ids'],c['forbidden_on_bits']))
    allowed,components=exact_component_count(len(opt),[(c['edge_ids'],[bool(v) for v in c['forbidden_on_bits']]) for c in encoded])
    total=1<<len(opt)
    document={'format':'FOLDBOUNDS_FROZEN_CROSSPROTEIN_GEOMETRIC_NOGOOD_V1','protein':name,
              'original_frozen_archive_sha256':ZIP_SHA,'reference_npz_sha256':npzsha,
              'epsilon_A':'1/20','cutoff_squared_A':'100','n_residues':len(coords),
              'contact_counts':list(expected[1:4]),'optional_edges':[list(e) for e in opt],
              'forbidden_patterns':encoded,'pair_clauses':pairs,'triple_clauses':len(encoded)-pairs,
              'components':components,'abstract_binary_assignments':str(total),
              'not_excluded_assignments':str(allowed),
              'excluded_fraction':str(F(total-allowed,total)),
              'scope':'Frozen reference geometry only. No NMR target labels accessed.'}
    raw=(json.dumps(document,sort_keys=True,separators=(',',':'))+'\n').encode()
    Path(outpath).write_bytes(raw)
    return {'protein':name,'optional':len(opt),'pairs':pairs,'triples':len(encoded)-pairs,
            'excluded_percent':float(100*(1-F(allowed,total))),
            'sha256':hashlib.sha256(raw).hexdigest()}

if __name__=='__main__':
    pa=argparse.ArgumentParser();pa.add_argument('archive');pa.add_argument('protein',choices=EXPECTED);pa.add_argument('output');a=pa.parse_args()
    print(json.dumps(produce(a.archive,a.protein,a.output),indent=2))