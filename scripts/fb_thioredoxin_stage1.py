#!/usr/bin/env python3
"""Reference-only topology/mechanics models; no target labels are accepted.

Weighted generalized ANM (a **known** class of model, Zheng 2008):
  H_a = sum_contacts (7/r)^2 [(1-a)u*u.T + a*I3]
        tensored with [[1,-1],[-1,1]], plus tau I.
At a=0 this exactly recovers the original inverse-square weighted ANM.
At a=1, this is an isotropic distance-weighted *ridge* GNM variant,
NOT the original unweighted pseudoinverse GNM comparator.

A deliberately simple secondary rank-fusion predictor mixes frozen pure
GNM resistance and directional weighted ANM ranks 3:1. This is NOT a novel
physical theory. It is included as a transparency/ablation comparator.
"""
from __future__ import annotations
import hashlib
import json
import numpy as np

CUTOFF_SQ=100.0
RIDGE=0.2
POWER=2
SCALE=7.0
ALPHA=0.5  # fixed reference-only physical interpolation parameter
BLEND_GNM=0.75  # secondary rank fusion only


def contacts(x, cutoff_sq=CUTOFF_SQ):
    x=np.asarray(x,dtype=np.float64)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all():
        raise ValueError('bad C-alpha coordinate matrix')
    edges=[]
    for i in range(len(x)):
        for j in range(i+1,len(x)):
            v=x[i]-x[j]; d2=float(v@v)
            if d2<=0:raise ValueError('coincident C-alpha residues')
            if d2<=cutoff_sq:edges.append((i,j,v,d2))
    if len(edges)<len(x)-1:raise ValueError('contact count too low')
    return edges


def rankdata(arr):
    x=np.asarray(arr,float)
    if not np.isfinite(x).all():raise ValueError('invalid score rank input')
    order=np.argsort(x,kind='stable')
    out=np.empty(len(x),dtype=float)
    k=0
    while k<len(order):
        j=k+1
        while j<len(order) and x[order[k]]==x[order[j]]:j+=1
        out[order[k:j]]=0.5*(k+j-1)+1
        k=j
    return out


def scalar_gnm(ref,ii,jj,edges=None):
    ref=np.asarray(ref,float)
    n=len(ref);edges=contacts(ref) if edges is None else edges
    L=np.zeros((n,n),dtype=np.float64)
    for i,j,v,d2 in edges:
        L[i,i]+=1;L[j,j]+=1;L[i,j]-=1;L[j,i]-=1
    # Fully specified scalar model; choose anchor using deterministic BFS
    adjacency=[set() for _ in range(n)]
    for i,j,_,_ in edges:adjacency[i].add(j);adjacency[j].add(i)
    dist=np.full((n,n),np.inf)
    for start in range(n):
        dist[start,start]=0;queue=[start]
        for curr in queue:
            for nb in sorted(adjacency[curr]):
                if not np.isfinite(dist[start,nb]):
                    dist[start,nb]=dist[start,curr]+1
                    queue.append(nb)
    if not np.isfinite(dist).all():raise ValueError('GNM network not connected')
    anchor=int(np.argmin(dist.sum(axis=1)))
    ids=[i for i in range(n) if i!=anchor];inv=np.linalg.inv(L[np.ix_(ids,ids)])
    ind=np.full(n,-1,int)
    for m,i in enumerate(ids):ind[i]=m
    scores=[]
    for i,j in zip(ii,jj):
        a,b=int(ind[i]),int(ind[j]);
        q=(inv[a,a] if a!=-1 else 0)+(inv[b,b] if b!=-1 else 0)-2*(inv[a,b] if a!=-1 and b!=-1 else 0)
        scores.append(q)
    return np.asarray(scores,float), {'anchor':anchor,'edges':len(edges)}


def block_anm(ref,ii,jj,alpha,weighted=True,edges=None):
    ref=np.asarray(ref,float)
    if not (0 <= alpha <= 1):raise ValueError('anisotropy mixture must be in [0,1]')
    n=len(ref);edges=contacts(ref) if edges is None else edges
    H=np.zeros((n*3,n*3),dtype=np.float64)
    for i,j,v,d2 in edges:
        u=v/np.sqrt(d2)
        weight=(SCALE**2/d2 if weighted else 1)
        K=weight*((1-alpha)*np.outer(u,u)+alpha*np.eye(3))
        a=slice(3*i,3*i+3);b=slice(3*j,3*j+3)
        H[a,a]+=K;H[b,b]+=K;H[a,b]-=K;H[b,a]-=K
    H.flat[::3*n+1]+=RIDGE
    # Only the proposal uses floating-point; exact GNM certification separate.
    cov=np.linalg.inv(H)
    r=ref[ii]-ref[jj]
    d=np.linalg.norm(r,axis=1)
    if np.any(d<=0):raise ValueError('coincident scored pair')
    u=r/d[:,None]
    blocks=cov.reshape(n,3,n,3)
    C=blocks[ii,:,ii,:]+blocks[jj,:,jj,:]-blocks[ii,:,jj,:]-blocks[jj,:,ii,:]
    score=np.einsum('pi,pij,pj->p',u,C,u)
    if not np.isfinite(score).all():raise ValueError('invalid inverse response')
    return score


def frozen_methods(ref,ii,jj):
    ref=np.asarray(ref,float);ii=np.asarray(ii,int);jj=np.asarray(jj,int)
    edges=contacts(ref)
    gnm,meta=scalar_gnm(ref,ii,jj,edges)
    weighted=block_anm(ref,ii,jj,0,True,edges)
    plain=block_anm(ref,ii,jj,0,False,edges)
    hybrid=block_anm(ref,ii,jj,ALPHA,True,edges)
    deg=np.zeros(len(ref),int)
    for i,j,_,_ in edges:deg[i]+=1;deg[j]+=1
    d=np.linalg.norm(ref[ii]-ref[jj],axis=1)
    fusion=BLEND_GNM*rankdata(gnm)/len(gnm)+(1-BLEND_GNM)*rankdata(weighted)/len(gnm)
    return {
        'GNM_resistance':gnm,
        'ANM_weighted_inverse_square':weighted,
        'ANM_ridge':plain,
        'packing_inverse':1/(deg[ii]+1)+1/(deg[jj]+1),
        'reference_distance':d,
        'sequence_distance':(jj-ii).astype(float),
        'GANM_half_isotropic':hybrid,
        'GNM_ANM_rank_fusion_75_25':fusion,
    },meta

# ---------- Independent blinded stage-one PyMOL adapter ----------
import csv
from pathlib import Path
import re
import zipfile

PROTOCOL_VERSION='FOLDBOUNDS_THIO_2TRX_STAGE1_20261105_V1'
REFERENCE_PDB='2TRX'
SEALED_TARGET='1XOA'  # DO NOT FETCH OR READ ANY TARGET DATA IN THIS STAGE.
MIN_SEQUENCE_SEPARATION=4
DISTANCE_RANGE=(5.,16.)
ORIGINAL_BASELINES=('GNM_resistance','ANM_weighted_inverse_square',
                   'ANM_ridge','packing_inverse','reference_distance','sequence_distance')
MIN_REFERENCE_RESIDUES=45
MIN_REFERENCE_PAIRS=100


def _digest(payload):return hashlib.sha256(payload).hexdigest()


def _source_contents():
    # PyMOL replaces the module-global __file__ with its own initialization path
    # under some run() modes. The function code object retains its origin.
    actual=Path(fb_thio_freeze.__code__.co_filename).resolve()
    if not actual.is_file() or actual.suffix.lower()!='.py':
        raise RuntimeError('Cannot establish true PyMOL source script path')
    blob=actual.read_bytes()
    if b'FOLDBOUNDS_THIO_2TRX_STAGE1_20261105_V1' not in blob or b'def fb_thio_freeze' not in blob:
        raise RuntimeError('Detected invalid PyMOL source file; no results produced')
    return blob


def _read_reference(cmd,obj):
    model=cmd.get_model('(%s) and chain A and name CA'%obj,state=1)
    coords={};names={}
    for atom in model.atom:
        label=str(atom.resi).strip()
        if not re.fullmatch(r'[0-9]+',label):continue
        k=int(label)
        if k in coords:
            raise ValueError('Duplicate C-alpha at residue %d; abstain'%k)
        coords[k]=[float(c) for c in atom.coord]
        names[k]=str(atom.resn).upper()
    residues=np.array(sorted(coords),dtype=np.int64)
    if len(residues)<MIN_REFERENCE_RESIDUES:
        raise ValueError('Only %d reference C-alphas: abstain'%len(residues))
    xyz=np.asarray([coords[int(k)] for k in residues],dtype=float)
    if not np.isfinite(xyz).all():raise ValueError('nonfinite reference coordinates')
    return residues, np.asarray([names[int(k)] for k in residues],dtype='U4'), xyz


def _eligible_pairs(x):
    ii=[];jj=[]
    for i in range(1,len(x)-1):
        for j in range(i+MIN_SEQUENCE_SEPARATION,len(x)-1):
            d=float(np.linalg.norm(x[i]-x[j]))
            if DISTANCE_RANGE[0]<=d<=DISTANCE_RANGE[1]:
                ii.append(i);jj.append(j)
    if len(ii)<MIN_REFERENCE_PAIRS:raise ValueError('only %d eligible reference pairs'%len(ii))
    return np.asarray(ii,dtype=np.int64),np.asarray(jj,dtype=np.int64)


def _frozen_protocol():
    return dict(
        version=PROTOCOL_VERSION,
        reference='2TRX crystal chain A state 1', sealed_target='1XOA NMR all conformers',
        prohibition='NEVER read/fetch target in stage one',
        primary='GANM_half_isotropic', secondary='GNM_ANM_rank_fusion_75_25',
        alpha_isotropic=ALPHA, inverse_square_contact_weights=True,
        contact_A=10.0,ridge=RIDGE,
        rank_fusion_GNM_weight=BLEND_GNM,
        baselines=list(ORIGINAL_BASELINES),
        eligible_pair_rule='exclude terminal CAs; at least 4 sequence indices apart; reference distance 5 to 16 A inclusive',
        strict_acceptance='primary Spearman >=0.50 AND primary minus best original baseline >=0.15',
        target_matching='same numbered and named residues; at least 85 percent pair retention; all target conformers >=10',
        limitations=['known G-ANM family, no new force law',
          'reference-only predictions do not validate real dynamics',
          'NMR deposited structures are refinement ensembles, not physical time series',
          'model settings frozen BEFORE access to target coordinates',
          '2TRX and 1XOA redox match should be verified; mismatch triggers abstention'],
        research_history='All prior ubiquitin, GB1, BPTI, CspB, crambin and RNase A targets are disclosed development data',
    )


def fb_thio_freeze():
    from pymol import cmd
    print('\n=== FOLDBOUNDS THIOREDOXIN 2TRX: REFERENCE-ONLY STAGE 1 ===')
    print('Sealed target 1XOA will NOT be accessed. Model alpha=0.5 fixed.')
    source=_source_contents()
    obj='fb_thio_2trx_reference'
    if obj in cmd.get_names('objects'):cmd.delete(obj)
    cmd.fetch(REFERENCE_PDB.lower(),name=obj,type='pdb',async_=0)
    if obj not in cmd.get_names('objects') or cmd.count_atoms(obj)<100:
        raise RuntimeError('2TRX reference could not be fetched')
    residues,names,xyz=_read_reference(cmd,obj)
    ii,jj=_eligible_pairs(xyz)
    methods,diagnostics=frozen_methods(xyz,ii,jj)
    if not all(len(v)==len(ii) and np.isfinite(v).all() for v in methods.values()):
        raise ValueError('invalid model score arrays; abstain')
    protocol=_frozen_protocol()
    src_hash=_digest(source)
    protocol['code_sha256']=src_hash
    output_root=Path.home()/'Downloads'
    if not output_root.is_dir():output_root=Path.home()
    parent=output_root/('FoldBounds_Thio_Stage1_'+PROTOCOL_VERSION)
    folder=parent;attempt=2
    while folder.exists():
        folder=Path(str(parent)+'_'+str(attempt));attempt+=1
    folder.mkdir(parents=True)
    npz=folder/'FROZEN_THIO_REFERENCE.npz'
    np.savez_compressed(npz,residues=residues,names=names,reference_coords=xyz,
        pair_i=ii,pair_j=jj,
        **{'prediction_'+k:np.asarray(v,float) for k,v in methods.items()})
    (folder/'FROZEN_PROTOCOL.json').write_text(json.dumps(protocol,indent=2,sort_keys=True)+'\n')
    (folder/'REFERENCE_DIAGNOSTICS.json').write_text(json.dumps(diagnostics,indent=2,sort_keys=True)+'\n')
    (folder/'fb_thioredoxin_stage1.py').write_bytes(source)
    (folder/'READ_ME_FIRST.txt').write_text(
        'This contains ONLY crystal reference 2TRX and frozen predictions. '\
        'DO NOT LOAD 1XOA until independent audit of this Stage1 ZIP.\n'\
        'Return entire Stage1 output archive to ChatGPT.\n',encoding='utf8')
    sums={p.name:_digest(p.read_bytes()) for p in sorted(folder.iterdir()) if p.is_file()}
    (folder/'SHA256SUMS.json').write_text(json.dumps(sums,indent=2,sort_keys=True)+'\n')
    destination=output_root/'FoldBounds_THIO_STAGE1_20261105_SEND_BACK.zip'
    if destination.exists():
        for number in range(2,1000):
            alternative=output_root/('FoldBounds_THIO_STAGE1_20261105_%02d_SEND_BACK.zip'%number)
            if not alternative.exists():destination=alternative;break
        else:raise RuntimeError('Could not name output ZIP safely')
    with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED,compresslevel=7) as z:
        for p in sorted(folder.iterdir()):z.write(p,p.name)
    checksum=_digest(destination.read_bytes())
    print('FROZEN:',len(residues),'CA residues,',len(ii),'eligible pairs,',diagnostics['edges'],'contacts')
    print('Source SHA256:',src_hash)
    print('Archive SHA256:',checksum)
    print('NO NMR TARGET ACCESSED. RETURN THE ENTIRE ZIP:',destination)
    return str(destination)

try:
    from pymol import cmd
    cmd.extend('fb_thio_freeze',fb_thio_freeze)
    print('FoldBounds thioredoxin Stage1 ready. Run: fb_thio_freeze')
except ImportError:
    pass

if __name__=='__main__':
    print('Load this Python file with PyMOL: run C:/.../fb_thioredoxin_stage1.py')