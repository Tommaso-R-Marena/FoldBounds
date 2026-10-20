"""Independent, model-aware π-ring geometry and fold-contact network separation.

This is purely structural analysis, NOT a π-interaction energy, predictive
fold-affinity estimator, or a thermodynamic population calculation.
"""
from __future__ import annotations
from collections import Counter
from itertools import combinations,product
from math import sqrt,dist
from statistics import median
from .ga_gb_nmr import NMRModel,contact_pairs,aromatic_centroid

PAIRS=((29,30),(29,45),(29,52),(30,45),(30,52),(45,52))

def _sub(a,b):return tuple(x-y for x,y in zip(a,b))
def _dot(a,b):return sum(x*y for x,y in zip(a,b))
def _cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def _norm(a):return sqrt(_dot(a,a))

def ring_normal(m:NMRModel,i:int):
    """Orientation-independent plane normal from carbon ring branches.

    Returns None for non-aromatics, absent atoms, or degenerate ring geometry.
    """
    if m.seq[i-1] not in 'FYW':return None
    a=m.sidechains.get(i,{})
    if not all(k in a for k in ('CG','CD1','CD2')):return None
    u=_sub(a['CD1'],a['CG']);v=_sub(a['CD2'],a['CG'])
    n=_cross(u,v);length=_norm(n)
    if length<1e-6:return None
    return tuple(x/length for x in n)

def classify_ring_pair(m:NMRModel,i:int,j:int)->dict|None:
    x=aromatic_centroid(m,i);y=aromatic_centroid(m,j)
    n1=ring_normal(m,i);n2=ring_normal(m,j)
    if x is None or y is None or n1 is None or n2 is None:return None
    r=_sub(y,x);d=_norm(r)
    alignment=min(1.0,abs(_dot(n1,n2)))
    axial=abs(_dot(r,n1));lateral=sqrt(max(0,d*d-axial*axial))
    # Purely geometric cutoffs; close rings without favorable orientation are
    # NOT automatically classified as stacked or energetically favorable.
    if d<=7.0 and alignment>=0.866025403784439 and 2.5<=axial<=5.5 and lateral<=4.5:
        label='parallel_plane_like'
    elif d<=7.0 and alignment<=0.5:
        label='perpendicular_plane_like'
    elif d<=7.0:
        label='close_other_orientation'
    else:
        label='separated'
    return {'centroid_distance_A':d,'abs_plane_normal_dot':alignment,
            'axial_ring_plane_separation_A':axial,'lateral_ring_offset_A':lateral,
            'orientation_class':label}

def _quantiles(items):
    x=sorted(items)
    if not x:raise ValueError('empty quantile series')
    return {'min':x[0],'median':median(x),'max':x[-1]}

def aromatic_summary(models:tuple[NMRModel,...])->list[dict]:
    out=[]
    if not models:raise ValueError('no models')
    for i,j in PAIRS:
        p=[classify_ring_pair(m,i,j) for m in models]
        real=[x for x in p if x is not None]
        c=Counter(x['orientation_class'] for x in real)
        out.append({'sites':[i,j],'valid_models':len(real),
                    'mean_centroid_distance_A':sum(x['centroid_distance_A'] for x in real)/len(real) if real else None,
                    'distance_range':_quantiles([x['centroid_distance_A'] for x in real]) if real else None,
                    'absolute_plane_normal_dot_range':_quantiles([x['abs_plane_normal_dot'] for x in real]) if real else None,
                    'orientation_class_counts':dict(sorted(c.items()))})
    return out

def jaccard(a:set,b:set)->float:
    d=len(a|b)
    return len(a&b)/d if d else 1.0

def network_discrimination(ga:tuple[NMRModel,...],gb:tuple[NMRModel,...])->dict:
    """Within-vs-between deposited-NMR-model contact-network Jaccard.

    Data are STRUCTURAL models of two already known folds; 400 cross-model
    comparisons are NOT prospective sequence-to-fold predictive validation.
    """
    if len(ga)<2 or len(gb)<2:raise ValueError('need independent models for both structures')
    amap=[contact_pairs(m) for m in ga];bmap=[contact_pairs(m) for m in gb]
    in_a=[jaccard(amap[i],amap[j]) for i,j in combinations(range(len(ga)),2)]
    in_b=[jaccard(bmap[i],bmap[j]) for i,j in combinations(range(len(gb)),2)]
    cross=[jaccard(x,y) for x,y in product(amap,bmap)]
    return {'within_GA95':_quantiles(in_a),'within_GB95':_quantiles(in_b),'between_folds':_quantiles(cross),
            'within_GA95_pairs':len(in_a),'within_GB95_pairs':len(in_b),'cross_fold_pairs':len(cross),
            'strict_model_level_separation':max(cross)<min(min(in_a),min(in_b)),
            'interpretation':'Retrospective deposited NMR contact-map separability, not mutation fold prediction.'}

def analyze_advanced(ga:tuple[NMRModel,...],gb:tuple[NMRModel,...])->dict:
    return {'GA95_aromatic_geometry':aromatic_summary(ga),
            'GB95_aromatic_geometry':aromatic_summary(gb),
            'contact_map_model_discrimination':network_discrimination(ga,gb)}