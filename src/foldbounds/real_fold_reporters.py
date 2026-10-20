"""Propose mutation-independent structural distance reporters for GA95/GB95.

Rank by WORST-CASE interval separation across all deposited NMR models,
not best-case mean difference; no trained model, estimated free energy,
experimental FRET prediction, or claim that chemical labeling is neutral.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import dist
from .ga_gb_nmr import GA95,GB95,MUTATIONS,NMRModel

@dataclass(frozen=True)
class DistanceReporter:
    i:int
    j:int
    aa_i:str
    aa_j:str
    GA95_min_A:float
    GA95_max_A:float
    GB95_min_A:float
    GB95_max_A:float
    worst_case_separation_A:float
    state_with_shorter_distance:str
    max_symmetric_CA_to_label_error_A:float

def reporter_candidates(a:tuple[NMRModel,...],b:tuple[NMRModel,...],
                        min_distance:float=12.0,max_distance:float=40.0,
                        min_seq_separation:int=4)->tuple[DistanceReporter,...]:
    """Only residues unchanged across GA95/GB95, with disjoint NMR CA ranges."""
    if not a or not b or any(x.seq!=GA95 for x in a) or any(x.seq!=GB95 for x in b):
        raise ValueError('real protein identity/mutation mismatch')
    n=len(GA95)
    excluded=set(MUTATIONS)
    out=[]
    for i in range(1,n+1):
        for j in range(i+min_seq_separation,n+1):
            if i in excluded or j in excluded:continue
            x=[dist(m.ca[i],m.ca[j]) for m in a]
            y=[dist(m.ca[i],m.ca[j]) for m in b]
            amin,amax=min(x),max(x)
            bmin,bmax=min(y),max(y)
            if min(amin,bmin)<min_distance or max(amax,bmax)>max_distance:continue
            if amax<bmin: gap=bmin-amax;shorter='GA95'
            elif bmax<amin:gap=amin-bmax;shorter='GB95'
            else:continue
            out.append(DistanceReporter(i,j,GA95[i-1],GA95[j-1],amin,amax,bmin,bmax,gap,shorter,gap/2))
    return tuple(sorted(out,key=lambda x:(-x.worst_case_separation_A,x.i,x.j)))

def fret_proxy(distance_A:float,forster_A:float)->float:
    if not 0<distance_A<1000 or not 0<forster_A<1000:
        raise ValueError('invalid lengths')
    return 1.0/(1.0+(distance_A/forster_A)**6)

def fret_ca_proxy_envelopes(reporter:DistanceReporter,forster_A:float=25.0,localization_A:float=0.0)->dict:
    if localization_A<0:raise ValueError('negative localization uncertainty')
    out={}
    for name in ('GA95','GB95'):
        lo=getattr(reporter,f'{name}_min_A');hi=getattr(reporter,f'{name}_max_A')
        low=max(1e-9,lo-localization_A);high=hi+localization_A
        # FRET transfer decreases with r.
        out[name]=[fret_proxy(high,forster_A),fret_proxy(low,forster_A)]
    return {'CA_distance_only_not_actual_dye_efficiency':out,
            'forster_radius_A_assumed':forster_A,
            'localization_bound_A_assumed':localization_A}

def _frac(x):
    if type(x) not in (int,Fraction):raise TypeError('use rational fractions, not unbounded floating measurements')
    return Fraction(x)

def invert_two_state_emission(y_lo,y_hi, a_lo,a_hi,b_lo,b_hi)->tuple[Fraction,Fraction]:
    """Sharp possible population of fold A from bounded two-state mixture.

    Assumptions: exactly A and B states (no unfolded state); emission y is
    a linear ensemble mixture; A emits more than B without overlapping bands;
    measured y in [y_lo,y_hi]. Other molecules/states invalidate these bounds.
    """
    ylo,yhi,alo,ahi,blo,bhi=map(_frac,(y_lo,y_hi,a_lo,a_hi,b_lo,b_hi))
    if not (0<=ylo<=yhi<=1 and 0<=blo<=bhi<alo<=ahi<=1):
        raise ValueError('invalid or overlapping ordered state emission bounds')
    if yhi<blo or ylo>ahi:raise ValueError('measurement incompatible with either state')
    pmin=max(Fraction(0),(ylo-bhi)/(ahi-bhi))
    pmax=min(Fraction(1),(yhi-blo)/(alo-blo))
    if pmin>pmax:raise ValueError('no two-state mixture compatible with measurement')
    return pmin,pmax

def invert_three_state_emission(y_lo,y_hi, a_lo,a_hi,b_lo,b_hi,u_lo,u_hi,u_max=Fraction(1)):
    """SHARP rational bounds on fold-A fraction with 3 states A/B/U.

    Given 0<=pA,pB,pU, sum=1, 0<=pU<=u_max, and independent
    state-conditional signal bands, find exact feasible min/max pA.
    This is a two-variable LP (pA,pU) solved by rational vertex enumeration.
    It does not infer an unfolded bound: u_max must be externally justified.
    """
    ylo,yhi,alo,ahi,blo,bhi,ulo,uhi,um=map(_frac,(y_lo,y_hi,a_lo,a_hi,b_lo,b_hi,u_lo,u_hi,u_max))
    if not (0<=ylo<=yhi<=1 and 0<=alo<=ahi<=1 and 0<=blo<=bhi<=1
            and 0<=ulo<=uhi<=1 and 0<=um<=1):
        raise ValueError('invalid observations, emissions, or unfolded-fraction bound')
    # a*x+b*u <= c; pB=1-pA-pU
    planes=[(-1,0,0),(0,-1,0),(1,1,1),(0,1,um),
            (alo-blo,ulo-blo,yhi-blo),
            (bhi-ahi,bhi-uhi,bhi-ylo)]
    vertices=set()
    for i in range(len(planes)):
        a,b,c=planes[i]
        for j in range(i+1,len(planes)):
            d,e,f=planes[j]
            det=a*e-b*d
            if det==0:continue
            x=(c*e-b*f)/det
            u=(a*f-c*d)/det
            if all(A*x+B*u<=C for A,B,C in planes):vertices.add((x,u))
    if not vertices:raise ValueError('no compatible three-state population under supplied constraints')
    return {'fold_A_lower':min(v[0] for v in vertices),
            'fold_A_upper':max(v[0] for v in vertices),
            'unfolded_lower':min(v[1] for v in vertices),
            'unfolded_upper':max(v[1] for v in vertices),
            'feasible_rational_vertices':len(vertices)}