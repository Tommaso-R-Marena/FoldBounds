"""Independent, exponentially costly 3D cubic lattice HP continuum certificates.

Research-only small-chain oracle. Not an all-atom protein force field.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from .ensemble_crossover import (
    probability, robust_bend_sign_certificate, unique_switch_bracket,
)

STEP_3D=((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))

@lru_cache(maxsize=None)
def geometry_3d(n:int):
    if type(n) is not int or not 2<=n<=10:
        raise ValueError('3D exact-oracle enumeration limited to n=2..10')
    path=[(0,0,0),(1,0,0)];occupied=set(path);out=[]
    def dfs(previous_direction,turns):
        if len(path)==n:
            pairs=tuple((i,j) for j in range(2,n) for i in range(j-1)
                        if sum(abs(path[i][k]-path[j][k]) for k in range(3))==1)
            out.append((pairs,turns));return
        x,y,z=path[-1]
        for direction in STEP_3D:
            p=(x+direction[0],y+direction[1],z+direction[2])
            if p not in occupied:
                path.append(p);occupied.add(p)
                dfs(direction,turns+(previous_direction!=direction))
                path.pop();occupied.remove(p)
    dfs((1,0,0),0)
    return tuple(out)

def spectrum_3d(sequence:str):
    if type(sequence) is not str or len(sequence)<2 or set(sequence)-{'H','P'}:
        raise ValueError('sequence must be H/P string length>=2')
    counter=Counter()
    for pair_list,turns in geometry_3d(len(sequence)):
        contacts=sum(sequence[i]=='H' and sequence[j]=='H' for i,j in pair_list)
        counter[(contacts,turns)]+=1
    return counter

def certified_3d_example():
    sequence='HPHHPHHH';threshold=1
    hist=spectrum_3d(sequence)
    count,negative,positive=robust_bend_sign_certificate(
        hist,threshold,Fraction(3,2),Fraction(1,4))
    if negative or not positive:raise AssertionError('3D positivity witness unavailable')
    qlo,qhi=Fraction(3,2),Fraction(2)
    blo,bhi=Fraction(1,4),Fraction(4)
    crossings={}
    for q in (qlo,qhi):
        lo,hi=unique_switch_bracket(hist,q,threshold,blo,bhi,bits=40)
        crossings[str(q)]=(lo,hi)
    return dict(sequence=sequence,model='3D anchored cubic-lattice HP',
                conformations=sum(hist.values()),histogram_bins=len(hist),
                coefficients_total=count,coefficients_positive=positive,
                interval=(probability(hist,qlo,blo,threshold),
                          probability(hist,qhi,bhi,threshold)),
                crossover_brackets=crossings)

if __name__=='__main__':
    print(certified_3d_example())