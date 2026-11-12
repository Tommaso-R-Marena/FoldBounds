"""Prior-building experiment (coordinator): for HP H->P mutation at residue i, with prefix-frontier regions
at depth m and EXACT region masses, how wide is the T1 coupled ratio enclosure due to ratio uncertainty only?
Width_m = (sum R_T a_T - sum r_T a_T)/Z with r_T=q^-max c_i, R_T = q^-min c_i over completions in T."""
import math, random, sys
from fractions import Fraction as F
from collections import defaultdict
DIRS=[(1,0),(-1,0),(0,1),(0,-1)]
def enum(n):
    walks=[]
    path=[(0,0),(1,0)]; occ={(0,0):0,(1,0):1}
    def rec():
        if len(path)==n: walks.append(list(path)); return
        x,y=path[-1]
        for dx,dy in DIRS:
            p=(x+dx,y+dy)
            if p in occ: continue
            occ[p]=len(path); path.append(p); rec(); path.pop(); del occ[p]
    rec(); return walks
def contacts(seq,w):
    pos={p:k for k,p in enumerate(w)}; k=0; ci=defaultdict(int)
    for a,p in enumerate(w):
        if seq[a]!='H': continue
        for dx,dy in DIRS:
            b=pos.get((p[0]+dx,p[1]+dy))
            if b is not None and b>a+1 and seq[b]=='H': k+=1; ci[a]+=1; ci[b]+=1
    return k,ci
n=int(sys.argv[1]) if len(sys.argv)>1 else 12
q=F(2); random.seed(1)
walks=enum(n); print("n",n,"walks",len(walks))
for trial in range(4):
    seq=''.join(random.choice('HP') for _ in range(n))
    for i in [1, n//2, n-2]:
        if seq[i]!='H': continue
        data=[]
        for w in walks:
            k,ci=contacts(seq,w); data.append((w,k,ci.get(i,0)))
        Z=sum(q**k for _,k,_ in data); ZB=sum(q**(k-c) for _,k,c in data)
        Qx=ZB/Z
        out=[]
        for m in range(2,n+1,max(1,(n-2)//5)):
            reg=defaultdict(lambda:[F(0),99,-1])
            for w,k,c in data:
                key=tuple(w[:m]); r=reg[key]; r[0]+=q**k; r[1]=min(r[1],c); r[2]=max(r[2],c)
            lo=sum(a*q**(-cmax) for a,cmin,cmax in reg.values())/Z
            hi=sum(a*q**(-cmin) for a,cmin,cmax in reg.values())/Z
            assert lo<=Qx<=hi
            out.append((m,len(reg),float(math.log(hi/lo))))
        triv=float(math.log(max(1,q**2)))
        print(f"seq={seq} i={i} exact dF={-math.log(Qx):.4f} trivial width={triv:.3f} | (m,#regions,width with EXACT masses):",
              ' '.join(f"({m},{nr},{wd:.3f})" for m,nr,wd in out))
