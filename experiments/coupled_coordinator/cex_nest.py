"""Brute-force check of T2.1 (enclosure) and T2.4 (nesting) on random pairwise MRFs on small graphs.
Perturbation: change unary factor at site 0 (A vs B). Separator radius R: interior = ball(R-1)... we use
interior I_R = {v: d(0,v) < R}, boundary S_R = {v : d(0,v) == R}, exterior = rest; valid since pairwise
factors only connect vertices at distance differing by <=1."""
import itertools, random
from fractions import Fraction as F
random.seed(3)
def bfs(adj,s):
    d={s:0}; fr=[s]
    while fr:
        nf=[]
        for u in fr:
            for w in adj[u]:
                if w not in d: d[w]=d[u]+1; nf.append(w)
        fr=nf
    return d
bad=0; checked=0
for trial in range(300):
    n=random.randint(4,8); k=random.choice([2,3])
    edges=set()
    # random connected graph
    for v in range(1,n): edges.add((random.randrange(v),v))
    for _ in range(random.randint(0,n)):
        a,b=random.sample(range(n),2); edges.add((min(a,b),max(a,b)))
    adj={v:set() for v in range(n)}
    for a,b in edges: adj[a].add(b); adj[b].add(a)
    pair={e:[[F(random.randint(1,9),random.randint(1,4)) for _ in range(k)] for _ in range(k)] for e in edges}
    unA=[[F(random.randint(1,9),3) for _ in range(k)] for v in range(n)]
    unB=[row[:] for row in unA]; unB[0]=[F(random.randint(1,9),2) for _ in range(k)]
    def w(conf,un):
        x=F(1)
        for v in range(n): x*=un[v][conf[v]]
        for (a,b),t in pair.items(): x*=t[conf[a]][conf[b]]
        return x
    confs=list(itertools.product(range(k),repeat=n))
    ZA=sum(w(c,unA) for c in confs); ZB=sum(w(c,unB) for c in confs); Q=ZB/ZA
    d=bfs(adj,0); maxd=max(d.values())
    prev=None
    for R in range(1,maxd+1):
        I=[v for v in range(n) if d[v]<R]; S=[v for v in range(n) if d[v]==R]
        # local weight: unaries of I, pair factors with an endpoint in I (others are exterior/boundary-only -> shared)
        def loc(conf,un):
            x=F(1)
            for v in I: x*=un[v][conf[v]]
            for (a,b),t in pair.items():
                if a in I or b in I: x*=t[conf[a]][conf[b]]
            return x
        rhos=[]
        for tau in itertools.product(range(k),repeat=len(S)):
            PA=PB=F(0)
            for sig in itertools.product(range(k),repeat=len(I)):
                conf=[0]*n
                for v,val in zip(I,sig): conf[v]=val
                for v,val in zip(S,tau): conf[v]=val
                PA+=loc(conf,unA); PB+=loc(conf,unB)
            rhos.append(PB/PA)
        lo,hi=min(rhos),max(rhos); checked+=1
        if not (lo<=Q<=hi): bad+=1; print("ENCLOSURE FAIL",trial,R,lo,Q,hi)
        if prev and not (prev[0]<=lo and hi<=prev[1]): bad+=1; print("NEST FAIL",trial,R,prev,(lo,hi))
        prev=(lo,hi)
    # full cover: R=maxd+1 means S empty -> exact
print("separator checks:",checked,"failures:",bad)
