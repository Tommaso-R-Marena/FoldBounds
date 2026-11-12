"""Check the soft-separator claim:
exterior couples to each boundary var v only via factors touching v alone (plus exterior vars),
with per-exterior-state log-range <= s_v  ==>  Q lies in [min,max] over vertex fields h_v(a) in {0,s_v},
and that the vertex extremes dominate random product fields and random mixture exteriors."""
import itertools, random, math
random.seed(7)

def Z(Phi, field):  # Phi: dict tau->weight; field: per v, list of multipliers per value
    tot=0.0
    for tau,w in Phi.items():
        m=1.0
        for v,a in enumerate(tau): m*=field[v][a]
        tot+=w*m
    return tot

viol=0; trials=0
for trial in range(3000):
    nb=random.randint(1,3); dom=[random.randint(2,3) for _ in range(nb)]
    states=list(itertools.product(*[range(d) for d in dom]))
    PhiA={t:random.uniform(0.1,3) for t in states}
    PhiB={t:random.uniform(0.1,3) for t in states}
    s=[random.uniform(0,2) for _ in range(nb)]
    # vertex fields: each value multiplier in {1, e^{s_v}}
    vals=[]
    choices=[list(itertools.product([1.0,math.exp(s[v])],repeat=dom[v])) for v in range(nb)]
    for fld in itertools.product(*choices):
        vals.append(Z(PhiB,fld)/Z(PhiA,fld))
    lo,hi=min(vals),max(vals)
    # random mixture exteriors: several exterior states, each with product factors with log-range <= s_v
    for _ in range(30):
        K=random.randint(1,5)
        numer=denom=0.0
        for k in range(K):
            g=random.uniform(0.01,5)
            fld=[[math.exp(random.uniform(0,s[v])) for a in range(dom[v])] for v in range(nb)]
            numer+=g*Z(PhiB,fld); denom+=g*Z(PhiA,fld)
        q=numer/denom; trials+=1
        if q<lo-1e-12 or q>hi+1e-12: viol+=1; print("VIOLATION",q,lo,hi)
print("soft-separator: trials",trials,"violations",viol)
