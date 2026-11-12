"""Independent counterexample search for spec T1.1, T1.3, T1.4, T1.5 (exact rationals)."""
import itertools, random
from fractions import Fraction as F
random.seed(12345)

def g_minus(reg, lam):
    return sum(min((r-lam)*L, (r-lam)*U) for (L,U,r,R) in reg)
def g_plus(reg, lam):
    return sum(max((R-lam)*L, (R-lam)*U) for (L,U,r,R) in reg)

def brute_min(reg):
    # min over vertices alpha in {L,U}^K, beta=r*alpha (linear-fractional => vertex)
    best=None
    for choice in itertools.product([0,1], repeat=len(reg)):
        a=[(U if c else L) for c,(L,U,r,R) in zip(choice,reg)]
        s=sum(a)
        if s<=0: continue
        q=sum(r*x for x,(L,U,r,R) in zip(a,reg))/s
        best=q if best is None or q<best else best
    return best
def brute_max(reg):
    best=None
    for choice in itertools.product([0,1], repeat=len(reg)):
        a=[(U if c else L) for c,(L,U,r,R) in zip(choice,reg)]
        s=sum(a)
        if s<=0: continue
        q=sum(R*x for x,(L,U,r,R) in zip(a,reg))/s
        best=q if best is None or q>best else best
    return best
def threshold_min(reg):
    order=sorted(range(len(reg)), key=lambda k: reg[k][2])
    best=None
    for t in range(len(reg)+1):
        a=[None]*len(reg)
        for pos,k in enumerate(order):
            L,U,r,R=reg[k]; a[k]= U if pos<t else L
        s=sum(a)
        if s<=0: continue
        q=sum(reg[k][2]*a[k] for k in range(len(reg)))/s
        best=q if best is None or q<best else best
    return best

def rnd():
    v=random.choice([0,0,1,2,3,5,7,F(1,2),F(1,3),F(5,4),10])
    return F(v)
fails=0; N=0
for trial in range(40000):
    K=random.randint(1,6)
    reg=[]
    for _ in range(K):
        L=rnd(); U=L+rnd(); r=rnd(); R=r+rnd()
        reg.append((L,U,r,R))
    if sum(U for L,U,r,R in reg)==0: continue
    N+=1
    lm=brute_min(reg); lp=brute_max(reg); tm=threshold_min(reg)
    assert lm==tm, (reg,lm,tm)
    # T1.1 exact test at lm and slightly above
    assert g_minus(reg,lm)>=0, ("g-(lam-)<0",reg,lm)
    assert g_minus(reg,lm)==0, ("g-(lam-)!=0",reg,lm,g_minus(reg,lm))
    eps=F(1,10**6)
    assert g_minus(reg,lm+eps)<0, ("g- not <0 above",reg,lm)
    assert g_plus(reg,lp)<=0 and g_plus(reg,lp-eps)>0
    # T1.4 dominance
    SU=sum(U for L,U,r,R in reg); SL=sum(L for L,U,r,R in reg)
    Sm=sum(r*L for L,U,r,R in reg)/SU
    assert lm>=Sm
    assert lm>=min(r for L,U,r,R in reg if U>0)
    if SL>0:
        Sp=sum(R*U for L,U,r,R in reg)/SL
        assert lp<=Sp
        # T1.5 upper iff
        cond=all(R*(U-L)==0 for L,U,r,R in reg)
        assert (lp==Sp)==cond, ("T1.5 upper", reg, lp, Sp)
    # T1.5 lower iff (when Sm>0)
    if Sm>0:
        cond=all(r*(U-L)==0 for L,U,r,R in reg)
        assert (lm==Sm)==cond, ("T1.5 lower",reg,lm,Sm)
print("T1 checks passed on", N, "random exact instances")
