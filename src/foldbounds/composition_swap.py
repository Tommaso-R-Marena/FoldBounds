"""Exact, adversarially robust isocompositional cubic-lattice HP switches.

Event: >=2 geometric contacts (8-mer), >=3 (10-mer). Not atomistic
fold switching. Uses exact 3D enumeration and exact rational Bernstein
coefficients from parent research PRs; Python-to-Lean refinement is open.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from .mutation_switch import compact_spectrum, _bernstein_one
from .mutation_switch_10 import exact_spectrum10
from .ensemble_uncertainty import lower_probability, upper_probability

DEFAULT=(F(8),F(12),F(3,4),F(5,4),F(5,6),F(6,5))

def mask(seq):
    if type(seq) is not str or len(seq) not in (8,10) or set(seq)-{'H','P'}:
        raise ValueError('an 8- or 10-residue H/P sequence is required')
    return sum(1<<i for i,x in enumerate(seq) if x=='H')

def sequence(bits,n):
    if type(bits) is not int or type(n) is not int or n not in (8,10) or not 0<=bits<1<<n:
        raise ValueError('invalid sequence mask')
    return ''.join('H' if bits>>i&1 else 'P' for i in range(n))

def swapped_sites(a,b):
    if len(a)!=len(b):raise ValueError('sequence lengths differ')
    x,y=mask(a),mask(b)
    diff=x^y
    if diff.bit_count()!=2 or x.bit_count()!=y.bit_count():
        raise ValueError('exactly one H/P exchange is required')
    i,j=[k for k in range(len(a)) if diff>>k&1]
    if i%2!=j%2:
        raise ValueError('exchange must preserve H counts on both chain parity sublattices')
    return i+1,j+1

@lru_cache(maxsize=2048)
def nominal_box(bits,n,box):
    if type(box) is not tuple or len(box)!=6 or not all(type(x) is F for x in box):
        raise ValueError('six exact rational endpoints required')
    qa,qb,ba,bb,fl,fu=box
    if not (0<qa<=qb and 0<ba<=bb and 0<fl<=fu):
        raise ValueError('positive, ordered parameter intervals required')
    seq=sequence(bits,n)
    spectrum=compact_spectrum(bits) if n==8 else exact_spectrum10(seq)
    a,b,k=_bernstein_one(spectrum,qa,qb,ba,bb)
    return a,b,k

def corrected_box(bits,n,box):
    lower,upper,_=nominal_box(bits,n,box)
    return lower_probability(lower,box[4],box[5]),upper_probability(upper,box[4],box[5])

@dataclass(frozen=True)
class SwapCertificate:
    parent:str
    mutant:str
    sites:tuple[int,int]
    parameters:tuple[F,...]
    parent_robust_upper:F
    mutant_robust_lower:F
    separation:F
    conformations:int

def certify_swap(parent,mutant,parameters=DEFAULT):
    sites=swapped_sites(parent,mutant)
    n=len(parent)
    pl,pu=corrected_box(mask(parent),n,parameters)
    ml,mu=corrected_box(mask(mutant),n,parameters)
    if not pu<F(1,2)<ml:
        raise ValueError('strict robust minority/majority switch not certified')
    return SwapCertificate(parent,mutant,sites,parameters,pu,ml,ml-pu,
                           13565 if n==8 else 308981)

def verify_swap(c):
    if not isinstance(c,SwapCertificate):return False
    try:return c==certify_swap(c.parent,c.mutant,c.parameters)
    except (TypeError,ValueError,ArithmeticError,AssertionError):return False

@dataclass(frozen=True)
class Landscape:
    residues:int
    parameters:tuple[F,...]
    sequences:int
    allowed_edges:int
    certified_switches:tuple[tuple[int,int,int,int],...]
    reversal_orbits:int

def certify_landscape(n,parameters=DEFAULT):
    if type(n) is not int or n not in (8,10):raise ValueError('length must be 8 or 10')
    prob=[corrected_box(bits,n,parameters) for bits in range(1<<n)]
    switches=[];edges=0
    for m in range(1<<n):
        for i in range(n):
            if not m>>i&1:continue
            for j in range(n):
                if i%2!=j%2 or m>>j&1:continue
                t=m^(1<<i)^(1<<j)
                if m>=t:continue  # exactly one inspection per undirected edge
                edges+=1
                if prob[m][1]<F(1,2)<prob[t][0]:switches.append((m,t,i+1,j+1))
                elif prob[t][1]<F(1,2)<prob[m][0]:switches.append((t,m,j+1,i+1))
    expected=2*((n//2)*(n//2-1)//2)*(1<<(n-2))
    if edges!=expected:raise AssertionError('incomplete mutation graph')
    def rev(m):return sum(((m>>i)&1)<<(n-1-i) for i in range(n))
    orbits={min((a,b),(rev(a),rev(b))) for a,b,_,_ in switches}
    return Landscape(n,parameters,1<<n,edges,tuple(sorted(switches)),len(orbits))

def verify_landscape(c):
    if not isinstance(c,Landscape):return False
    try:return c==certify_landscape(c.residues,c.parameters)
    except (TypeError,ValueError,ArithmeticError,AssertionError):return False

@dataclass(frozen=True)
class ErrorBudgetCertificate:
    parent:str
    mutant:str
    q_lower:F
    q_upper:F
    bend_lower:F
    bend_upper:F
    max_symmetric_multiplier_squared:F
    limiting_condition:str

def certify_error_budget(parent,mutant,q_lower=F(8),q_upper=F(12),
                         bend_lower=F(3,4),bend_upper=F(5,4)):
    swapped_sites(parent,mutant)
    n=len(parent)
    box=q_lower,q_upper,bend_lower,bend_upper,F(1),F(1)
    a,b,_=nominal_box(mask(parent),n,box)
    c,d,_=nominal_box(mask(mutant),n,box)
    if not (0<b<1 and 0<c<1):raise ValueError('degenerate nominal bounds')
    parent_limit=(1-b)/b
    mutant_limit=c/(1-c)
    threshold=min(parent_limit,mutant_limit)
    if threshold<=1:raise ValueError('strict switch is not provable')
    return ErrorBudgetCertificate(parent,mutant,q_lower,q_upper,bend_lower,
               bend_upper,threshold,'parent' if parent_limit<=mutant_limit else 'mutant')

def verify_error_budget(c):
    if not isinstance(c,ErrorBudgetCertificate):return False
    try:return c==certify_error_budget(c.parent,c.mutant,c.q_lower,c.q_upper,
                                        c.bend_lower,c.bend_upper)
    except (TypeError,ValueError,ArithmeticError,AssertionError):return False
