"""Independent exhaustive portfolio membership and optimization oracles."""
from fractions import Fraction
from itertools import product
from dataclasses import replace
import pytest

from foldbounds.lattice import ROOT,_neighbors,CertificateError
from foldbounds.shell_portfolio import (PortfolioCertificate, optimize_shell_portfolio,
                                        verify_shell_portfolio)


def complete_walks(n):
    paths=[ROOT]
    for _ in range(2,n):
        paths=[w for p in paths for w in _neighbors(p)]
    return paths


def selected(full, leaf):
    if full[:len(leaf)]!=leaf:
        return False
    r=len(leaf)
    if r==len(full):
        return True
    x,y=leaf[-1]
    for sx,sy in product((-1,1),repeat=2):
        phi=lambda site:sx*site[0]+sy*site[1]
        if phi((x,y))==max(phi(t) for t in leaf) and all(
                phi(full[i])-phi(full[i-1])==1 for i in range(r,len(full))):
            return True
    return False


def contact(seq,p):
    return sum(seq[i]=='H' and seq[j]=='H' and
               abs(p[i][0]-p[j][0])+abs(p[i][1]-p[j][1])==1
               for i in range(len(seq)) for j in range(i+2,len(seq)))


def brute_opt(seq,p,q,t,maxdepth,budget,which):
    def portfolios(path):
        yield (path,)
        if len(path)==maxdepth:
            return
        children=_neighbors(path)
        if not children:
            return
        pools=[tuple(portfolios(c)) for c in children]
        for parts in product(*pools):
            allparts=sum(parts,())
            if len(allparts)<=budget:
                yield allparts
    candidates=[]
    for leaves in portfolios(p):
        if len(leaves)>budget:continue
        cert=PortfolioCertificate(seq,p,q,t,tuple(leaves))
        m=verify_shell_portfolio(cert)
        candidates.append((getattr(m,which),m.leaf_count))
    return max(v for v,_ in candidates)

@pytest.mark.parametrize('n',[3,4,5,6,7])
@pytest.mark.parametrize('seqtype',['allH','mixed','allP'])
@pytest.mark.parametrize('event',[0,1,2])
def test_portfolio_cert_exact_small_oracle(n,seqtype,event):
    seq=('H'*n if seqtype=='allH' else 'P'*n if seqtype=='allP' else ('HPH'*n)[:n]);q=Fraction(3,2)
    v=optimize_shell_portfolio(seq,ROOT,q,event,leaf_budget=12,max_core_length=min(n,5))
    mass=verify_shell_portfolio(v.certificate)
    assert mass==v.mass
    walks=complete_walks(n)
    selected_walks=[p for p in walks if any(selected(p,leaf) for leaf in v.certificate.leaves)]
    e=sum((q**contact(seq,p) for p in selected_walks if contact(seq,p)>=event),Fraction(0))
    c=sum((q**contact(seq,p) for p in selected_walks if contact(seq,p)<event),Fraction(0))
    assert (mass.event,mass.complement,mass.distinct_paths)==(e,c,len(selected_walks))

@pytest.mark.parametrize('seq',["HHHHHH","HPHHPH","PPPPPP"])
@pytest.mark.parametrize('objective',['event','complement'])
def test_tree_knapsack_optimality_against_all_portfolios(seq,objective):
    q=Fraction(2); t=1
    plan=optimize_shell_portfolio(seq,ROOT,q,t,leaf_budget=8,max_core_length=5,objective=objective)
    assert getattr(plan.mass,objective)==brute_opt(seq,ROOT,q,t,5,8,objective)

@pytest.mark.parametrize('budget',[1,2,3,6,12,20])
def test_budget_monotonicity(budget):
    seq='HPHHPHHPHHHH'; p=optimize_shell_portfolio(seq,ROOT,Fraction(3,2),2,
                           leaf_budget=budget,max_core_length=8)
    assert p.mass.leaf_count<=budget
    if budget>1:
        prev=optimize_shell_portfolio(seq,ROOT,Fraction(3,2),2,
                              leaf_budget=budget-1,max_core_length=8)
        assert p.mass.event>=prev.mass.event


def test_invalid_and_adversarial_cover_rejected():
    seq='HHHHHH';p=optimize_shell_portfolio(seq,ROOT,Fraction(2),1,leaf_budget=7,max_core_length=5)
    cert=p.certificate
    with pytest.raises(CertificateError):
        verify_shell_portfolio(replace(cert,leaves=cert.leaves+cert.leaves[:1]))
    with pytest.raises(CertificateError):
        verify_shell_portfolio(replace(cert,leaves=(ROOT,)+cert.leaves))
    with pytest.raises(CertificateError):
        verify_shell_portfolio(replace(cert,leaves=(ROOT+((1,1),),)))
    with pytest.raises(CertificateError):
        verify_shell_portfolio(replace(cert,leaves=(ROOT+((0,0),),)))
    with pytest.raises(ValueError):
        optimize_shell_portfolio(seq,ROOT,Fraction(2),1,leaf_budget=True,max_core_length=5)
    with pytest.raises(ValueError):
        optimize_shell_portfolio(seq,ROOT,Fraction(2),1,leaf_budget=5,max_core_length=True)
    with pytest.raises(ValueError):
        optimize_shell_portfolio(seq,ROOT,Fraction(2),True,leaf_budget=5,max_core_length=5)


def test_independent_complete_twelve_residue_oracle_for_25_leaf_optimum():
    """All 30,073 canonical n=12 SAWs, rescored without portfolio shortcuts."""
    seq='HPHHPHHPHHHH'; q=Fraction(3,2); t=2
    opt=optimize_shell_portfolio(seq,ROOT,q,t,
                                  leaf_budget=25,max_core_length=9,
                                  objective='event')
    walks=complete_walks(len(seq))
    assert len(walks)==30073
    chosen=[p for p in walks if any(selected(p,leaf)
                                    for leaf in opt.certificate.leaves)]
    event=sum((q**contact(seq,p) for p in chosen if contact(seq,p)>=t),Fraction(0))
    other=sum((q**contact(seq,p) for p in chosen if contact(seq,p)<t),Fraction(0))
    assert (len(chosen),event,other)==(4510,Fraction(4455,4),Fraction(4464))
    assert (opt.mass.distinct_paths,opt.mass.event,opt.mass.complement)==(
        len(chosen),event,other)