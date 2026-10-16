"""Exact rational finite-memory transfer ablation against matched settings."""
import csv,sys
from fractions import Fraction
from time import perf_counter
from foldbounds.lattice import build_frontier,exact_spectrum
from foldbounds.lookahead import verify_lookahead

def run():
    fields=('sequence','n','q','threshold','memory','elapsed_ms',
            'z_lower','z_upper','uncertainty_ratio','e_lower','e_upper',
            'other_lower','other_upper','p_lower','p_upper','p_width',
            'exact_z_if_known','exact_p_if_known')
    w=csv.DictWriter(sys.stdout,fieldnames=fields);w.writeheader()
    cases=(('H'*8,Fraction(2),1),('H'*12,Fraction(2),1),
           ('H'*16,Fraction(2),1),('H'*18,Fraction(2),1),
           ('HPHHPHHPHHHH',Fraction(3,2),2))
    for seq,q,t in cases:
        cert=build_frontier(seq,q,t,2)
        if len(seq)<=9:
            hist=exact_spectrum(seq)
            z=sum((n*q**k for k,n in hist.items()),Fraction(0))
            e=sum((n*q**k for k,n in hist.items() if k>=t),Fraction(0))
            p=e/z
        else:z=p=None
        prev=None
        for memory in (2,4,6):
            start=perf_counter()
            b=verify_lookahead(cert,lookahead=0,use_matching=True,
                use_transfer=True,use_monotone=True,transfer_memory=memory)
            ms=1000*(perf_counter()-start)
            if prev is not None:
                assert b.partition_upper<=prev.partition_upper
                assert b.event_upper<=prev.event_upper
                assert b.complement_upper<=prev.complement_upper
                assert b.partition_lower==prev.partition_lower
                assert b.probability_lower>=prev.probability_lower
                assert b.probability_upper<=prev.probability_upper
            if z is not None:
                assert b.partition_lower<=z<=b.partition_upper
                assert b.probability_lower<=p<=b.probability_upper
            prev=b
            w.writerow(dict(sequence=seq,n=len(seq),q=q,threshold=t,
                memory=memory,elapsed_ms=f'{ms:.3f}',
                z_lower=b.partition_lower,z_upper=b.partition_upper,
                uncertainty_ratio=b.partition_upper/b.partition_lower,
                e_lower=b.event_lower,e_upper=b.event_upper,
                other_lower=b.complement_lower,other_upper=b.complement_upper,
                p_lower=b.probability_lower,p_upper=b.probability_upper,
                p_width=b.probability_upper-b.probability_lower,
                exact_z_if_known=z if z is not None else '',
                exact_p_if_known=p if p is not None else ''))

if __name__=='__main__':run()
