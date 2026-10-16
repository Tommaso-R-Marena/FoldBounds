"""Compare exact-rational monotone family lower bounds at fixed other settings.

No claim of algorithmic novelty, experimentally realistic proteins, or
a Lean-verified Python implementation. Larger-model rows have no exact oracle.
"""
import csv
import sys
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier,exact_spectrum
from foldbounds.lookahead import verify_lookahead


def main():
    cols=('sequence','n','q','event_at_least','frontier_depth','witness_nodes',
          'monotone','milliseconds','z_lower','z_upper','event_lower','other_lower',
          'p_lower','p_upper','p_width','exact_z_or_blank','exact_p_or_blank')
    w=csv.DictWriter(sys.stdout,fieldnames=cols)
    w.writeheader()
    cases=(
        ('H'*8,Fraction(2),1),
        ('H'*12,Fraction(2),1),
        ('H'*16,Fraction(2),1),
        ('H'*18,Fraction(2),1),
        ('HPHHPHHPHHHH',Fraction(3,2),2),
    )
    for seq,q,t in cases:
        cert=build_frontier(seq,q,t,depth=2)
        if len(seq)<=10:
            histogram=exact_spectrum(seq)
            z=sum((count*q**k for k,count in histogram.items()),Fraction(0))
            a=sum((count*q**k for k,count in histogram.items() if k>=t),Fraction(0))
            p=a/z
        else:
            z=p=None
        for nodes in (0,64):
            before=None
            for enabled in (False,True):
                start=perf_counter()
                bounds=verify_lookahead(cert,lookahead=0,use_matching=True,
                    use_transfer=True,witness_nodes=nodes,use_monotone=enabled)
                elapsed=(perf_counter()-start)*1000
                if before:
                    assert bounds.partition_lower>=before.partition_lower
                    assert bounds.event_lower>=before.event_lower
                    assert bounds.complement_lower>=before.complement_lower
                    assert bounds.partition_upper<=before.partition_upper
                    assert bounds.probability_lower>=before.probability_lower
                    assert bounds.probability_upper<=before.probability_upper
                if z is not None:
                    assert bounds.partition_lower<=z<=bounds.partition_upper
                    assert bounds.probability_lower<=p<=bounds.probability_upper
                before=bounds
                w.writerow(dict(sequence=seq,n=len(seq),q=q,event_at_least=t,
                    frontier_depth=2,witness_nodes=nodes,monotone=enabled,
                    milliseconds=f'{elapsed:.3f}',z_lower=bounds.partition_lower,
                    z_upper=bounds.partition_upper,event_lower=bounds.event_lower,
                    other_lower=bounds.complement_lower,
                    p_lower=bounds.probability_lower,
                    p_upper=bounds.probability_upper,
                    p_width=bounds.probability_upper-bounds.probability_lower,
                    exact_z_or_blank=z if z is not None else '',
                    exact_p_or_blank=p if p is not None else ''))


if __name__=='__main__':
    main()
