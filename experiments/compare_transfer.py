"""Exact-rational internal ablation of prefix-obstacle transfer upper bounds."""
import csv
import sys
from fractions import Fraction
from time import perf_counter
from foldbounds.lattice import build_frontier,exact_spectrum
from foldbounds.lookahead import verify_lookahead


def main():
    writer=csv.DictWriter(sys.stdout,fieldnames=[
        'sequence','q','event_at_least','depth','lookahead','method',
        'exact_z','exact_p','z_lower','z_upper','p_lower','p_upper',
        'p_width','milliseconds'])
    writer.writeheader()
    for seq,q,threshold in [
        ('HHHHHHHH',Fraction(2),1),
        ('HPHHPHHP',Fraction(3,2),1),
        ('HHPHHHPHHH',Fraction(2),2),
        ('HHHHHHHHHHHH',Fraction(2),1),
    ]:
        sp=exact_spectrum(seq)
        z=sum((n*q**k for k,n in sp.items()),Fraction(0))
        e=sum((n*q**k for k,n in sp.items() if k>=threshold),Fraction(0))
        p=e/z
        for depth in (2,4):
            cert=build_frontier(seq,q,threshold,depth)
            for lookahead in (0,2):
                previous=None
                for mode,transfer in (('matching+perimeter',False),
                                      ('matching+perimeter+transfer',True)):
                    start=perf_counter()
                    b=verify_lookahead(cert,lookahead=lookahead,use_matching=True,
                                       use_transfer=transfer)
                    ms=1000*(perf_counter()-start)
                    assert b.partition_lower<=z<=b.partition_upper
                    assert b.probability_lower<=p<=b.probability_upper
                    if previous is not None:
                        assert previous.partition_lower<=b.partition_lower
                        assert b.partition_upper<=previous.partition_upper
                        assert previous.probability_lower<=b.probability_lower
                        assert b.probability_upper<=previous.probability_upper
                    previous=b
                    writer.writerow(dict(sequence=seq,q=q,event_at_least=threshold,
                        depth=depth,lookahead=lookahead,method=mode,
                        exact_z=str(z),exact_p=str(p),z_lower=str(b.partition_lower),
                        z_upper=str(b.partition_upper),p_lower=str(b.probability_lower),
                        p_upper=str(b.probability_upper),
                        p_width=str(b.probability_upper-b.probability_lower),
                        milliseconds=f'{ms:.3f}'))


if __name__=='__main__':
    main()
