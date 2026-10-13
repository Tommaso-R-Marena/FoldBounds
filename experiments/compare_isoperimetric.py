"""Compare exactly the same certificate verifier with/without perimeter caps.

Report wall time and the exact rational bounds; sample input is fixed.
These are comparisons against the immediately preceding FoldBounds
implementation, NOT against published methods or physical experiments.
"""
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead


def main():
    for seq,q,threshold in [
        ('HHHHHHHH',Fraction(2),1),
        ('HHPHHHPHHH',Fraction(2),2),
        ('HHHHHHHHHHHH',Fraction(3,2),2),
    ]:
        spectrum=exact_spectrum(seq)
        exact_z=sum((count*q**k for k,count in spectrum.items()),Fraction(0))
        exact_e=sum((count*q**k for k,count in spectrum.items()
                     if k>=threshold),Fraction(0))
        print(f'CASE seq={seq} q={q} event>={threshold} exact_Z={exact_z} exact_P={exact_e/exact_z}')
        for lookahead in (0,1,2):
            cert=build_frontier(seq,q,threshold,2)
            bounds=[]
            for perimeter in (False,True):
                started=perf_counter()
                res=verify_lookahead(cert,lookahead=lookahead,use_matching=True,
                                     use_perimeter=perimeter)
                elapsed_ms=1000*(perf_counter()-started)
                assert res.partition_lower<=exact_z<=res.partition_upper
                assert res.probability_lower<=exact_e/exact_z<=res.probability_upper
                bounds.append(res)
                print(f'  lookahead={lookahead} perimeter={perimeter} '
                      f'Z=[{res.partition_lower},{res.partition_upper}] '
                      f'P=[{res.probability_lower},{res.probability_upper}] '
                      f'elapsed_ms={elapsed_ms:.3f}')
            before,after=bounds
            assert before.partition_lower<=after.partition_lower<=exact_z
            assert exact_z<=after.partition_upper<=before.partition_upper
            assert before.probability_lower<=after.probability_lower
            assert after.probability_upper<=before.probability_upper


if __name__=='__main__':
    main()
