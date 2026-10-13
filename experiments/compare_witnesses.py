"""Reproducible exact-rational witness-lower-bound ablation.

All true values are obtained independently from exact enumeration.
Reports both CPU duration and exact event/complement lower masses.
No claim of novelty or scaling breakthrough.
"""
import csv
import sys
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier, exact_spectrum
from foldbounds.lookahead import verify_lookahead


def main():
    columns=('sequence','q','event_at_least','frontier_depth','lookahead',
             'witness_nodes','time_ms','z_exact','z_lower','z_upper',
             'event_lower','complement_lower','p_exact','p_lower','p_upper',
             'probability_width')
    writer=csv.DictWriter(sys.stdout,fieldnames=columns)
    writer.writeheader()
    for sequence,q,threshold in (
        ('HHHHHHHH',Fraction(2),1),
        ('HPHHPHHP',Fraction(3,2),1),
        ('HHPHHHPHHH',Fraction(2),2),
    ):
        spectrum=exact_spectrum(sequence)
        exact_z=sum((count*q**k for k,count in spectrum.items()),Fraction(0))
        exact_e=sum((count*q**k for k,count in spectrum.items()
                     if k>=threshold),Fraction(0))
        exact_p=exact_e/exact_z
        for depth in (2,4):
            for lookahead in (0,2):
                previous=None
                for budget in (0,16,64,256):
                    cert=build_frontier(sequence,q,threshold,depth)
                    start=perf_counter()
                    bound=verify_lookahead(cert,lookahead=lookahead,
                        use_matching=True,use_transfer=True,
                        witness_nodes=budget)
                    time_ms=1000*(perf_counter()-start)
                    assert bound.event_lower<=exact_e<=bound.event_upper
                    assert bound.complement_lower<=exact_z-exact_e<=bound.complement_upper
                    assert bound.partition_lower<=exact_z<=bound.partition_upper
                    assert bound.probability_lower<=exact_p<=bound.probability_upper
                    if previous is not None:
                        assert previous.event_lower<=bound.event_lower
                        assert previous.complement_lower<=bound.complement_lower
                        assert previous.partition_lower<=bound.partition_lower
                        assert bound.partition_upper<=previous.partition_upper
                    previous=bound
                    writer.writerow(dict(sequence=sequence,q=q,
                        event_at_least=threshold,frontier_depth=depth,
                        lookahead=lookahead,witness_nodes=budget,
                        time_ms=f'{time_ms:.3f}',z_exact=str(exact_z),
                        z_lower=str(bound.partition_lower),
                        z_upper=str(bound.partition_upper),
                        event_lower=str(bound.event_lower),
                        complement_lower=str(bound.complement_lower),
                        p_exact=str(exact_p),
                        p_lower=str(bound.probability_lower),
                        p_upper=str(bound.probability_upper),
                        probability_width=str(bound.probability_upper-bound.probability_lower)))


if __name__ == '__main__':
    main()
