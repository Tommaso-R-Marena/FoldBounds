"""Bound width and CPU scaling without a full-conformation oracle.

Only candidate-model diagnostics; absence of an exact oracle means we
must not present these rows as independent validation of soundness.
"""
import csv
import sys
from fractions import Fraction
from time import perf_counter

from foldbounds.lattice import build_frontier
from foldbounds.lookahead import verify_lookahead


def main():
    fields=('length','sequence','q','event_at_least','witness_nodes',
            'z_lower','z_upper','upper_lower_ratio','p_lower','p_upper',
            'p_width','elapsed_ms')
    writer=csv.DictWriter(sys.stdout,fieldnames=fields)
    writer.writeheader()
    for n in (10,12,14,16,18):
        sequence='H'*n
        q=Fraction(2)
        threshold=1
        cert=build_frontier(sequence,q,threshold,depth=2)
        for budget in (0,64,256):
            start=perf_counter()
            b=verify_lookahead(cert,lookahead=0,use_matching=True,
                               use_transfer=True,witness_nodes=budget)
            ms=1000*(perf_counter()-start)
            assert 0<=b.partition_lower<=b.partition_upper
            assert 0<=b.probability_lower<=b.probability_upper<=1
            ratio=b.partition_upper/b.partition_lower if b.partition_lower else None
            writer.writerow(dict(length=n,sequence=sequence,q=q,
                event_at_least=threshold,witness_nodes=budget,
                z_lower=str(b.partition_lower),z_upper=str(b.partition_upper),
                upper_lower_ratio=(str(ratio) if ratio is not None else 'undefined'),
                p_lower=str(b.probability_lower),
                p_upper=str(b.probability_upper),
                p_width=str(b.probability_upper-b.probability_lower),
                elapsed_ms=f'{ms:.3f}'))


if __name__ == '__main__':
    main()
