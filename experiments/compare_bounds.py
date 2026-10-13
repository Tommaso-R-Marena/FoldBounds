"""Reproducible comparison with the original reference bounds; not SOTA."""
from fractions import Fraction
from time import perf_counter
from foldbounds.lattice import build_frontier, exact_spectrum, verify_certificate
from foldbounds.lookahead import verify_lookahead


def main():
    for seq in ('HHHHHHHH', 'HPHHPHHP', 'HHPHHHPHHH'):
        q, threshold = Fraction(2), 1
        cert = build_frontier(seq, q, threshold, 2)
        density = exact_spectrum(seq)
        z = sum((count * q**score for score,count in density.items()), Fraction(0))
        ev = sum((count * q**score for score,count in density.items() if score>=threshold), Fraction(0))
        p = ev / z
        print(f'\nSequence: {seq}; q={q}; k>={threshold}; exact Z={z}; exact P={p}')
        for name, verifier in [('baseline', lambda:verify_certificate(cert))] + [
            (f'lookahead-{k}', lambda k=k:verify_lookahead(cert,lookahead=k)) for k in (0,2,4)
        ]:
            start = perf_counter()
            b = verifier()
            elapsed = perf_counter() - start
            assert b.partition_lower <= z <= b.partition_upper
            assert b.probability_lower <= p <= b.probability_upper
            print(f'{name:13} Z=[{b.partition_lower},{b.partition_upper}] '
                  f'P=[{float(b.probability_lower):.6f},{float(b.probability_upper):.6f}] '
                  f'seconds={elapsed:.5f}')


if __name__ == '__main__':
    main()
