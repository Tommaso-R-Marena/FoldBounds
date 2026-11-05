"""Benchmark all HP sequence partition functions against direct atlas queries.

Starts from shared independently audited geometry; compares the cost of
one-time polynomial expansion PLUS all queries with all direct contact-map
queries. Geometry construction is shared and excluded from both sides.
Results are machine-dependent and correctness is asserted before printing.

Examples:
  python experiments/universal_landscape_benchmark.py
  python experiments/universal_landscape_benchmark.py --lengths 14
"""
import argparse
from fractions import Fraction
from time import perf_counter

from foldbounds.contact_map_atlas import audit_atlas, build_atlas
from foldbounds.universal_landscape import build_landscape


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--lengths', type=int, nargs='+', default=[10, 12])
    parser.add_argument('--q', type=Fraction, default=Fraction(2))
    args = parser.parse_args()
    if args.q < 1 or any(not 2 <= n <= 14 for n in args.lengths):
        parser.error('q>=1 and lengths within [2,14] required')
    print('n,walks,maps,sequences,degree,transform_plus_query_s,direct_query_s,speedup')
    for n in args.lengths:
        atlas = build_atlas(n)
        if not audit_atlas(atlas, max_length=14):
            raise RuntimeError('independent geometry audit failed')
        sequences = tuple(''.join('H' if x & (1 << i) else 'P' for i in range(n))
                          for x in range(1 << n))
        t0 = perf_counter()
        landscape = build_landscape(atlas, require_audit=False)
        optimized = tuple(landscape.partition(s, args.q) for s in sequences)
        optimized_seconds = perf_counter() - t0
        t0 = perf_counter()
        baseline = tuple(atlas.partition(s, args.q) for s in sequences)
        direct_seconds = perf_counter() - t0
        assert optimized == baseline, 'universal transform disagreed with atlas oracle'
        print(f'{n},{atlas.conformations},{len(atlas.multiplicities)},{len(sequences)},'
              f'{len(landscape.polynomial("H" * n))-1},{optimized_seconds:.6f},'
              f'{direct_seconds:.6f},{direct_seconds / optimized_seconds:.2f}', flush=True)


if __name__ == '__main__':
    main()
