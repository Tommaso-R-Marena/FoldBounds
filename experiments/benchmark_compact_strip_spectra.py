"""Reproduce compact-strip exact-rational lower-mass improvements."""
from __future__ import annotations

from fractions import Fraction
from foldbounds.compact_strip_spectra import (
    build_compact_strip_certificate, verify_compact_strip_certificate,
    compact_plus_directed_mass,
)
from foldbounds.partial_directed import partially_directed_spectrum


def main():
    print('All-H: existing PDSAW lower vs verified disjoint compact+PDSAW lower')
    for n, height in ((32,8),(48,6),(64,8)):
        sequence='H'*n
        reference=dict(partially_directed_spectrum(sequence))
        for q in (2,10,100):
            prev=sum((count*Fraction(q)**k for k,count in reference.items()),Fraction(0))
            separate=verify_compact_strip_certificate(build_compact_strip_certificate(
                sequence,height=height,disjoint_from_partially_directed=True),Fraction(q))
            hybrid=compact_plus_directed_mass(sequence,Fraction(q),height=height)
            assert hybrid.partition_lower==prev+separate.partition_lower
            print(f'n={n:2} tile_h={height} q={q:3} ' 
                  f'compact={separate.witnessed_paths:,} '
                  f'combined/old={float(hybrid.partition_lower/prev):.5f}x')
    print('\nHeterogeneous sequence exact rational subset spectrum:')
    seq=('HPHHPHPH' * 5)[:40]
    cert=build_compact_strip_certificate(seq,height=5)
    bound=verify_compact_strip_certificate(cert,Fraction(7,3),1)
    print(f'n={len(seq)}, h=5, distinct paths={bound.witnessed_paths:,}, '
          f'contact classes={len(cert.spectrum)}, Z_sub={bound.partition_lower}')

if __name__=='__main__':
    main()
