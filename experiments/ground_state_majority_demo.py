"""Reproduce exact all-H ground-state coefficients and strong-coupling odds."""
from __future__ import annotations
import argparse
from fractions import Fraction
from time import perf_counter

from foldbounds.ground_state_transfer import (
    build_ground_state_certificate, verify_ground_state_certificate,
    certify_ground_population, directed_ground_union,
)


def main():
    parser=argparse.ArgumentParser(description="Exact 2D all-H ground-state enumeration")
    parser.add_argument("--full",action="store_true",help="include exact 8x8 run")
    args=parser.parse_args()
    for side in ((2,3,4,5,6,7,8) if args.full else (2,3,4,5,6)):
        start=perf_counter()
        c=build_ground_state_certificate(side)
        assert verify_ground_state_certificate(c)
        print(f"{side}x{side}: n={side*side}, exact max score={c.max_contacts}, "
              f"ground states={c.degeneracy:,}, time={perf_counter()-start:.2f}s")
        if side==8:
            q=Fraction(10)
            p=certify_ground_population(c,q,30)
            robust=certify_ground_population(c,q,30,
                correction_lower=Fraction(1,2),correction_upper=Fraction(2))
            combined=directed_ground_union(c,q)
            print("count upper",p.walk_count_cap)
            print("P(K>=30) >= ",float(p.probability_lower))
            print("robust [1/2,2] P(K>=30) >= ",float(robust.probability_lower))
            print("exact directed overlap:",combined.directed_overlap)
            print("union / directed:",float(combined.union_lower/combined.directed_lower))
            old=combined.directed_lower+105840*q**49
            print("union / old compact+directed:",float(combined.union_lower/old))
            assert c.degeneracy==1511656642
            assert combined.directed_overlap==2
            assert p.probability_lower>Fraction(96,100)
            assert robust.probability_lower>Fraction(87,100)


if __name__=="__main__":
    main()
