"""Reproduce robust shape-rank upper-mass ablations (no full SAW enumeration)."""
from fractions import Fraction
from time import perf_counter
from foldbounds.shape_rank_upper import build_shape_rank_certificate, verify_shape_rank_certificate
from foldbounds.saw_count_cap import optimal_saw_count_cap
from foldbounds.horizon_certificate import ROOT,_score_interval

def main():
    for n in (12,16,20,24):
        seq='H'*n
        t=perf_counter()
        c=build_shape_rank_certificate(seq,ROOT,suffix_memory=4)
        for q in (Fraction(1),Fraction(3,2),Fraction(2),Fraction(3)):
            z=verify_shape_rank_certificate(c,q)
            baseline=optimal_saw_count_cap(n-2).count_upper*q**_score_interval(seq,ROOT)[1]
            print({'n':n,'q':str(q),'gain':str(baseline/z.partition_upper),
                   'score_classes':z.score_classes,'ranked_count':z.ranked_paths,
                   'partition_upper':str(z.partition_upper),
                   'elapsed_s':round(perf_counter()-t,3)},flush=True)

if __name__=='__main__':main()
