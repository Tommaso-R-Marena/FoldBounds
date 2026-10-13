"""Compare anytime refinement policies on a matched 2D HP toy model.

All policies get the same number of prefix expansions, but their time
costs may differ greatly. Exact enumeration is used only as an oracle
after selecting each fixed test case, never by the search policies.
These results are NOT comparisons against published state of the art.
"""
from csv import DictWriter
from fractions import Fraction
from time import perf_counter

from foldbounds.adaptive import adaptive_refine, verify_refinement_run
from foldbounds.lattice import build_frontier, exact_spectrum


def main() -> None:
    cases = [
        ('HHHHHHHH', Fraction(2), 1),
        ('HPHHPHHP', Fraction(3,2), 1),
        ('HHPHHHPHHH', Fraction(2), 2),
    ]
    headers=('sequence','q','event_at_least','policy','step_budget','steps_executed',
             'partition_width','probability_width','elapsed_ms','exact_partition',
             'exact_probability','checked')
    writer=DictWriter(__import__('sys').stdout,fieldnames=headers)
    writer.writeheader()
    for sequence,q,threshold in cases:
        spectrum=exact_spectrum(sequence)
        z=sum((count*q**k for k,count in spectrum.items()),Fraction(0))
        event=sum((count*q**k for k,count in spectrum.items() if k>=threshold),Fraction(0))
        exact_p=event/z
        for budget in (4,8,16,32):
            for policy in ('breadth','mass','probability'):
                cert=build_frontier(sequence,q,threshold,2)
                start=perf_counter()
                result=adaptive_refine(cert,steps=budget,policy=policy,
                                       lookahead=0,use_matching=True)
                checked=verify_refinement_run(result)
                duration=(perf_counter()-start)*1000
                assert checked.partition_lower<=z<=checked.partition_upper
                assert checked.probability_lower<=exact_p<=checked.probability_upper
                writer.writerow(dict(sequence=sequence,q=q,event_at_least=threshold,
                    policy=policy,step_budget=budget,
                    steps_executed=len(result.history)-1,
                    partition_width=str(checked.partition_upper-checked.partition_lower),
                    probability_width=str(checked.probability_upper-checked.probability_lower),
                    elapsed_ms=f'{duration:.3f}',exact_partition=str(z),
                    exact_probability=str(exact_p),checked=True))


if __name__=='__main__':
    main()
