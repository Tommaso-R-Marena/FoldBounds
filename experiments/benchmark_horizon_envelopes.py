"""Standalone reproducible model-relative horizon upper-envelope benchmark."""
from fractions import Fraction
from time import perf_counter
from foldbounds.horizon_certificate import (
    build_horizon_certificate,verify_horizon_certificate,
    select_best_horizon_certificate,_component_possible)


def main():
    for n in (10,12,20,30):
        q=Fraction(2); sequence='H'*n
        t=perf_counter()
        coarse=verify_horizon_certificate(build_horizon_certificate(sequence,horizon=0),q)
        selected=select_best_horizon_certificate(sequence,q=q,horizon=8)
        fine=verify_horizon_certificate(selected,q)
        print({'length':n,'q':str(q),'horizon_selected':selected.horizon,
               'partition_upper_base':str(coarse.partition_upper),
               'partition_upper_refined':str(fine.partition_upper),
               'improvement_factor':str(coarse.partition_upper/fine.partition_upper),
               'dead_leaves':fine.dead_leaves,'elapsed_s':round(perf_counter()-t,4)})
    # A geometrically trapped prefix with a free neighboring cell.
    trapped=((0,0),(1,0),(2,0),(2,1),(2,2),(1,2),(0,2),(-1,2),(-1,1),(0,1))
    assert not _component_possible(trapped,14)
    c=build_horizon_certificate('H'*14,trapped,horizon=0)
    b=verify_horizon_certificate(c,Fraction(2))
    assert b.partition_upper==0
    print({'nonimmediately_trapped_prefix':trapped,'complete_extensions_upper':str(b.partition_upper)})

if __name__=='__main__':main()
