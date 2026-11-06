"""Reproducible anchored HP subset/thermodynamic benchmark (not protein data)."""
from __future__ import annotations

from fractions import Fraction as Q
from time import perf_counter

from foldbounds.lattice import build_frontier
from foldbounds.monotone_spectral import monotone_spectral_histogram
from foldbounds.partial_directed import partially_directed_spectrum
from foldbounds.thermal_response import thermal_response_bounds


def main() -> None:
    print("Model: 2D anchored HP square lattice, all-H, q=1")
    print("n,monotone_paths,directed_paths,path_gain,known_Z_gain,"
          "variance_lower_old,variance_lower_new,mean_upper_old,mean_upper_new,"
          "P_contacts>=1_lower_new,seconds_directed")
    root=((0,0),(1,0))
    for n in (12,20,30,40):
        seq='H'*n
        original_spectrum=monotone_spectral_histogram(seq,root)
        t=perf_counter()
        directed_spectrum=partially_directed_spectrum(seq,root)
        elapsed=perf_counter()-t
        original_paths=sum(c for _,c in original_spectrum)
        directed_paths=sum(c for _,c in directed_spectrum)
        cert=build_frontier(seq,Q(1),1,2)
        previous=thermal_response_bounds(cert,use_monotone_spectra=True,
                                         use_saw_count_cap=True)
        advanced=thermal_response_bounds(cert,use_partially_directed_spectra=True,
                                         use_saw_count_cap=True,
                                         use_geometric_contact_cap=True)
        assert advanced.partially_directed_paths==directed_paths
        assert advanced.variance_lower>previous.variance_lower
        assert advanced.partition_upper<=previous.partition_upper
        assert advanced.variance_upper<=previous.variance_upper
        print(','.join(map(str,(n,original_paths,directed_paths,
                                 f'{directed_paths/original_paths:.3f}',
                                 f'{float(advanced.partition_lower/previous.partition_lower):.3f}',
                                 str(previous.variance_lower),str(advanced.variance_lower),
                                 str(previous.mean_upper),str(advanced.mean_upper),
                                 str(advanced.probability_lower),
                                 f'{elapsed:.3f}'))))


if __name__=="__main__":
    main()
