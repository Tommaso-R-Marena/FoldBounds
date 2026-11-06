"""Reproduce exact kinetics-vs-equilibrium certificates; no random simulation."""
from __future__ import annotations

from fractions import Fraction
from time import perf_counter

from foldbounds.kinetic_barriers import propose_kinetic_basin, propose_optimal_kinetic_cut
from foldbounds.kinetic_barrier_verifier import verify_kinetic_basin


def main() -> None:
    q, horizon = Fraction(10000), 10**9
    print("n,states,native_equilibrium_prob,basin_size,cut_flux,barrier,"
          "basin_start_hit_upper,fixed_start_hit_upper,optimality_audited,seconds")
    for n in (9, 11):
        start = perf_counter()
        c = propose_optimal_kinetic_cut("H" * n, q)
        r = verify_kinetic_basin(c, q, horizon)
        assert r.optimal_cut_verified_at_q
        print(f"{n},{r.total_states},{float(r.native_equilibrium_probability):.12g},"
              f"{r.basin_states},{r.boundary_flux},{r.barrier_order},"
              f"{float(r.exit_by_steps_upper):.12g},"
              f"{float(r.deterministic_seed_exit_by_steps_upper):.12g},"
              f"yes,{perf_counter()-start:.3f}")
        if n == 11:
            baseline = verify_kinetic_basin(
                propose_kinetic_basin("H"*n, q, min_barrier=4), q, horizon)
            print(f"# unoptimized cut flux={baseline.boundary_flux}, "
                  f"mincut flux={r.boundary_flux}, "
                  f"improvement={baseline.boundary_flux/r.boundary_flux}x")
    start = perf_counter()
    frozen = verify_kinetic_basin(
        propose_kinetic_basin("H"*12, q, allow_closed=True), q, horizon)
    assert frozen.closed_class
    print(f"12,{frozen.total_states},{float(frozen.native_equilibrium_probability):.12g},"
          f"{frozen.basin_states},0,infinite,0,0,not-applicable,"
          f"{perf_counter()-start:.3f}")
    print("# probabilities are exact Fractions internally; CSV is rounded for display")
    print("# one step = one uniformly drawn bead/compass proposal (including rejection)")
    print("# all results concern a 2D HP toy model, not real protein kinetic time")


if __name__ == "__main__":
    main()
