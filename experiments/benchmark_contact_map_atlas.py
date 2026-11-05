"""Reproducible repeated-sequence benchmark against production exact_spectrum.

No timing threshold is asserted: timing varies by machine. Correctness is
asserted for every sequence before the results are printed.
"""
from __future__ import annotations

from itertools import product
from time import perf_counter

from foldbounds.contact_map_atlas import build_atlas
from foldbounds.lattice import exact_spectrum


def main() -> None:
    print("n,sequences,walks,distinct_maps,atlas_seconds,repeated_exact_seconds,speedup")
    for n in (8, 9):
        sequences = ("".join(x) for x in product("HP", repeat=n))
        t0 = perf_counter()
        atlas = build_atlas(n)
        atlas_values = {s: atlas.spectrum(s) for s in sequences}
        t_atlas = perf_counter() - t0
        t0 = perf_counter()
        direct_values = {s: exact_spectrum(s) for s in atlas_values}
        t_naive = perf_counter() - t0
        assert atlas_values == direct_values, "atlas disagrees with reference enumerator"
        print(f"{n},{len(atlas_values)},{atlas.conformations},{len(atlas.multiplicities)},"
              f"{t_atlas:.6f},{t_naive:.6f},{t_naive / t_atlas:.2f}")


if __name__ == "__main__":
    main()
