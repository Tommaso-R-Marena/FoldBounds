# External SOTA comparison

This document records **real external benchmarks** (not in-repo toys) and how
FoldBounds compares to published state of the art. Machine-readable inputs:
`experiments/sota/benchmarks.json`. Latest numbers:
`experiments/sota/results/SUMMARY.md`.

## Benchmark suites adopted

| Suite | Source | What is scored |
|---|---|---|
| OEIS A001411 | Jensen / OEIS | Exact n-step SAW counts on Z²; FoldBounds uses A001411(n−1)/4 (anchored first bond) |
| S1-1 … S1-8 | Shmygelska–Hoos 2005; Unger–Moult | Classic 2D HP ground-state energies E* |
| Fibonacci S13/S21/S34 | AIP Advances 2015 | Exact enumeration targets for HP thermodynamics |

## SOTA baselines we run

| Method | Role | Trust |
|---|---|---|
| **MINWALK** (Giaquinta–Pozzi 2013; vendored `third_party/minwalk`) | Exact ground-state SOTA on 2D square HP | Unverified C++ |
| **hpenum** | Naive exhaustive enumerator | Unverified Python |
| **foldkernel** | Certified FOLDBOUNDS-KERNEL 1 checker | Lean-proved soundness |
| REMC / PERM / ACO (literature) | Heuristic GS SOTA | Cited, not re-run every time |

## Competitive claims (honest)

1. **Match published ground-state optima** on the community S1/Unger–Moult suite
   using the exact SOTA solver (MINWALK), with FoldBounds certificates accepted
   on the same instances.
2. **Exact match to OEIS A001411** (anchored) for all tested n; foldkernel’s
   in-process enumeration is competitive with / faster than the Python baseline
   at moderate n.
3. **Exact contact histograms** on Fibonacci S13 agreeing between foldkernel and
   exhaustive oracle.
4. **Unique vs heuristic SOTA (REMC/PERM/ACO):** FoldBounds returns
   machine-checked rational enclosures of Z and P(k≥t), not sampled energies.

Transfer-matrix SAW enumeration (Jensen to length 79) remains ahead for pure
counting at extreme n; FoldBounds’ differentiator is the verified checker and
HP thermodynamics certificates on the same lattice model.
