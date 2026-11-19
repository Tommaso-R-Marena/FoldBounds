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
   using the exact SOTA solver (MINWALK). Vendored MINWALK now dumps an anchored
   `PATH` so FoldBounds can certify the same conformation.
2. **Certified eLo>0** on S1-1…S1-5 and S1-8 (event `k ≥ E*` nonempty) via
   MINWALK-path witness spines checked by the Lean kernel — the heuristic REMC
   SOTA (Thachuk 2007) reports energies without certificates.
3. **Exact match to OEIS A001411** (anchored) for all tested n; foldkernel
   ~10× vs Python hpenum at n≈16.
4. **Exact contact histograms** on Fibonacci S13/S21; S21 ~7× vs exhaustive
   Python.

Transfer-matrix SAW enumeration (Jensen to length 79) remains ahead for pure
counting at extreme n. S1-6/S1-7 exact MINWALK is still open at short harness
budgets (long runs in progress). REMC stays faster when no certificate is
required; FoldBounds’ differentiator is machine-checked rational enclosures.
