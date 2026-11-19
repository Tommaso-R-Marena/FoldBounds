# External SOTA comparison

This document records **real external benchmarks** (not in-repo toys) and how
FoldBounds compares to published state of the art. Machine-readable inputs:
`experiments/sota/benchmarks.json`. Latest numbers:
`experiments/sota/results/SUMMARY.md`.
Internal matched-budget kernel vs foldcheck2: `experiments/kernel/bench/results/SUMMARY.md`.

## Symmetry / anchoring (non-negotiable)

FoldBounds uses anchored first bond `(0,0)→(1,0)` and does **not** identify
reflections. Conformation count for length-`n` walks =
`A001411(n−1)/4`. Papers that report 8-fold reduced counts (e.g. Lau–Dill
`n=10 → 2034`) must be converted before claiming agreement
(`16268/8 = 2034`; FoldBounds-anchored = `4067`).

## Benchmark suites adopted

| Suite | Source | What is scored |
|---|---|---|
| OEIS A001411 | Jensen / OEIS | Exact n-step SAW counts on Z² (anchored) |
| S1-1 … S1-8 + 20merB | Shmygelska–Hoos; Unger–Moult; Istrail/Tortilla | Classic 2D HP ground-state energies E* |
| Fibonacci S13/S21/S34 | AIP Advances 2015 | Exact enumeration targets for HP thermodynamics |
| Literature REMC | Thachuk–Shmygelska–Hoos 2007 | Heuristic success % / mean runtime (cited) |

## SOTA baselines we run / cite

| Method | Role | Trust |
|---|---|---|
| **MINWALK** (Giaquinta–Pozzi 2013; vendored `third_party/minwalk`) | Exact ground-state SOTA; `PATH` dump for FoldBounds witnesses | Unverified C++ |
| **hpenum** | Naive exhaustive enumerator | Unverified Python |
| **foldkernel** | Certified FOLDBOUNDS-KERNEL 1 checker | Lean-proved soundness |
| REMC / PERM / ACO (literature) | Heuristic GS SOTA | Cited, not re-run every time |
| Jensen / Lee / Beaton TM (literature) | Exact SAW / ISAW Ω(K) oracles | Cited; next oracle integration |

## Competitive claims (honest)

1. **Match published ground-state optima** on the community S1/Unger–Moult /
   Istrail suite using MINWALK. Vendored MINWALK dumps an anchored `PATH` so
   FoldBounds can certify the same conformation.
2. **Certified eLo>0** on S1-1…S1-5, S1-8, and 20merB (event `k ≥ E*` nonempty)
   via MINWALK-path witness spines — REMC reports energies without certificates.
3. **Exact match to OEIS A001411** (anchored) for all tested n; foldkernel
   ~10× vs Python hpenum at n≈16.
4. **Exact contact histograms** on Fibonacci S13/S21; S21 ~7× vs exhaustive
   Python.
5. **Internal matched-budget bench:** 0 disagreements vs refs; foldkernel
   median ~60× vs foldcheck2 wall; S5 held-out wins 25/25 on certified `p_width`.

Transfer-matrix SAW enumeration (Jensen to length 79) remains ahead for pure
counting at extreme n. S1-6/S1-7 exact MINWALK is still open at short harness
budgets (long `--ci` runs in progress; paper itself calls S1-6 “huge”). REMC
stays faster when no certificate is required; FoldBounds’ differentiator is
machine-checked rational enclosures.

## Ranking: compete now vs aspirational

| Benchmark | Now | Later |
|---|---|---|
| OEIS SAW counts (anchored) | Yes | — |
| S1-1…S1-5, S1-8, 20merB certified eLo>0 | Yes | — |
| Fib S13/S21 exact histograms | Yes | S34 |
| Homopolymer ISAW Ω(K)/Z vs Lee–Beaton | Partial (all-H via OEIS total) | Exact tables N→42 |
| S1-6/S1-7 proven E* + cert | In progress | Done |
| S1-9…S1-11 ground energy | No | Aspirational |
| Full DOS for 2D100 | No | WL / flatPERM territory |

## Next adoption pack

1. Beaton `sq_isaw_counts.ser` / Lee Ω(K) tables with explicit symmetry adapter.
2. Finish S1-6/S1-7 MINWALK `--ci` + PATH certs.
3. Optional: Wüst Wang–Landau for approximate DOS cross-check (not ground truth).
