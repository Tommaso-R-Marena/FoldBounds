# Certified frontier session status (resumable checkpoint)

Last updated: 2026-11-18 20:20 America/New_York (kernel grid + soundness + benches).

## Git state

* Work branch: `feature/certified-frontier-kernel-grid-9cd6`
  (branched from `feature/fable52-certified-frontier-v1-9cd6`, itself based on
  draft PR #32 head `fc3de05`).
* Deliverable: Mathlib-free `FOLDBOUNDS-KERNEL 1` checker with kernel-checked
  string-level soundness, grid enumerator proved equal to the list enumerator,
  Python reference/oracle/adversarial corpus, matched-budget benches, docs.

## What is green (EXECUTED)

| Check | Result |
|---|---|
| `lake build foldkernel` | 0 · binary ~5.0 MB |
| `lake build FoldBounds.FrontierKernelSound` | 0 · includes `enumFast_eq` |
| `lake build FoldBounds.FrontierKernelTests` | 0 · `decide +kernel` suite |
| `lake build FoldBounds.FrontierKernelModelTestsHHHH` | 0 · `P=4/11` bridge |
| `lake build FoldBounds.FrontierKernelModelTestsH8Interval` | 0 · enclosure bridge |
| `lake build FoldBounds.FrontierKernelAudit` | 0 · 323+403+102+3+3 decls, standard axioms only; `enumFast_eq` listed |
| `scripts/check_kernel_gate_attack.py` | 0 · planted sorry / native_decide / look-alikes rejected |
| fixtures accept/reject via `foldkernel check` | 30/30 accept, 182/182 reject |
| Python differential (`experiments/kernel/results/differential_summary.json`) | 757 certs, 0 disagreements; 3069 mutants, 0 binary disagreements |
| Matched-budget benches | `experiments/kernel/bench/results/SUMMARY.md` (speedups vs foldcheck2 / hpenum; S1 exactness 25/25) |

## Design notes

* Executable kernel (`FrontierKernel.lean`) stays Mathlib-free; soundness
  (`FrontierKernelSound.lean`) reuses the existing frontier semantic library.
* `enumFast` (grid) is what `runTree` runs; `enumFast_eq` ties it to `enumGo`.
* Model bridge theorems that need Sound are split across modules so each
  `decide +kernel` elaborates alone (a full root-`X` decide for `n=8` while
  Sound is loaded exceeds the memory cgroup; exact `n=8` is covered by the
  compiled binary, the Mathlib-free hybrid test, and the enclosure bridge).

## Next (if continuing)

1. Keep the full `--full` bench warm; refresh SUMMARY if a longer run finishes.
2. Draft PR against `feature/fable52-certified-frontier-v1-9cd6` (or stack on #32).
3. Optional: further histogram / parser micro-opts; do not weaken soundness.
