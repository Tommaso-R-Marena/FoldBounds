# Certified frontier session status (resumable checkpoint)

Last updated: 2026-11-17 19:20 EST (checkpoint 1, reconnaissance complete).

## Git state

* Work branch: `feature/fable52-certified-frontier-v1-9cd6`
  (the session's branch-naming policy requires the `feature/…-9cd6` form; the
  prompt's suggested name `research/fable52-certified-frontier-v1` is recorded
  here for traceability).
* Base: `fc3de052d67f54e0316d56ef7577c272d4afd02d`, the head of draft PR #32
  (`integration/aristotle-frontier-v2-20261020`, "Aristotle v2"), which is
  stacked on draft PR #26. `main` is `e373a422eb785f5cb69163b0e114343b150abf08`.
* Why not `main`: the prompt's P0 target (a non-enumerative Lean checker with a
  kernel-checked acceptance theorem) already exists on draft PRs #26/#32
  (`FrontierCheck.check_sound`, `FrontierNext.check_sound`). Re-implementing it
  on `main` would duplicate that work. This branch depends on its semantic
  library (`Frontier.ext`, `sum_ext_split`, `card_ext_le_pow_three`, degree and
  geometric caps, families). PR #32's head is an ancestor of PR #91's head and
  the frontier modules are byte-identical in both.
* No other branch or PR has been modified.

## Environment (EXECUTED)

* Lean 4.28.0 (elan 4.3.0), Lake 5.0.0, Mathlib v4.28.0
  (`8f9d9cff6bd728b17a24e163c9402775d9e6a365`, `lake exe cache get`: 8007 files).
* Python 3.12.3, venv `.venv`, `pip install -e '.[test]'` (pytest 9.2.0).
* 4 cores, 15 GB RAM, Linux 6.12.

## Baseline results (EXECUTED)

`main` @ e373a42:

| Command | Exit | Result |
|---|---|---|
| `lake build` | 0 | 8041 jobs, 1m15s with the Mathlib cache |
| `lake env lean FoldBounds/AxiomAudit.lean` | 0 | 50 reports |
| `python3 scripts/check_lean_sources.py <log>` | 0 | "50 declarations, standard axioms only" |
| `lake env lean FoldBounds/ExecutableTests.lean` | 0 | |
| `lake env lean FoldBounds/LatticeModelTests.lean` | 0 | |
| `PYTHONPATH=src python3 scripts/cross_language_diff.py` | 0 | 27 groups |
| `python -m pytest -q` | 0 | 1260 passed in 26.6 s |

PR #32 head @ fc3de05 (this branch's base):

| Command | Exit | Result |
|---|---|---|
| `lake build` | 0 | 8061 jobs, 1m58s |
| `lake build foldcheck foldcheck2` | 0 | 16089 jobs, 8m04s wall, 25 CPU-min; each binary is 222 MB (Mathlib linked) |

Replication of the remaining PR #32 gates and a speed baseline of
`foldcheck`/`foldcheck2` are running (results under
`experiments/kernel/results/baseline/`).

## Decision record

1. P0 exists on drafts → do not duplicate; verify it and close its documented
   gaps instead.
2. Gaps selected (from PR #32's own reports and the PR #96 audit):
   * the compiled checkers are specification-style and slow, and link Mathlib;
   * the text parser is outside every theorem;
   * certificates are checked for one `(q, threshold)` at a time;
   * complement-mass caps ignore the threshold (`q^K` instead of `q^min(K,t-1)`).
3. Plan: a Mathlib-free executable kernel `FoldBounds/FrontierKernel.lean`
   (format `FOLDBOUNDS-KERNEL 1`, spec in `CERTIFIED_FRONTIER_KERNEL_SPEC.md`)
   with a soundness proof in `FoldBounds/FrontierKernelSound.lean` that reuses
   the existing semantic lemmas, a theorem about the string-level function, and
   an independent Python reference, producer, adversarial tests and benchmarks.

## Next concrete tasks

1. Write and compile `FoldBounds/FrontierKernel.lean` (executable, no Mathlib).
2. Prove `FoldBounds/FrontierKernelSound.lean` (state invariant, step lemma,
   enumeration, bound entries, top-level theorem); strict axiom audit.
3. Runner `foldkernel`, differential tests against the Python reference and
   against `foldcheck2`, benchmarks, report, draft PR.
